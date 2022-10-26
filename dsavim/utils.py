"""Simple utilities for T3 imaging.
"""
import subprocess
import json, os, glob
import re
import yaml
from astropy.time import Time
from dsautils import cnf
import dsamfs.utils as pu
from antpos import utils
import dsacalib.constants as ct
import astropy.units as u
import numpy as np
from pkg_resources import resource_filename
from astropy.coordinates import SkyCoord

RFCFILE = resource_filename('dsavim', 'data/rfc_2022b_cat.txt')

def rfc_lookup(pos: 'astropy.coordinates.SkyCoord', dra: float = 15., ddec: float = 1.5, thresh: float = 0.04, findNearest: bool = False) -> dict:
    """look up nearby RFC calibrators"""

    cod,nam1,jname,rhr,rmn,rsec,ddeg,dmn,dsec,D_alp,D_Del,Corr,Obs,Sband,sflux,Cband,sflux,Xband,xflux,Uband,uflux,Kband,kflux,Type,Cat = np.genfromtxt(RFCFILE,dtype=str).transpose()

    ra = []
    dec = []
    for i in np.arange(len(rhr)):
        ra.append(f"{rhr[i]}:{rmn[i]}:{rsec[i]}")
        dec.append(f"{ddeg[i]}:{dmn[i]}:{dsec[i]}")
    catpos = SkyCoord(ra,dec,unit=(u.hourangle, u.deg),frame='icrs')
    flux = []
    for i in np.arange(len(rhr)):
        flux.append(np.max(np.asarray([float(Sband[i]),float(Cband[i]),float(Xband[i])])))
    flux = np.asarray(flux)

    # cross-match
    ppos = SkyCoord([pos.ra],[pos.dec])
    idx, d2d, d3d = catpos.match_to_catalog_sky(ppos)

    # limit by max sep
    wrs = np.asarray((np.where(d2d.deg<np.sqrt(dra*dra+ddec*ddec)))[0]).astype('int')
    catposes = catpos[wrs]
    fluxes = flux[wrs]
    jnames = jname[wrs]
    
    # now run through and make / print strip matches
    oposes = []
    ofluxes = []
    ojnames = []
    for i,poses in enumerate(catposes):
        delra, deldec = poses.spherical_offsets_to(pos)
        if np.abs(delra.deg)<dra:
            if np.abs(deldec.deg)<ddec:
                if fluxes[i]>thresh:
                    oposes.append(poses); ofluxes.append(fluxes[i]); ojnames.append(jnames[i])
                    #print(f"{jnames[i]}: {poses.to_string('hmsdms')}, flux (Jy) = {fluxes[i]}, delRA = {delra}, delDEC = {deldec}")

    seps = []
    for i,poses in enumerate(oposes):
        sep = poses.separation(pos)
        seps.append(sep.deg)

    # do the findNearest if needed
    if findNearest:
        ooposes = []
        oofluxes = []
        oojnames = []
        oseps = []

        iseps = np.argsort(np.abs(np.asarray(seps)))
        found = 0
        ii = 0
        while found != 2:
            i = iseps[ii]
            ii += 1
            delra, deldec = oposes[i].spherical_offsets_to(pos)
            if found==0:
                ooposes.append(oposes[i])
                oofluxes.append(ofluxes[i])
                oojnames.append(ojnames[i])
                oseps.append(seps[i])
                if delra.deg<0.:
                    found=-1
                else:
                    found=1
            if found==1:
                if delra.deg<0.:
                    ooposes.append(oposes[i])
                    oofluxes.append(ofluxes[i])
                    oojnames.append(ojnames[i])
                    oseps.append(seps[i])
                    found=2
            if found==-1:
                if delra.deg>0.:
                    ooposes.append(oposes[i])
                    oofluxes.append(ofluxes[i])
                    oojnames.append(ojnames[i])
                    oseps.append(seps[i])
                    found=2
                
        return {"jname":oojnames, "position":ooposes, "flux":oofluxes, "sep":oseps}

    return {"jname":ojnames, "position":oposes, "flux":ofluxes, "sep":seps}
    

def load_params(paramfile: str) -> dict:
    """Load parameters for voltage correlation from a yaml file."""
    with open(paramfile) as yamlf:
        voltage_corr_params = yaml.load(yamlf, Loader=yaml.FullLoader)['voltage_corr']
    conf = cnf.Conf(use_etcd=False)
    corrconf = conf.get('corr')
    mfsconf = conf.get('fringe')

    voltage_corr_params['ch0'] = corrconf['ch0']
    voltage_corr_params['f0_GHz'] = corrconf['f0_GHz']
    voltage_corr_params['antennas'] = list(corrconf['antenna_order'].values())[:63]
    ant_od = corrconf['antenna_order']
    antenna_order = [int(ad) for ad in list(ant_od.values())]
    voltage_corr_params['outrigger_delays'] = mfsconf['outrigger_delays']
    df = utils.get_itrf(
        latlon_center=(ct.OVRO_LAT * u.rad, ct.OVRO_LON *
                       u.rad, ct.OVRO_ALT * u.m)
    )
    ant_itrf = np.array([df['dx_m'], df['dy_m'], df['dz_m']]).T
    nants_telescope = max(df.index)
    voltage_corr_params['ant_itrf'] = ant_itrf
    voltage_corr_params['nants_telescope'] = nants_telescope
    voltage_corr_params['snapdelays'] = pu.get_delays(antenna_order, nants_telescope)    

    return voltage_corr_params


def get_tstart_from_json(headername: str) -> 'astropy.time.Time':
    """Extract the start time from the header file."""
    with open(headername) as jsonf:
        metadata = json.load(jsonf)
    metadata = metadata[list(metadata.keys())[0]]
    tstart = Time(metadata['mjds'], format='mjd')
    return tstart


def get_DM_from_json(headername: str) -> float:
    """Extract the dispersion measure from the header file."""
    with open(headername) as jsonf:
        metadata = json.load(jsonf)
    metadata = metadata[list(metadata.keys())[0]]
    return metadata['dm']


def find_beamformer_weights(candtime: 'astropy.time.Time', bfdir: str = '/data/dsa110/T3/calibs/') -> str:
    """Find the beamformer weights that were in use at a time `candtime`.
    
    In `/data/dsa110/T3/calibs/`, the times in the beamformer weight names are the times when they were
    uploaded to the correlator nodes. Therefore, we want the most recent calibration files that were
    created before `candtime`.
    """
    isot_string = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]"
    isot_pattern = re.compile(isot_string)
    avail_calibs = sorted(
        [
            isot_pattern.findall(calib_path)[0] for calib_path
            in glob.iglob(f"{bfdir}/beamformer_weights_{isot_string}.yaml")],
        reverse=True)
    for avail_calib in avail_calibs:
        if avail_calib < isot_pattern.findall(candtime.isot)[0]:
            return avail_calib


def rsync_file(infile: str, outfile: str) -> str:
    """Rsyncs a file from the correlator machines.

    Parameters
    ----------
    infile : str
        The sourcefile string, e.g. 'corr01.sas.pvt:/home/user/data/fl_out.1.5618974'
    outfile : str
        The destination string, e.g. '/home/user/data/fl_out.1.5618974'

    Returns
    -------
    str
        The full path to the rsynced file in its destination.
    """
    command = f"rsync -avvP --inplace {infile} {outfile}"
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True
    )
    proc_stdout = str(process.communicate()[0].strip())
    print(proc_stdout)
    if process.returncode != 0:
        return None
    return outfile
