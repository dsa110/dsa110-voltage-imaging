from collections import namedtuple
import astropy.units as u
from pkg_resources import resource_filename

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
T3PARAMS = load_params(PARAMFILE)
BURST_START_S = 1907*262.144e-6

def pipeline_component(targetfn, inqueue, outqueue=None):
    """Generate a component of the pipeline."""
    def inner():
        done = False
        while not done:
            try:
                item = inqueue.get()
            except queue.Empty:
                time.sleep(10)
                continue
            if item == 'END':
                done = True
            else:
                item = targetfn(item)
            if outqueue is not None:
                outqueue.put(item)
    return inner

def generate_rsync_component(rsync: bool) -> "Callable":
    """Generate an rsync function."""

    def rsync_all_files(item):
        """Rsync or copy a file."""
        srcfile, vfile = item
        rsync_done = False
        if not os.path.exists(vfile):
            if rsync:
                rsync_file(srcfile, vfile)
            else:
                os.symlink(srcfile, vfile)
        return vfile
    return rsync_all_files

def generate_correlate_component(
        ntint: int, corr_ch0: dict, full_pol: bool, ncorrfiles: "Manager().Value",
        ncorrfiles_lock: "Manager().Lock") -> "Callable":
    """Generate a correlator function."""

    def correlate(vfile):
        """Correlate a file."""
        if ncorrfiles.value > 2:
            time.sleep(10)

        corr = re.findall('corr\d\d', vfile)[0]
        if not os.path.exists('{0}.corr'.format(vfile)):
            first_channel_MHz = corr_ch0[corr]
            command = (
                '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/toolkit_dev '
                f'-i {vfile} -o {vfile}.corr -t {ntint} -c {first_channel_MHz} '
                f'-d delays.dat {"" if full_pol else "-a"}')
            print(command)
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True)
            proc_stdout = str(process.communicate()[0].strip())
            print(proc_stdout)

        corrfile = '{0}.corr'.format(vfile)

        with ncorrfiles_lock:
            ncorrfiles.value += 1

        return corrfile

    return correlate

def generate_uvh5_component(
        candname: str, declination: "Manager().Value", vis_params: dict, start_offset: int,
        end_offset: int, ncorrfiles: "Manager().Value", ncorrfiles_lock: "Manager().Lock") -> "Callable":
    """Generate a uvh5 writer."""
    def write_uvh5(corrfile):
        uvh5name = generate_uvh5(
            '{0}/{1}'.format(T3PARAMS['corrdir'], candname),
            declination.value*u.deg,
            corrfile=corrfile,
            vis_params=vis_params,
            start_offset=start_offset,
            end_offset=end_offset
        )

        os.remove(corrfile)
        with ncorrfiles_lock:
            ncorrfiles.value -= 1

    return write_uvh5

def process_join(targetfn):
    def inner():
        process = Process(
            target=targetfn,
            daemon=True)
        process.start()
        process.join()
        return process
    return inner

def generate_declination_component(
        declination: "Manager().Value", declination_lock: "Manager().Lock",
        tstart: "astropy.time.Time") -> "Callable":
    def get_declination_etcd():
        with declination_lock:
            if declination.value is None:
                declination.value = get_declination(
                    get_elevation(tstart)
                ).to_value(u.deg)
    return get_declination_etcd

def generate_delay_table(headername, declination):
    # TODO: Do this with python, instead of with Vikram's script
    gen_delay_script = ('/home/ubuntu/anaconda3/envs/dana/bin/python '
                        '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/gen_delays.py')
    command = f'{gen_delay_script} {headername} {declination}'
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        shell=True)
    proc_stdout = str(process.communicate()[0].strip())
    print(proc_stdout)

def initialize_system():
    SystemSetup = namedtuple(
        "SystemSetup", "T3dir T3archivedir corrdir msdir start_time_offset ref_freq corr_ch0_MHz")
    start_time_offset = BURST_START_S*u.s
    corr_ch0_MHz = {key: 1e3*vis_params['fobs'][value] for key, value in vis_params['corr_ch0'].items()}

    system_setup = SystemSetup()
    return system_setup

def initialize_candidate(candname, datestring, system_setup):
    Candidate = namedtuple('Candidate', 'name time voltagefiles local')
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    t3dir = system_setup.T3dir
    archivedir = system_setup.archivedir

    if datestring == 'current':
        local = False
        headerfile = f'{system_setup.t3dir}/{candname}.json'
        voltagefiles = [f'{corr}.sas.pvt:/home/ubuntu/data/{candname}_data.out' for corr in corrlist]

    else:
        local = True
        headerfile = f'{system_setup.archivedir}/{datestring}/{candname}.json'
        voltagefiles = [f'{system_setup.archivedir}/{datestring}/{corr}_{candname}_data.out'
                        for corr in corrlist]

    tstart = get_tstart_from_json(headerfile)
    dispersion_measure = get_DM_from_json(headername)
    if dispersion_measure < 1:
        dispersion_measure = None

    cand = Candidate(candname, start_time, voltagefiles, local)

def initialize_correlater(fullpol, ntint, cand, system_setup):
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    reftime = cand.time + system_setup.start_time_offset
    npol = 4 if full_pol else 2
    nfint = 1 if full_pol else 8
    corrfiles = [f'{system_setup.corrdir}/{corr}_{cand.name}_data.out' for corr in CORR_LIST]

    CorrelaterParameters = namedtuple('Correlator', 'reftime npol nfint ntint files')
    correlater_params = CorrelaterParameters(reftime, npol, nfint, ntint, corrfiles)

def initialize_uvh5(cand, system_setup):
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    UVH5Parameters = namedtuple('UVH5', 'files')
    uvh5files = [f'{system_setup.corrdir}/{cand.name}_{corr}.hdf5' for corr in corrlist]
    uvh5_params = UVH5Parameters(uvh5files)

def initialize_vis_params(corrparams, cand, system_setup)
    vis_params = parse_visibility_parameters(T3PARAMS, cand.time, corrparams.ntint)
    vis_params['tref'] = corrparams.reftime
    vis_params['npol'] = corrparams.npol
    vis_params['nfint'] = corrparams.nfint
    vis_params['ntint'] = corrparams.ntint
    return vis_params


    corr_ch0_MHz_safe = MappingProxyType(corr_ch0_MHz) # Thread-safe mapping
    vis_params_safe = MappingProxyType(vis_params)

# Candidate parameters
# name time dm pointing_dec_deg voltagefiles local

# System parameters
# t3dir archivedir corrdir msdir start_time_offset ref_freq corr_ch0_MHz

# Correlater parameters
# reftime npol nfint ntint corrfiles

# UVH5 parameters
# ??? uvh5files

# MS parameters
# msname