#!/usr/bin/env python
"""End-to-end FRB localization driver (pythonized batch.bash).

Runs, for one M8-era candidate (archive layout under
/dataz/dsa110/candidates/<name>/):

  1. measure   — form a full-band coherent filterbank with the bbproc
                 toolkit using PRE-correlation dedispersion (-m), measure
                 the burst time bin, and derive ntint/startoffset/
                 stopoffset for the science MS (no more hand-tuning).
  2. burstms   — voltages_to_ms.py science MS (<cand>.ms), one ntint bin
                 containing the burst.
  3. fullms    — voltages_to_ms.py full-file MS (<cand>_full.ms,
                 ntint=4096) for the burst-field image.
  4. fieldms   — create_ms.py: field MS from the ±2 h calibration-archive
                 HDF5s hard-linked by C3 (<cand_dir>/calibration/).
  5. bpcalms   — bandpass-calibrator MS for --bpcal-source (default
                 2253+161 = 3C454.3): copy the latest per-transit MS the
                 realtime calibration system writes to
                 /dataz/dsa110/operations/calibration (copied because the
                 calibrate stage modifies flags); only if none exists,
                 build one from /dataz/dsa110/operations/correlator HDF5s
                 at the latest transit. Skipped if --bpcal-ms points at an
                 existing MS.
  6. calibrate — bp_gain_cal.py: bandpass on the calibrator, apply to the
                 field and burst MSs (field.ms, frb_bp.ms).
  7. icrs      — create_all_icrs.py: MSs on <=10 RFC VLBI calibrators
                 transiting around the burst.
  8. icrsimg   — calibrate_image_icrs.py: calibrate + wsclean the RFC MSs.
  9. image     — wsclean the burst (frb) and field (fimg) images.

After this completes, run the localization notebook
(notebooks/localize.ipynb) for the astrometric tie and final position.

Derived parameters are written to <workdir>/<cand>/localization_run.json.

Run inside the casa38 environment on h23:
  /home/ubuntu/anaconda3/envs/casa38/bin/python run_localization.py <cand>
"""

import argparse
import glob
import json
import os
import re
import struct
import subprocess
import sys

import numpy as np

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
BBPROC_DIR = "/home/ubuntu/proj/dsa110-shell/dsa110-bbproc"
TOOLKIT = f"{BBPROC_DIR}/toolkit"
CORE_ANTENNAS = f"{BBPROC_DIR}/config/core_antennas.txt"
ARCHIVE_ROOT = "/dataz/dsa110/candidates"
CORRELATOR_DIR = "/dataz/dsa110/operations/correlator"
OPS_CALIBRATION_DIR = "/dataz/dsa110/operations/calibration"
WORKROOT = "/media/ubuntu/ssd/localization_processing"

NATIVE_SAMPLE_S = 32.768e-6
FULL_NTINT = 4096

STAGES = ["measure", "fullms", "burstms", "fieldms", "bpcalms",
          "calibrate", "icrs", "icrsimg", "image"]


def run(cmd, cwd=None):
    print(f"+ {' '.join(str(c) for c in cmd)}", flush=True)
    subprocess.run([str(c) for c in cmd], cwd=cwd, check=True)


def pyscript(name, *args):
    run([sys.executable, os.path.join(SCRIPTS_DIR, name), *args],
        cwd=SCRIPTS_DIR)


# --------------------------------------------------------------------
# stage 1: measure the burst bin from a pre-dedispersed filterbank
# --------------------------------------------------------------------

def read_sigproc(path):
    """Minimal SIGPROC reader -> (header dict, data [ntime, nchan])."""
    hdr = {}
    with open(path, "rb") as fh:
        def rstr():
            (n,) = struct.unpack("<i", fh.read(4))
            return fh.read(n).decode()

        assert rstr() == "HEADER_START"
        while True:
            key = rstr()
            if key == "HEADER_END":
                break
            if key in ("rawdatafile", "source_name"):
                hdr[key] = rstr()
            elif key in ("telescope_id", "machine_id", "data_type",
                         "nchans", "nbits", "nifs", "barycentric",
                         "nsamples"):
                (hdr[key],) = struct.unpack("<i", fh.read(4))
            else:
                (hdr[key],) = struct.unpack("<d", fh.read(8))
        data = np.fromfile(fh, dtype=np.float32)
    nchan = hdr["nchans"]
    data = data[: (data.size // nchan) * nchan].reshape(-1, nchan)
    return hdr, data


def measure_burst(burst, workdir, args):
    """Toolkit filterbank with -m (pre-correlation dedispersion, ref
    1530 MHz — the same shifts voltages_to_ms applies before
    correlation), so the measured peak bin maps directly onto the
    correlated time axis."""
    cdir = burst.candidate_dir
    volt_dir = cdir / "Level2" / "voltages"
    fil = os.path.join(workdir, f"{burst.name}_predm.fil")

    cmd = [TOOLKIT, "-D", volt_dir, "-E", burst.name, "-P", fil,
           "--l", repr(burst.l_rad), "--m", repr(burst.m_rad),
           "--dec-deg", repr(burst.pointing_dec_deg),
           "-m", repr(burst.dm),
           "--core", CORE_ANTENNAS,
           "--tscrunch", str(args.measure_tscrunch),
           "--mjd", repr(burst.file_start.mjd),
           "--gpu", str(args.gpu)]
    cals = sorted(glob.glob(str(cdir / "filterbank" / "cal"
                                / "beamformer_weights_sb00_*.dat")))
    if cals:
        cmd += ["-w", cals[-1]]
    else:
        print("WARNING: no cal snapshot under filterbank/cal — "
              "measuring without calibration")
    run(cmd)

    hdr, data = read_sigproc(fil)
    # per-channel bandpass normalization before collapsing, else the
    # bandpass ripple dominates the off-pulse variance
    med_c = np.median(data, axis=0)
    mad_c = np.median(np.abs(data - med_c), axis=0) * 1.4826
    good = mad_c > 0
    norm = np.zeros_like(data)
    norm[:, good] = (data[:, good] - med_c[good]) / mad_c[good]
    ts = norm.mean(axis=1)
    raw_peak = int(np.argmax(ts))
    off = np.ones(ts.size, dtype=bool)
    off[max(0, raw_peak - 100):raw_peak + 101] = False
    med = np.median(ts[off])
    mad = np.median(np.abs(ts[off] - med)) * 1.4826 or 1.0
    # boxcar width search (bursts are often scattered over many bins —
    # twmx spans ~4 ms at DM 1702)
    snr, peak, width = -np.inf, raw_peak, 1
    for w in (1, 2, 4, 8, 16, 32):
        conv = np.convolve(ts - med, np.ones(w) / np.sqrt(w),
                           "same") / mad
        b = int(np.argmax(conv))
        if conv[b] > snr:
            snr, peak, width = float(conv[b]), b, w
    tsamp = hdr.get("tsamp", args.measure_tscrunch * NATIVE_SAMPLE_S)
    t_peak_s = peak * tsamp  # from file start, at the 1530 MHz reference
    print(f"measured burst: bin {peak} (boxcar width {width} = "
          f"{width * tsamp * 1e3:.2f} ms), t = {t_peak_s:.4f} s from "
          f"file start, S/N {snr:.1f}")
    if snr < args.min_snr:
        raise RuntimeError(
            f"burst S/N {snr:.1f} < {args.min_snr} in the dedispersed "
            "filterbank — check DM/(l,m)/cal before correlating")
    return t_peak_s, snr


def derive_offsets(burst, t_peak_s, ntint):
    """Science-MS offsets: the single ntint bin containing the burst."""
    bin_f = t_peak_s / (ntint * NATIVE_SAMPLE_S)
    start = int(np.floor(bin_f))
    frac = bin_f - start
    if min(frac, 1.0 - frac) < 0.1:
        print(f"WARNING: burst sits {frac:.2f} into bin {start} — close "
              "to a bin edge; consider a different --ntint")
    nbins_total = burst.n_native_samples // ntint
    if not 0 <= start < nbins_total:
        raise RuntimeError(
            f"derived bin {start} outside file (0..{nbins_total - 1})")
    return start, start + 1


# --------------------------------------------------------------------
# stage 5: bandpass-calibrator MS from the operations correlator archive
# --------------------------------------------------------------------

def resolve_calibrator(name):
    """J2000 position from the dsacalib VLA calibrator list."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import dsacalib

    path = os.path.join(os.path.dirname(dsacalib.__file__), "data",
                        "vlacalibrators.txt")
    pat = re.compile(
        rf"^{re.escape(name)}\s+J2000\s+\S\s+(\S+)\s+(\S+)", re.M)
    with open(path) as fh:
        m = pat.search(fh.read())
    if m is None:
        raise ValueError(f"calibrator {name} not found (J2000) in {path}")
    ra, dec = m.group(1), m.group(2).replace("'", "m").replace('"', "s")
    return SkyCoord(ra, dec, unit=(u.hourangle, u.deg))


def latest_transit_files(calpos, duration, filelength):
    """Datetime strings of correlator HDF5s covering the most recent
    transit of `calpos` present in CORRELATOR_DIR."""
    from astropy.time import Time
    import astropy.units as u
    from dsacalib.routines import cal_in_datetime
    import dsacalib.constants as ct

    ovro_lon = ct.OVRO_LON * u.rad

    files = sorted(glob.glob(f"{CORRELATOR_DIR}/*_sb00.hdf5"))
    if not files:
        raise RuntimeError(f"no sb00 HDF5s under {CORRELATOR_DIR}")
    datetimes = [os.path.basename(f)[:19] for f in files]
    t_last = Time(datetimes[-1], format="isot")

    # step back one sidereal day at a time until files cover a transit
    sidereal_day = 0.9972695663 * u.day
    lst = t_last.sidereal_time("apparent", longitude=ovro_lon)
    dt = (calpos.ra - lst).wrap_at("180d").to_value(u.deg) / 360.98565
    transit = t_last + dt * u.day
    if transit > t_last:
        transit = transit - sidereal_day
    for _ in range(30):
        hits = [d for d in datetimes
                if cal_in_datetime(d, transit, duration, filelength)]
        if hits:
            return hits, hits[0][:10], transit
        transit = transit - sidereal_day
    raise RuntimeError(
        "no correlator HDF5s covering any recent transit of the "
        "bandpass calibrator")


def find_bpcal_ms(args, workdir):
    """Bandpass-calibrator MS, cheapest source first: (1) already in the
    workdir, (2) the latest per-transit MS the realtime calibration
    system writes to OPS_CALIBRATION_DIR — copied into the workdir since
    the calibrate stage rewrites its flags, (3) built from the
    operations/correlator HDF5s."""
    name = args.bpcal_source
    if not args.remake_bpcal:
        local = sorted(glob.glob(f"{workdir}/*_{name}.ms"))
        if local:
            print(f"bpcal MS exists: {local[-1]}")
            return local[-1]
    ops = sorted(glob.glob(f"{OPS_CALIBRATION_DIR}/*_{name}.ms"))
    if ops:
        src = ops[-1]
        dst = os.path.join(workdir, os.path.basename(src))
        print(f"copying realtime-cal MS {src} -> {dst}")
        run(["rsync", "-a", "--delete", src + "/", dst + "/"])
        return dst
    print(f"no *_{name}.ms under {OPS_CALIBRATION_DIR} — building from "
          "correlator HDF5s")
    return build_bpcal_ms(args, workdir)


def build_bpcal_ms(args, workdir):
    import astropy.units as u
    from dsacalib.ms_io import convert_calibrator_pass_to_ms
    from dsacalib.utils import generate_calibrator_source
    import dsacalib.config as configuration

    name = args.bpcal_source
    calpos = resolve_calibrator(name)
    print(f"{name}: {calpos.to_string('hmsdms')}")
    hits, date, transit = latest_transit_files(
        calpos, duration=args.bpcal_duration_min * u.min,
        filelength=5 * u.min)
    print(f"latest {name} transit {transit.isot} — {len(hits)} files")

    cal = generate_calibrator_source(name, ra=calpos.ra, dec=calpos.dec)
    config = configuration.Configuration()
    convert_calibrator_pass_to_ms(
        cal, date, hits, msdir=workdir, hdf5dir=CORRELATOR_DIR,
        refmjd=config.refmjd)
    made = sorted(glob.glob(f"{workdir}/*_{name}.ms"))
    if not made:
        raise RuntimeError(f"convert_calibrator_pass_to_ms produced no "
                           f"*_{name}.ms in {workdir}")
    return made[-1]


# --------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("candname", help="candidate name (M8 archive)")
    p.add_argument("--dm", type=float, default=None,
                   help="dispersion measure; default: C2 dm_median")
    p.add_argument("--ntint", type=int, default=512,
                   help="native samples per science-MS integration")
    p.add_argument("--startoffset", type=int, default=None,
                   help="override the measured science-MS start bin")
    p.add_argument("--stopoffset", type=int, default=None,
                   help="override the measured science-MS stop bin")
    p.add_argument("--bpcal-source", default="2253+161",
                   help="bandpass calibrator (dsacalib VLA list name)")
    p.add_argument("--bpcal-ms", default=None,
                   help="use this existing bandpass-calibrator MS "
                        "instead of building one")
    p.add_argument("--bpcal-flux-jy", type=float, default=15.0,
                   help="assumed calibrator flux density "
                        "(3C454.3 ~ 15 Jy at 1.4 GHz, variable)")
    p.add_argument("--bpcal-duration-min", type=float, default=10.0,
                   help="minutes of calibrator transit to convert (HDF5 "
                        "fallback only)")
    p.add_argument("--remake-bpcal", action="store_true",
                   help="refresh the workdir bpcal MS even if one exists")
    p.add_argument("--refant", default="104",
                   help="reference antenna for the CASA solves")
    p.add_argument("--measure-tscrunch", type=int, default=8,
                   help="filterbank time integration for the burst-bin "
                        "measurement (262 us at 8)")
    p.add_argument("--min-snr", type=float, default=7.0,
                   help="abort if the dedispersed-filterbank burst S/N "
                        "is below this")
    p.add_argument("--time-shift-s", type=float, default=0.0,
                   help="seconds ADDED to the archived absolute "
                        "timestamps. The C2/search labels are LATE by a "
                        "per-restart constant (wall-clock latch; measure "
                        "it from the NVSS RA offset of the full-image, "
                        "or the dump-manifest block arithmetic) — pass "
                        "the NEGATIVE of the measured lateness, e.g. "
                        "-4.90.")
    p.add_argument("--gpu", type=int, default=0, help="toolkit GPU")
    p.add_argument("--no-dispfix", action="store_true",
                   help="skip the per-channel dispersion-rotation "
                        "correction (dsavim.dispfix) on the burst/full "
                        "MSs. The pre-correlation dedispersion leaves "
                        "each channel imaged with the label-time fringe "
                        "rotation, biasing positions west by "
                        "~omega*<tau(DM)>*cos(dec) (DM-dependent; ~0.7\" "
                        "at DM 158, ~7\" at DM 1700) and smearing the "
                        "burst along RA. Default: apply.")
    p.add_argument("--stages", default=",".join(STAGES),
                   help="comma-separated stages to run "
                        f"(of: {','.join(STAGES)})")
    return p.parse_args()


def main():
    args = parse_args()
    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    bad = set(stages) - set(STAGES)
    if bad:
        sys.exit(f"unknown stage(s): {sorted(bad)}")

    from dsavim.metadata import load_burst
    burst = load_burst(args.candname, dm=args.dm,
                       time_shift_s=args.time_shift_s)
    if burst.edge_clamped:
        print("WARNING: (l, m) is clamped to the image-grid edge — the "
              "archived burst position is untrustworthy")
    workdir = os.path.join(WORKROOT, args.candname)
    os.makedirs(workdir, exist_ok=True)

    record = {"candname": args.candname, "dm": burst.dm,
              "ntint": args.ntint,
              "dispfix": not args.no_dispfix,
              "time_shift_s": args.time_shift_s,
              "pointing_dec_deg": burst.pointing_dec_deg,
              "ra_deg": burst.ra_deg, "dec_deg": burst.dec_deg,
              "file_start_mjd": burst.file_start.mjd,
              "n_native_samples": burst.n_native_samples,
              "stages": stages}

    start, stop = args.startoffset, args.stopoffset
    if "measure" in stages:
        t_peak_s, snr = measure_burst(burst, workdir, args)
        record["t_peak_measured_s"] = t_peak_s
        record["filterbank_snr"] = snr
        m_start, m_stop = derive_offsets(burst, t_peak_s, args.ntint)
        print(f"derived science offsets: {m_start}..{m_stop} "
              f"(ntint {args.ntint})")
        if start is None:
            start, stop = m_start, m_stop
    if start is None or stop is None:
        if {"burstms"} & set(stages):
            sys.exit("no science offsets: run the measure stage or pass "
                     "--startoffset/--stopoffset")
    record["startoffset"], record["stopoffset"] = start, stop

    # 2026-07-21: per-channel dispersion-rotation correction applied
    # in place right after each voltage MS is written (see
    # dsavim.dispfix; the raw MSs carry only a DATA column).
    def _dispfix(ms_path):
        if args.no_dispfix:
            return
        from dsavim.dispfix import fix_ms
        fix_ms(ms_path, burst.dm, ms_out=ms_path, column="DATA")
        print(f"dispfix (DM {burst.dm}) applied in place -> {ms_path}")

    # fullms MUST precede burstms: both write <cand>.ms, and the rename
    # to <cand>_full.ms happens here.
    if "fullms" in stages:
        nbins = burst.n_native_samples // FULL_NTINT
        pyscript("voltages_to_ms.py", args.candname,
                 "--ntint", FULL_NTINT, "--startoffset", 0,
                 "--stopoffset", nbins, "--dm", burst.dm,
                 "--time-shift-s", args.time_shift_s)
        full = os.path.join(workdir, f"{args.candname}_full.ms")
        if os.path.exists(full):
            run(["rm", "-rf", full])
        os.rename(os.path.join(workdir, f"{args.candname}.ms"), full)
        _dispfix(full)
        print(f"full-file MS -> {full}")

    if "burstms" in stages:
        pyscript("voltages_to_ms.py", args.candname,
                 "--ntint", args.ntint, "--startoffset", start,
                 "--stopoffset", stop, "--dm", burst.dm,
                 "--time-shift-s", args.time_shift_s)
        _dispfix(os.path.join(workdir, f"{args.candname}.ms"))

    if "fieldms" in stages:
        pyscript("create_ms.py", args.candname)

    bpcal_ms = args.bpcal_ms
    if "bpcalms" in stages and bpcal_ms is None:
        bpcal_ms = find_bpcal_ms(args, workdir)
    record["bpcal_ms"] = bpcal_ms

    if "calibrate" in stages:
        if bpcal_ms is None:
            sys.exit("calibrate stage needs --bpcal-ms or the bpcalms "
                     "stage")
        pyscript("bp_gain_cal.py", "--candname", args.candname,
                 "--bpcal", bpcal_ms, "--flux-jy", args.bpcal_flux_jy,
                 "--refant", args.refant)

    if "icrs" in stages:
        pyscript("create_all_icrs.py", args.candname)

    if "icrsimg" in stages:
        pyscript("calibrate_image_icrs.py", args.candname)

    if "image" in stages:
        imdir = os.path.join(workdir, "images")
        os.makedirs(imdir, exist_ok=True)
        run(["wsclean", "-weight", "briggs", "0.5", "-name", "frb",
             "-size", "10000", "10000", "-scale", "2asec",
             "-minuv-l", "200", "-niter", "20",
             "-auto-threshold", "4.0", "-theoretic-beam",
             "../frb_bp.ms"], cwd=imdir)
        run(["wsclean", "-weight", "briggs", "0.5", "-name", "fimg",
             "-size", "6000", "6000", "-scale", "2asec",
             "-minuv-l", "200", "-niter", "1000",
             "-auto-threshold", "4.0", "-theoretic-beam",
             "../field.ms"], cwd=imdir)

    with open(os.path.join(workdir, "localization_run.json"), "w") as fh:
        json.dump(record, fh, indent=2)
    print(f"done — parameters in {workdir}/localization_run.json")
    print("next: run notebooks/localize.ipynb")


if __name__ == "__main__":
    main()
