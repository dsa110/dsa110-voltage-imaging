"""Burst metadata from the dsa110-rt (M8-era) candidate archive.

The 2026+ realtime system writes, per candidate under
``/dataz/dsa110/candidates/<name>/``:

* ``Level3/<name>.json`` — C2 cluster record: ``c2.t_peak_mjd`` (pulse
  arrival at the BOTTOM of the processed band, 1.311 GHz),
  ``c2.dm_median``, ``c2.snr_max``, ``c2.l_median``/``c2.m_median``
  (image-plane direction cosines from the dec-fringe-stopped meridian,
  radians), ``c2.peak_event_specnum``;
* ``Level2/voltages/<name>_sbNN.json`` — per-subband voltage-dump
  manifests: ``mjd_target``, ``target_block_n``, ``block_n_first/last``,
  ``n_blocks_written``, ``dropped_block_ns`` (block = 0.134217728 s =
  4096 native samples of 32.768 us);
* ``filterbank/filterbank.json`` — includes ``dec_deg`` (the pointing
  declination used for the F21 fold) when C3 generated filterbanks.

This module rolls those into one :class:`Burst` object that the
localization pipeline consumes — replacing the T3-era
``T2_<name>.json`` (``mjds``/``dm``/``ra``/``dec``/``specnum``) reader.

Conventions (validated 2026-07-15):

* ``RA = LST_ovro(t_peak) + l/cos(dec)`` — the +l→east sign pinned
  against the edge-clamped B1933+16 pulse 260715qwfe;
* ``Dec = dec_pointing + m`` — F21/injection convention (m multiplies
  the antenna northings);
* positions here are PHASE-CENTER quality (~arcmin; LST via GMST
  approximation). The final FRB position always comes from imaging +
  the RFC astrometric tie, never from this.
* NOTE the image grid spans ~±1 deg: ``l_median``/``m_median`` values
  at exactly ±0.017453 rad (±1.0 deg) are edge-clamped and NOT
  trustworthy (260715qwfe).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from astropy.time import Time

ARCHIVE_ROOT = Path("/dataz/dsa110/candidates")

#: fada block length (s) and native sample (s) — dsa110-rt constants.
BLOCK_S = 0.134217728
NATIVE_SAMPLE_S = 32.768e-6
NATIVE_PER_BLOCK = 4096

OVRO_LON_DEG = -118.2817

SUBBANDS = [f"sb{i:02d}" for i in range(16)]


def _lst_deg(mjd: float) -> float:
    """Approximate local apparent sidereal time at OVRO (deg).

    GMST polynomial good to ~1 s of time (~15 arcsec) — phase-center
    quality, which is all this module promises."""
    d = mjd - 51544.5
    gmst_h = (18.697374558 + 24.06570982441908 * d) % 24.0
    return (gmst_h * 15.0 + OVRO_LON_DEG) % 360.0


@dataclass(frozen=True)
class Burst:
    """Everything the localization pipeline needs about one candidate."""

    name: str
    candidate_dir: Path
    # burst
    t_peak: Time              # pulse arrival at the band bottom (1.311 GHz)
    dm: float                 # pc/cc (C2 dm_median unless overridden)
    snr: float
    l_rad: float
    m_rad: float
    pointing_dec_deg: float
    ra_deg: float             # phase-center quality (see module docstring)
    dec_deg: float
    # voltage files
    voltage_files: List[Path]
    file_start: Time          # first sample of every fragment
    n_native_samples: int     # per fragment (n_blocks_written * 4096)
    burst_offset_s: float     # t_peak - file_start, seconds
    dropped_blocks: dict      # sb -> list of dropped block_ns

    @property
    def calibration_dir(self) -> Path:
        return self.candidate_dir / "calibration"

    @property
    def edge_clamped(self) -> bool:
        """True if the (l, m) sits on the ±1 deg image-grid edge."""
        lim = math.radians(1.0)
        return (abs(abs(self.l_rad) - lim) < 1e-6
                or abs(abs(self.m_rad) - lim) < 1e-6)


def load_burst(candname: str,
               archive_root: Path = ARCHIVE_ROOT,
               dm: Optional[float] = None,
               pointing_dec_deg: Optional[float] = None,
               time_shift_s: float = 0.0) -> Burst:
    """Read the M8 candidate archive into a :class:`Burst`.

    ``dm``/``pointing_dec_deg`` override the archived values (the DM
    used for pre-correlation dedispersion may legitimately be refined
    beyond the search's dm_median).

    ``time_shift_s`` is ADDED to the LABEL-derived timestamps (t_peak,
    and file_start when falling back to mjd_target anchoring). Since
    2026-07-19 ``file_start`` is anchored on the manifests'
    ``block_mjd_first`` (capture-arm chain: armed_mjd + block
    arithmetic — never had the wall-clock-latch bug), so for M8
    manifests with armed provenance ``time_shift_s`` no longer moves
    the voltage-MS time axis. It remains the workaround for older
    archives, where the C2/search wall-clock latch (search_compute
    mjd_at_specnum_0) made every label LATE by a per-restart constant
    (measured +4.9 s for the 2026-07-15 morning epoch, +6.3 s for the
    14:18 epoch). Pass the NEGATIVE of the measured lateness there.
    """
    cdir = Path(archive_root) / candname
    l3path = cdir / "Level3" / f"{candname}.json"
    with l3path.open() as fh:
        c2 = (json.load(fh) or {}).get("c2") or {}
    if not c2:
        raise ValueError(f"no c2 block in {l3path}")

    t_peak_mjd = float(c2["t_peak_mjd"]) + time_shift_s / 86400.0
    l_rad = float(c2.get("l_median") or 0.0)
    m_rad = float(c2.get("m_median") or 0.0)
    if dm is None:
        dm = float(c2["dm_median"])

    # pointing dec: filterbank provenance (exact value C3 used), else
    # explicit argument.
    if pointing_dec_deg is None:
        fbmeta = cdir / "filterbank" / "filterbank.json"
        if fbmeta.is_file():
            with fbmeta.open() as fh:
                dd = (json.load(fh) or {}).get("dec_deg")
            if dd is not None:
                pointing_dec_deg = float(dd)
    if pointing_dec_deg is None:
        raise ValueError(
            f"{candname}: pointing declination unavailable (no "
            "filterbank/filterbank.json dec_deg) — pass pointing_dec_deg")

    # file timing from the first available voltage manifest
    vdir = cdir / "Level2" / "voltages"
    manifest = None
    dropped = {}
    for sb in SUBBANDS:
        mpath = vdir / f"{candname}_{sb}.json"
        if mpath.is_file():
            with mpath.open() as fh:
                doc = json.load(fh)
            dropped[sb] = list(doc.get("dropped_block_ns") or [])
            if manifest is None:
                manifest = doc
    file_start_from_arm: Optional[float] = None
    if manifest is not None:
        n_blocks = int(manifest["n_blocks_written"])
        pre_blocks = int(manifest["target_block_n"]) - int(
            manifest["block_n_first"])
        # 2026-07-19 (B1933 ladder analysis): anchor the file start on
        # ``block_mjd_first`` — the corr-side capture-arm chain
        # (armed_mjd + block arithmetic, same base as the slow-vis
        # HDF5s) — NOT on the C2 label. The label (mjd_target =
        # c2.t_peak_mjd) carries per-event cluster-centroid jitter of
        # ~0.1 s (~1.5" of RA in the voltage MS) even after the
        # 18b0d0c wall-latch fix. ``block_mjd_first`` records the
        # correct time for its file's first byte under BOTH manifest
        # versions (v1's 1-based labels and v2's 0-based labels use
        # matching formulas), so no version switch is needed here.
        # The armed chain never had the wall-clock-latch bug, so
        # ``time_shift_s`` deliberately does NOT apply to it.
        bmf = manifest.get("block_mjd_first")
        if bmf is not None and float(bmf) > 40000.0:
            file_start_from_arm = float(bmf)
        mjd_target = float(manifest.get("mjd_target") or 0.0)
        if mjd_target < 40000.0:
            # pre-fix manifests recorded 0.0; use the C2 record
            # (t_peak_mjd already carries time_shift_s)
            mjd_target = t_peak_mjd
        else:
            mjd_target += time_shift_s / 86400.0
    else:
        # Manifest-less fallback (events before 2026-07-15, when C3
        # started collecting the staging manifests): block count from
        # the fragment size; the pre-window from the era of the dump
        # (n_pre was 8 before the 2026-07-15 ~14:18 UTC corr restart
        # that deployed the sweep-aware 14/8 window).
        frag = next((f for f in
                     (vdir / f"{candname}_{sb}_data.out" for sb in SUBBANDS)
                     if f.is_file()), None)
        if frag is None:
            raise FileNotFoundError(
                f"{candname}: no voltage fragments under {vdir}")
        block_bytes = 301989888
        n_blocks = int(frag.stat().st_size // block_bytes)
        pre_blocks = 8 if t_peak_mjd < 61236.59861 else 14
        mjd_target = t_peak_mjd
        import warnings
        warnings.warn(
            f"{candname}: no staging manifests; reconstructed timing "
            f"(n_blocks={n_blocks}, pre_blocks={pre_blocks} by era)")
    if file_start_from_arm is not None:
        file_start_mjd = file_start_from_arm
    else:
        # label-anchored fallback (manifest-less events and manifests
        # predating armed_mjd provenance, 2026-07-16)
        file_start_mjd = mjd_target - pre_blocks * BLOCK_S / 86400.0

    # sky position (phase-center quality)
    dec_deg = pointing_dec_deg + math.degrees(m_rad)
    ra_deg = (_lst_deg(t_peak_mjd)
              + math.degrees(l_rad) / math.cos(math.radians(dec_deg))
              ) % 360.0

    voltage_files = [vdir / f"{candname}_{sb}_data.out" for sb in SUBBANDS]

    return Burst(
        name=candname,
        candidate_dir=cdir,
        t_peak=Time(t_peak_mjd, format="mjd"),
        dm=dm,
        snr=float(c2.get("snr_max") or 0.0),
        l_rad=l_rad,
        m_rad=m_rad,
        pointing_dec_deg=pointing_dec_deg,
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        voltage_files=voltage_files,
        file_start=Time(file_start_mjd, format="mjd"),
        n_native_samples=n_blocks * NATIVE_PER_BLOCK,
        burst_offset_s=(t_peak_mjd - file_start_mjd) * 86400.0,
        dropped_blocks=dropped,
    )


def burst_bin(burst: Burst, ntint: int,
              reffreq_ghz: float = 1.53) -> float:
    """Post-correlation time-bin index of the (dedispersed) burst.

    With pre-correlation dedispersion referenced to ``reffreq_ghz``
    (the legacy toolkit convention, 1530 MHz), every channel's burst
    power lands at the arrival time extrapolated to that frequency:

        t_ref = t_peak − τ(reffreq → 1.311 GHz, DM)

    measured from the file start, divided by the integration time.
    """
    nu_bot = 1.311281  # GHz, processed band bottom (t_peak reference)
    tau = 4.15e-3 * burst.dm * (nu_bot ** -2 - reffreq_ghz ** -2)
    t_ref = burst.burst_offset_s - tau
    return t_ref / (ntint * NATIVE_SAMPLE_S)
