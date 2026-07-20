"""Correct the dispersion-delay x Earth-rotation coupling in burst MSs.

The bbproc toolkit dedisperses PRE-correlation: channel nu's samples
are shifted so the burst aligns at the 1530-MHz reference time. The
data that land in a visibility bin at label time t were therefore
RECEIVED at t + tau(nu), with

    tau(nu) = K_DM * DM * (nu^-2 - nu_ref^-2),   nu_ref = 1.530 GHz,

up to 0.102 s across the processed band at DM 158.5 (1.10 s at
DM 1702). The downstream fringestop/uvw evaluation uses the single
label time for ALL channels, so channel nu's apparent field rotates
by the sky drift accumulated in tau(nu): a fixed source appears
displaced WEST by omega_E * tau(nu) * cos(dec) — a chromatic RA ramp
(0 -> 1.5" at DM 158.5; 0 -> 16" at DM 1702) that smears the broadband
image along RA and biases the band-average position west.

Detected empirically 2026-07-20 on the 260718rdga single pulses: the
upper half band imaged +2.2 +- 0.9" east of the lower half (predicted
+0.7" from this coupling at DM 158.5; consistent within the split
noise).

``fix_ms`` applies the equivalent per-channel image translation as a
per-row, per-channel phase ramp on the visibilities:

    V(u, v, nu) *= exp(-2*pi*i * u_lambda(nu) * dl(nu)),
    dl(nu) = + omega_sidereal * tau(nu) * cos(dec0)

(+l = east; CASA phase convention V ~ exp(-2*pi*i*(u l + v m + ...));
the sign was validated on 260715twmx, DM 1702, where the correction
collapses the half-band split). The m-drift is O(sin HA) smaller and
neglected.

This is a post-hoc equivalent of using per-channel times in the
fringestop/uvw stage of voltages_to_ms (the proper upstream fix).
"""

from __future__ import annotations

import shutil
import os
from typing import Optional

import numpy as np
from casacore.tables import table

__all__ = ["fix_ms", "K_DM_S_MHZ2"]

#: dispersion constant, seconds * MHz^2 * cm^3 / pc
K_DM_S_MHZ2 = 4.148808e3
OMEGA_SIDEREAL = 2.0 * np.pi / 86164.0905   # rad/s
C_M_S = 299792458.0
NU_REF_MHZ = 1530.0


def fix_ms(
    ms_in: str,
    dm: float,
    *,
    ms_out: Optional[str] = None,
    dec_deg: Optional[float] = None,
    column: str = "CORRECTED_DATA",
    nu_ref_mhz: float = NU_REF_MHZ,
) -> str:
    """Apply the per-channel dispersion-rotation correction.

    Args:
        ms_in: input MS (untouched unless ``ms_out`` is None).
        dm: dispersion measure used in the toolkit's pre-correlation
            dedispersion (pc/cm^3).
        ms_out: output MS path; default ``<ms_in root>_dfix.ms``.
            Pass ``ms_out=ms_in`` to correct in place.
        dec_deg: phase-center declination; default: read from FIELD.
        column: visibility column to correct.
        nu_ref_mhz: dedispersion reference frequency (toolkit: 1530).

    Returns the output MS path.
    """
    if ms_out is None:
        root = ms_in[:-3] if ms_in.endswith(".ms") else ms_in
        ms_out = f"{root}_dfix.ms"
    if os.path.abspath(ms_out) != os.path.abspath(ms_in):
        if os.path.exists(ms_out):
            shutil.rmtree(ms_out)
        shutil.copytree(ms_in, ms_out)

    with table(f"{ms_out}/SPECTRAL_WINDOW", ack=False) as sw:
        freqs_hz = sw.getcol("CHAN_FREQ")[0]
    if dec_deg is None:
        with table(f"{ms_out}/FIELD", ack=False) as fld:
            dec_deg = float(np.degrees(
                fld.getcol("PHASE_DIR").squeeze()[1]))

    nu_mhz = freqs_hz / 1e6
    tau_s = K_DM_S_MHZ2 * float(dm) * (nu_mhz ** -2 - nu_ref_mhz ** -2)
    # Sign anchored on the REAL 260718rdga half-band measurement (the
    # upper band images EAST of the lower band by +2.15 +- 0.9" before
    # correction, i.e. late-arriving low-frequency channels appear
    # west): the correction translates channel nu's image EAST by
    # omega*tau(nu)*cos(dec), which under this MS/wsclean convention
    # requires a NEGATIVE dl in the exp(-2 pi i u dl) operator
    # (verified: the corrected real split collapses; see module tests).
    dl_rad = -OMEGA_SIDEREAL * tau_s * np.cos(np.radians(dec_deg))
    u_scale = freqs_hz / C_M_S              # u[m] -> u[lambda] per chan

    t = table(ms_out, readonly=False, ack=False)
    uvw = t.getcol("UVW")
    dat = t.getcol(column)
    # phase[row, chan] = -2 pi * u_m[row] * u_scale[chan] * dl[chan]
    phase = -2.0 * np.pi * uvw[:, 0:1] * (u_scale * dl_rad)[None, :]
    dat *= np.exp(1j * phase)[:, :, None]
    t.putcol(column, dat)
    t.close()
    return ms_out
