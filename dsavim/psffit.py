"""Snapshot burst position estimator (constrained peak + 2-D quadratic).

Motivation and validation (2026-07-20, B1933+16 per-pulse astrometry
study on the 260718 SNR-ladder voltage dumps; 32 single-pulse images):

* The DSA-110 single-bin snapshot beam is a ~13-28" core on a ~68"
  half-power plateau (PA ~105 deg east of north). Single-pulse position
  errors are thermal with sigma_RA ~= 40"/SNR, sigma_Dec ~= 15"/SNR
  (the plateau's curvature scale, not the core, sets the RA figure).

* For a point source in white visibility noise under natural(-ish)
  weighting, the dirty image IS the matched filter and its noise
  covariance is proportional to the PSF, so the whitened template fit
  degenerates to the CONSTRAINED ARGMAX OF THE IMAGE ITSELF. Fitting
  or correlating the PSF template against the image applies the beam
  twice and measurably WORSENS the scatter (validated: 2.8" -> 3.3-4.9"
  per event) — do not "PSF-fit" these images.

* The Cramer-Rao bound (~19"/SNR in RA) is NOT reachable for single
  pulses at SNR 10-30: the gap is the estimation-theory threshold
  effect on the plateau's ambiguity structure, not estimator
  inefficiency. Weighting changes don't help either (briggs 0/-1,
  uniform, and u-only density flattening all measured worse than
  briggs 0.5 — the S/N cost always beats the curvature gain).

* What does help, slightly: replacing the 1-D 3-point parabola with a
  2-D quadratic least-squares refinement over a +-3 px neighbourhood
  of the constrained argmax: mean RA scatter 2.79" -> 2.58" (-8%),
  Dec 1.03" -> 0.89" (-13%) across the 32-pulse validation set, never
  worse than the parabola on any event. That is this module. The
  validated error model for its output (ML fit over the 32 pulses;
  note beta and sys are partially degenerate at this sample size):

      sigma_RA  ~ hypot(20"/SNR, 1.8")
      sigma_Dec ~ hypot(12"/SNR, 0.25")     (briggs 0.5 images)

Usage:
    from dsavim.psffit import fit_position
    res = fit_position("frb-image.fits", expected_coord)
    res.coord, res.snr, res.dra_arcsec, res.ddec_arcsec
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Union

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

__all__ = ["PsfFitResult", "fit_position"]


@dataclass
class PsfFitResult:
    coord: SkyCoord           #: fitted position
    dra_arcsec: float         #: (fit - expected) * cos(dec), arcsec
    ddec_arcsec: float        #: (fit - expected), arcsec
    snr: float                #: peak / clipped off-source rms
    peak: float               #: image value at the integer peak
    curvature_ok: bool        #: 2-D quadratic had a proper maximum
    sigma_ra_est: float       #: 37"/SNR error estimate (arcsec)
    sigma_dec_est: float      #: 13"/SNR error estimate (arcsec)


def _clipped_rms(a: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    cut = np.percentile(np.abs(a), 95)
    return float(np.std(a[np.abs(a) < cut]))


def fit_position(
    image: Union[str, "np.ndarray"],
    expected: SkyCoord,
    *,
    wcs: Optional[WCS] = None,
    search_asec: float = 25.0,
    quad_halfwidth: int = 3,
) -> PsfFitResult:
    """Constrained-argmax + 2-D quadratic position fit.

    Args:
        image: wsclean image FITS path (restored or dirty), or a 2-D
            array with ``wcs`` supplied.
        expected: prior position (the C2/C1 detection position for
            FRBs). The argmax is constrained to ``search_asec`` of it —
            essential on the plateau-dominated snapshot beam, where an
            unconstrained peak can lock onto a distant noise bump.
        search_asec: half-width of the prior box. Match it to the
            upstream position uncertainty (~half a detection-grid cell
            for C2 positions).
        quad_halfwidth: half-width (pixels) of the 2-D quadratic
            refinement neighbourhood. 3 validated best at the 2"/px
            production scale.
    """
    if isinstance(image, str):
        hdu = fits.open(image)[0]
        img = hdu.data.squeeze().astype(np.float64)
        wcs = WCS(hdu.header).celestial
    else:
        img = np.asarray(image, dtype=np.float64)
        if wcs is None:
            raise ValueError("wcs required when passing arrays")

    scale = abs(wcs.wcs.cdelt[1]) * 3600.0
    x, y = wcs.world_to_pixel(expected)
    xi, yi = int(round(float(x))), int(round(float(y)))
    s = int(round(search_asec / scale))
    q = int(quad_halfwidth)
    if not (s + q < xi < img.shape[1] - s - q
            and s + q < yi < img.shape[0] - s - q):
        raise ValueError("expected position too close to the image edge")

    cut = img[yi - s:yi + s + 1, xi - s:xi + s + 1]
    iy, ix = np.unravel_index(np.nanargmax(cut), cut.shape)
    cy, cx = yi - s + iy, xi - s + ix

    yy, xx = np.mgrid[-q:q + 1, -q:q + 1]
    patch = img[cy - q:cy + q + 1, cx - q:cx + q + 1]
    A = np.column_stack([np.ones(xx.size), xx.ravel(), yy.ravel(),
                         xx.ravel() ** 2, yy.ravel() ** 2,
                         (xx * yy).ravel()])
    coef, *_ = np.linalg.lstsq(A, patch.ravel(), rcond=None)
    H = np.array([[2 * coef[3], coef[5]], [coef[5], 2 * coef[4]]])
    curv_ok = bool(np.all(np.linalg.eigvalsh(H) < 0))
    if curv_ok:
        dx, dy = np.linalg.solve(H, [-coef[1], -coef[2]])
        dx = float(np.clip(dx, -q, q))
        dy = float(np.clip(dy, -q, q))
    else:
        dx = dy = 0.0

    coord = wcs.pixel_to_world(cx + dx, cy + dy)

    # off-source rms from bands well away from the search box
    R = max(4 * s, 50)
    y0, y1 = max(0, yi - R), min(img.shape[0], yi + R)
    band = np.concatenate([
        img[y0:y1, max(0, xi - R):max(0, xi - 2 * s)].ravel(),
        img[y0:y1, min(img.shape[1], xi + 2 * s):min(img.shape[1], xi + R)].ravel(),
    ])
    rms = _clipped_rms(band) if band.size else np.nan
    pk = float(img[cy, cx])
    snr = pk / rms if rms and np.isfinite(rms) and rms > 0 else np.inf

    dra = (coord.ra - expected.ra).to_value("arcsec") * np.cos(
        expected.dec.rad)
    ddec = (coord.dec - expected.dec).to_value("arcsec")
    return PsfFitResult(
        coord=coord, dra_arcsec=float(dra), ddec_arcsec=float(ddec),
        snr=snr, peak=pk, curvature_ok=curv_ok,
        sigma_ra_est=(float(np.hypot(20.0 / snr, 1.8))
                      if np.isfinite(snr) else np.nan),
        sigma_dec_est=(float(np.hypot(12.0 / snr, 0.25))
                       if np.isfinite(snr) else np.nan),
    )
