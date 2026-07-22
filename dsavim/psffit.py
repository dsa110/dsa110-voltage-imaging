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
  worse than the parabola on any event. That is this module.

* Errors are FORMAL, computed per image as Cov = H^-1 C_g H^-1 (see
  ``fit_position``), with the signal Hessian measured from the actual
  PSF core and the noise-gradient covariance from off-source pixels.
  Validated 2026-07-22 on 8 fixed-chain B1933 pulses (S/N 24-45):
  measured scatter RA 0.68" / Dec 0.30" vs mean formal 0.93" / 0.38"
  — tracks the thermal scatter, mildly conservative. At S/N <~ 15 the
  plateau threshold effect adds error the formal figure cannot see.

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
    #: FORMAL 1-sigma centroid errors (arcsec), from
    #: Cov = H^-1 C_g H^-1 with H = peak x (PSF Hessian at its core,
    #: measured from the psf image when supplied, else from the image
    #: quadratic) and C_g = the measured off-source noise-GRADIENT
    #: covariance. For a Gaussian beam in beam-correlated noise this
    #: reduces to the classic FWHM/(2.355 SNR) ~ 0.42 FWHM/SNR. Purely
    #: statistical: excludes the low-S/N threshold effect (sidelobe /
    #: plateau ambiguity) and all calibration/tie systematics.
    sigma_ra_est: float
    sigma_dec_est: float
    #: correlation coefficient of the (RA, Dec) errors
    sigma_corr: float = float("nan")
    #: True when the PSF image supplied the signal Hessian
    formal_from_psf: bool = False


def _quad_fit(img: np.ndarray, cy: int, cx: int, q: int):
    """LS 2-D quadratic about (cy, cx); returns (coef[6], H[2,2])."""
    yy, xx = np.mgrid[-q:q + 1, -q:q + 1]
    patch = img[cy - q:cy + q + 1, cx - q:cx + q + 1]
    A = np.column_stack([np.ones(xx.size), xx.ravel(), yy.ravel(),
                         xx.ravel() ** 2, yy.ravel() ** 2,
                         (xx * yy).ravel()])
    coef, *_ = np.linalg.lstsq(A, patch.ravel(), rcond=None)
    H = np.array([[2 * coef[3], coef[5]], [coef[5], 2 * coef[4]]])
    return coef, H


def _grad_cov(strips, rms):
    """Covariance of the noise gradient from off-source 2-D strips.

    Pixels brighter than 5 x rms (field sources) are masked before the
    sample covariance so a bright confusing source cannot inflate it.
    """
    gxs, gys = [], []
    for s in strips:
        if s.shape[0] < 8 or s.shape[1] < 8:
            continue
        gy, gx = np.gradient(s)
        ok = np.abs(s) < 5 * rms
        ok &= np.isfinite(gx) & np.isfinite(gy)
        gxs.append(gx[ok])
        gys.append(gy[ok])
    if not gxs:
        return None
    gx = np.concatenate(gxs)
    gy = np.concatenate(gys)
    if gx.size < 500:
        return None
    return np.cov(np.vstack([gx, gy]))


def _clipped_rms(a: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    cut = np.percentile(np.abs(a), 95)
    return float(np.std(a[np.abs(a) < cut]))


def fit_position(
    image: Union[str, "np.ndarray"],
    expected: SkyCoord,
    *,
    wcs: Optional[WCS] = None,
    psf: Optional[Union[str, "np.ndarray"]] = None,
    search_asec: float = 25.0,
    quad_halfwidth: int = 3,
) -> PsfFitResult:
    """Constrained-argmax + 2-D quadratic position fit with FORMAL errors.

    Args:
        image: wsclean image FITS path (restored or dirty), or a 2-D
            array with ``wcs`` supplied.
        expected: prior position (the C2/C1 detection position for
            FRBs). The argmax is constrained to ``search_asec`` of it —
            essential on the plateau-dominated snapshot beam, where an
            unconstrained peak can lock onto a distant noise bump.
        psf: matching ``-psf.fits`` path or array. When supplied, the
            signal Hessian in the error covariance comes from the PSF
            core (exact, noise-free); otherwise from the image's own
            quadratic fit. ALWAYS pass it for cleaned/restored images:
            there the source has the restoring-Gaussian shape (much
            steeper than the snapshot dirty beam whose matched-filter
            peak actually set the position), so the image-quadratic
            fallback underestimates sigma several-fold (measured 5x in
            RA on the B1933 pulse set). The fallback is only trustworthy
            on niter=0 / dirty images.
        search_asec: half-width of the prior box. Match it to the
            upstream position uncertainty (~half a detection-grid cell
            for C2 positions).
        quad_halfwidth: half-width (pixels) of the 2-D quadratic
            refinement neighbourhood. 3 validated best at the 2"/px
            production scale (use ~10 at 0.3"/px).

    Formal error model: the fitted position solves grad(I) = 0 near the
    peak, with I = A B(x - x0) + n. Linearising, the error covariance is

        Cov(x0_hat) = H^-1 C_g H^-1,
        H   = A grad^2 B(0)       (signal-peak Hessian),
        C_g = Cov[grad n]         (noise-gradient covariance),

    with both measured from the data: H from the PSF core (times the
    fitted amplitude) and C_g from off-source pixels of THIS image, so
    the actual noise correlation (weighting, uv coverage, cleaning
    state) is accounted for without assumptions. Gaussian-beam +
    beam-correlated-noise limit: sigma = FWHM / (2.355 SNR).
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

    coef, H_img = _quad_fit(img, cy, cx, q)
    curv_ok = bool(np.all(np.linalg.eigvalsh(H_img) < 0))
    if curv_ok:
        dx, dy = np.linalg.solve(H_img, [-coef[1], -coef[2]])
        dx = float(np.clip(dx, -q, q))
        dy = float(np.clip(dy, -q, q))
    else:
        dx = dy = 0.0

    coord = wcs.pixel_to_world(cx + dx, cy + dy)

    # off-source strips well away from the search box: rms + noise-
    # gradient covariance for the formal errors
    R = max(4 * s, 50)
    y0, y1 = max(0, yi - R), min(img.shape[0], yi + R)
    strips = [
        img[y0:y1, max(0, xi - R):max(0, xi - 2 * s)],
        img[y0:y1, min(img.shape[1], xi + 2 * s):min(img.shape[1], xi + R)],
    ]
    band = np.concatenate([s_.ravel() for s_ in strips])
    rms = _clipped_rms(band) if band.size else np.nan
    pk = float(img[cy, cx])
    snr = pk / rms if rms and np.isfinite(rms) and rms > 0 else np.inf

    # ---- formal centroid covariance:  Cov = H^-1 C_g H^-1 ----
    sig_x = sig_y = corr = float("nan")
    from_psf = False
    if np.isfinite(rms):
        H_sig = None
        if psf is not None:
            p = (fits.open(psf)[0].data.squeeze().astype(np.float64)
                 if isinstance(psf, str) else
                 np.asarray(psf, dtype=np.float64))
            pcy, pcx = np.unravel_index(np.nanargmax(p), p.shape)
            if (q < pcx < p.shape[1] - q and q < pcy < p.shape[0] - q):
                _, H_psf = _quad_fit(p, pcy, pcx, q)
                if np.all(np.linalg.eigvalsh(H_psf) < 0):
                    H_sig = pk * H_psf / p[pcy, pcx]
                    from_psf = True
        if H_sig is None and curv_ok:
            H_sig = H_img
        C_g = _grad_cov(strips, rms)
        if H_sig is not None and C_g is not None:
            try:
                Hinv = np.linalg.inv(H_sig)
                cov = Hinv @ C_g @ Hinv          # pixels^2
                if cov[0, 0] > 0 and cov[1, 1] > 0:
                    sig_x = float(np.sqrt(cov[0, 0]) * scale)
                    sig_y = float(np.sqrt(cov[1, 1]) * scale)
                    corr = float(cov[0, 1] /
                                 np.sqrt(cov[0, 0] * cov[1, 1]))
            except np.linalg.LinAlgError:
                pass

    dra = (coord.ra - expected.ra).to_value("arcsec") * np.cos(
        expected.dec.rad)
    ddec = (coord.dec - expected.dec).to_value("arcsec")
    return PsfFitResult(
        coord=coord, dra_arcsec=float(dra), ddec_arcsec=float(ddec),
        snr=snr, peak=pk, curvature_ok=curv_ok,
        sigma_ra_est=sig_x, sigma_dec_est=sig_y,
        sigma_corr=corr, formal_from_psf=from_psf,
    )
