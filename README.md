# dsa110-voltage-imaging

Interferometric localization of DSA-110 FRBs from M8-era raw voltage dumps.

Correlates the candidate archive's voltage fragments with the
[dsa110-bbproc](https://github.com/dsa110/dsa110-bbproc) toolkit
(pre-correlation dedispersion), builds calibrated measurement sets, ties the
astrometry to RFC VLBI calibrators, and produces a corrected FRB position.

The pre-M8 (T3-era) version of this repo — `toolkit_dev` correlator,
`T2_<cand>.json` metadata, `Level2/calibration` layout, `batch.bash` — is
preserved on the `legacy` branch.

## Layout

- `dsavim/` — python package. `metadata.py` reads the M8 candidate archive
  (`/dataz/dsa110/candidates/<name>/`: `Level3/<name>.json` C2 record,
  `Level2/voltages/` fragments + manifests, `filterbank/` provenance,
  `calibration/` ±2 h correlator HDF5s hard-linked by C3) into a `Burst`
  object. `voltages_to_ms.py`, `generate_uvh5.py`, `uvh5_to_ms.py` do the
  correlation and MS conversion.
- `scripts/run_localization.py` — **the driver** (pythonized batch.bash).
- `scripts/voltages_to_ms.py`, `create_ms.py`, `bp_gain_cal.py`,
  `create_all_icrs.py`, `calibrate_image_icrs.py` — pipeline stages invoked
  by the driver (each runnable standalone).
- `notebooks/localize.ipynb` — the localization analysis, run after the
  driver.

## Usage (h23, `casa38` environment)

```bash
conda activate casa38
cd scripts
python run_localization.py <candname>
```

The driver measures the burst time bin itself (a coherent filterbank formed
with pre-correlation dedispersion), derives the science-MS offsets, and runs
every stage through the wsclean images; see `python run_localization.py -h`
for stage selection and overrides. Key options:

- `--dm` — refine beyond the search's `dm_median`.
- `--ntint` — science-MS integration (native 32.768 µs samples; default 512).
- `--bpcal-source` — bandpass calibrator, by VLA-calibrator-list name
  (default `2253+161` = 3C454.3). The latest per-transit MS from the
  realtime calibration system (`/dataz/dsa110/operations/calibration`) is
  copied in; if none exists it is built from the
  `/dataz/dsa110/operations/correlator` HDF5s. Or pass an existing MS with
  `--bpcal-ms`.
- `--bpcal-flux-jy`, `--refant` — calibration knobs (defaults 15 Jy, ant 104).

Outputs land in `/media/ubuntu/ssd/localization_processing/<candname>/`,
including `localization_run.json` (all derived parameters). RFC calibrator
counts are limited to 10 everywhere.

Then open `notebooks/localize.ipynb`, set `trigname`, and run all cells for
the bandpass QC, NVSS astrometry check, VLBI offset fit, and the corrected
FRB position.
