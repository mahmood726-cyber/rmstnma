Data Ingestion and Standard Format

Goal: ensure all NMA datasets follow a common schema so they can be consumed consistently by rmstnma and other tooling.

Standard columns (required):
- study_id: character, study identifier
- treatment: character, treatment/arm label
- time: numeric, time in months (>= 0)
- survival: numeric, Kaplan–Meier survival in [0, 1]

Optional columns (recommended):
- source: character (e.g., km, ipdfromkm)
- arm_id: character/integer per-study arm identifier
- n_risk: numeric number at risk

Where to put raw files:
- Place CSV/TSV files in data-raw/raw/. One file per dataset.
- Column names will be auto-standardized if you use common alternatives like study, trt, arm, t, S, surv.

Build pipeline:
1) Add raw files to data-raw/raw/
2) Run: source("data-raw/build_datasets.R")
3) On success, standardized .rda files are written to data/

Validation:
- Validation checks presence of required columns, types, value ranges, duplicate keys, and basic connectivity.
- Strict mode also checks monotone non-increasing survival per arm and survival≈1 near time≈0.

Tips for harmonization:
- Ensure time is in months; convert beforehand if needed.
- Ensure no duplicate (study_id, treatment, time) rows.
- Use consistent treatment labels across studies.

