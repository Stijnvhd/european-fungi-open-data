# Open-data analyses of European fungal fruiting, nitrogen deposition and species naming

Code and result tables supporting three manuscripts:

1. Later every autumn: four decades of fungal fruiting phenology in north-western Europe, corrected for observer effort and compared with satellite canopy senescence.
2. Cleaner air, returning fungi: the fruiting of nitrogen-sensitive ectomycorrhizal genera tracks declining nitrogen deposition across Europe, 1990 to 2023.
3. Two thirds unnamed: the pace of fungal species description and the named share of the ITS reference library.

Author: Stijn Van Hijfte, Galactic Automation BV, Zottegem, Belgium (stijn@galacticautomation.com).

## Layout

- `code/` Python scripts. Run from the repository root; they read `data/` and write `out/`.
- `results/` the output tables and JSON summaries that the manuscripts report (copies of `out/`).

## Data sources (all open)

- GBIF occurrence search API (https://api.gbif.org/v1/occurrence/search), faceted counts only, queried 24 October 2026. Dataset keys for the national recording schemes are in `results/dataset_names.json`.
- EMEP MSC-W rv5.6 2025 trend runs (type 2), total nitrogen deposition 1990 to 2023, via the MET Norway THREDDS server.
- GADM 4.1 administrative boundaries.
- NASA POWER daily point data (T2M, PRECTOTCORR), 1981 to 2025, at the 108 sites in `results/sites.csv`.
- MODIS MCD12Q2 v6 land surface phenology via the ORNL DAAC subset API, 2001 to 2023.
- Species Fungorum Plus, April 2024, ChecklistBank dataset 2073 (Catalogue of Life data package).
- UNITE general FASTA release for Fungi, version 10 (4 April 2024), classification table via the Swedish Biodiversity Data Infrastructure.

## Pipeline

| Step | Script | Writes |
| --- | --- | --- |
| Download fungal and plant counts by ten-day period | `pull_dekads.py`, `pull_effort.py`, `pull_records.py` | `data/dekads_*.parquet` |
| Season metrics, effort correction, trends | `phenology.py` | `out/phenology_results.json`, `out/phenology_species_slopes.csv`, `out/phenology_metrics.parquet` |
| Weather drivers (NASA POWER) | `pull_sites.py`, `climate_drivers.py` | `out/climate_results.json`, `out/country_climate.csv` |
| Satellite phenology and coupling | `pull_modis.py`, `decoupling.py` | `out/decoupling_results.json`, `out/host_coupling.json`, `out/country_leaf_phenology.csv` |
| Indicator genera by district and year | `pull_genus.py`, `pull_ngroups*.py` | `data/ngroups_gadm2*.parquet` |
| EMEP deposition | `pull_emep.py` | `data/emep_ndep_halfdeg.parquet` |
| Dose curve, trajectories, exceedance, projections | `nitrogen.py`, `nitrogen2.py`, `nitrogen_post.py` | `out/core_nitrogen_results.json`, `out/*ndep*`, `out/exceedance_fixed.json` |
| Legacy test | `nitrogen_legacy.py` | `out/legacy_results.json`, `out/legacy_table.csv` |
| Species description rates and UNITE named share | `naming.py`, `naming2.py` | `out/naming_summary.json`, `out/naming_scenarios.json`, `out/unite_summary.json`, `out/descriptions_per_year.csv` |
| Figures | `fig_*.py`, `figstyle.py` | `figures/*.png` |

## Software

Python 3.13; pandas 3.0.5, numpy 2.5.3, pyarrow 26.0.0, statsmodels 0.15.0, pyfixest 0.60.0, scipy 1.18.1, matplotlib 3.11.2, netCDF4 1.7.5, shapely 2.2.0, requests 2.34.2.

## Licence

Code: MIT. Result tables: CC BY 4.0. Source data remain under the licences of their providers.
