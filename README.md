# Intracellular INT-formazan is a dye-dependent optical modulator of nucleic-acid fluorescence in single marine prokaryoplankton cells

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22757989-blue.svg)](https://doi.org/10.5281/zenodo.22757989)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.txt)

Analysis pipeline and per-cell datasets for Villegas-Mendoza *et al.*

**DOI (this release, v1.0.0):** [10.5281/zenodo.22757989](https://doi.org/10.5281/zenodo.22757989) — the version cited in the manuscript.
**DOI (all versions):** [10.5281/zenodo.22757988](https://doi.org/10.5281/zenodo.22757988) — always resolves to the latest release.
**Raw image data:** BioImage Archive, accession **S-BIAD4069**

---

## What this repository is — and what it is not

This repository contains the image-analysis pipeline and the per-cell tables from which every figure and statistical table of the associated manuscript can be reproduced.

**It does not measure the respiration rate of individual cells.** The experimental design deliberately separates two scales: community respiration is measured independently for each biochemical oxygen demand (BOD) bottle with planar oxygen optodes, whereas fluorescence descriptors are extracted from individual segmented cells. In the manuscript's own words, the analysis asks

> whether a bottle-level physiological state is reflected in single-cell optical phenotypes, without assuming that the respiration rate of any individual cell has been measured.

What the pipeline quantifies is a set of **per-cell optical descriptors** — intensity, heterogeneity, spatial-pattern and shape — on RAW confocal images and on their MSSR-transformed counterparts, and fits mixed-effects models with bottle as a random effect. The central result is that intracellular INT-formazan shifts nucleic-acid fluorescence in **opposite directions** depending on the dye (down with DAPI, up with SYBR Green I), and that the apparent SYBR Green I response **reverses again** between confocal and STED. Any predictive use of these descriptors as a respiration proxy remains an open question, not a result of this work.

---

## Dataset at a glance

| | |
|---|---|
| BOD bottles in the experiment | 9 (0.39–2.69 µM O₂ h⁻¹, two thermal clusters) |
| Bottles passing image-level quality control | 7 (0.39–2.41 µM O₂ h⁻¹) |
| Confocal acquisitions (`.oib`) | 557 |
| Source images entering the per-cell table | 468 |
| Segmented cells | 1,507 |
| Descriptors per cell | 37 |
| Mixed-effects models (37 descriptors × 2 fluorochromes) | 74, of which **73 converge** |

Confocal imaging: inverted **Olympus FluoView FV1000** on an IX81 body, Laboratorio Nacional de Microscopía Avanzada (LNMA, IBT-UNAM), UPLSAPO 60× (NA = 1.30), 800 × 800 px at 12-bit, native `.oib`, lateral pixel size 52.9 nm.
STED imaging: inverted **Leica TCS SP8 STED 3X-FALCON**, SLN lab, ICFO (Castelldefels, Barcelona), 456 × 456 px at 8-bit, native `.lif`, lateral pixel size 31.9 nm.

Full acquisition parameters are given in the manuscript (§1.2, §1.3) and in Supplementary Methods.

---

## Repository layout

```
MSSR_PIPELINE/          analysis package — one module per stage
  config.py               paths, calibration, parameters, experimental design
                          ← the only file that normally needs editing
  mssr_core.py            MSSR algorithm (spatial and temporal)
  run_pipeline.py         orchestrator
  stage0_convert.py       batch .oib → .tif conversion, archived by condition
  stage1_segmentation.py  Otsu segmentation, 96 × 96 px crop extraction
  stage2_mssr.py          MSSR applied to the crops
  stage3_parameters.py    FWHM (Gaussian fit), peak detection, nearest neighbour
  stage4_statistics.py    descriptive statistics and models
  stage5_heatmaps.py      intensity heatmaps on a common global scale
  README.md               detailed guide to the package (in Spanish)

notebooks/
  MSSR_Pipeline_260327.ipynb   Colab runner for the 260327 session; all logic
                               lives in MSSR_PIPELINE/, the notebook only calls it

data/                   per-cell tables (see the data dictionary below)
  features.csv            1,507 cells × 49 columns
  bottles.csv             9 bottles: temperature and measured respiration rate
  int_samples.csv         mapping of each INT sample to its BOD bottle of origin
  feature_results.csv     74 fitted models (coefficients, p-values, FDR q-values)
```

Stages run in sequence **0 → 1 → 2 → 3 → {4, 5}**. Stage 0 is only needed when starting from raw `.oib` files; if TIFFs are already organised in the expected hierarchy, start at stage 1.

---

## Data dictionary

### `data/features.csv` — one row per segmented cell

| Group | Columns |
|---|---|
| Identifiers | `session`, `image_uid`, `image_id`, `cell_id` |
| Shape | `area_px`, `area_um2`, `perimeter_um`, `eq_diameter_nm`, `major_axis_nm`, `minor_axis_nm`, `aspect_ratio`, `solidity`, `eccentricity` |
| RAW intensity | `raw_mean`, `raw_median`, `raw_sum`, `raw_max`, `raw_min`, `raw_std` |
| RAW heterogeneity | `raw_cv`, `raw_entropy`, `raw_skew`, `raw_kurtosis` |
| RAW spatial pattern | `raw_n_local_maxima`, `raw_n_local_minima`, `raw_frac_dark`, `raw_frac_bright` |
| MSSR counterparts | the same families, prefixed `mssr_` |

Intensities are in raw 12-bit grey levels, with **no normalisation or rescaling** at any point in the pipeline.

### `data/bottles.csv`

`bottle`, `temperature_C`, `respiration_uM_O2_h` (optode-measured community respiration), `is_control`, `notes`.

### `data/int_samples.csv`

`int_sample`, `bottle_bod` (the BOD bottle each INT aliquot came from), `notes`.

### `data/feature_results.csv` — one row per fitted model

`feature`, `fluorophore`, `n_cells`, `n_bottles`, `converged`, `beta_treatment`, `p_treatment`, `q_treatment`, `beta_respiration`, `p_respiration`, `beta_interaction`, `p_interaction`, `q_interaction`.

Models have the form `feature ~ C(treatment) * respiration + (1 | bottle)`. The `q_` columns are Benjamini–Hochberg FDR-adjusted p-values. One of the 74 models does not converge (`converged = False`) and is excluded from the reported analysis, leaving 73.

---

## Raw images

The raw confocal and STED image data are **not** in this repository. They are deposited in the **BioImage Archive** under accession **S-BIAD4069**: 557 native Olympus confocal images, the two Leica STED containers, and the derived analysis tables.

The deposit is private until its release date and can be shared with reviewers on request.

---

## Reproducing the analysis

```bash
git clone https://github.com/jvillegas18/mssr-confocal-respiration
cd mssr-confocal-respiration
pip install -r requirements.txt
```

Requires **Python 3.10+**.

**From the per-cell tables** (no images needed) — this reproduces every figure and statistical table of the manuscript:

```bash
cd MSSR_PIPELINE
python run_pipeline.py 4 5     # statistics + heatmaps
```

**From raw images**, after downloading S-BIAD4069: edit `BASE_DIR` and `RAW_OIB_DIR` in `config.py`, then run a dry pass of stage 0 before committing to the conversion:

```python
import stage0_convert as s0
s0.convert_all(dry_run=True)    # classifies and prints the shape of each .oib, writes nothing
s0.convert_all(dry_run=False)
```

```bash
python run_pipeline.py 1 2 3 4 5
```

`MSSR_PIPELINE/README.md` documents the expected folder hierarchy and every parameter in `config.py`.

---

## Citation

If you use this software or these data, please cite the manuscript and the archived release. Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).

> Villegas-Mendoza J, Tiznado-Ramos F, Maske-Rubach H, Cajal-Medrano R, Marsal-Terés M, Loza-Álvarez P, Guerrero A. *Intracellular INT-formazan is a dye-dependent optical modulator of nucleic-acid fluorescence in single marine prokaryoplankton cells.* `[PENDING: journal reference once accepted]`

> Villegas-Mendoza J, Tiznado-Ramos F, Maske-Rubach H, Cajal-Medrano R, Marsal-Terés M, Loza-Álvarez P, Guerrero A. (2026). *Analysis pipeline for: Intracellular INT-formazan is a dye-dependent optical modulator of nucleic-acid fluorescence in single marine prokaryoplankton cells* (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.22757989

## Authors

| | Affiliation | ORCID |
|---|---|---|
| **Josué Villegas-Mendoza**\* | Facultad de Ciencias Marinas, UABC, Ensenada, Mexico | [0000-0001-5614-3375](https://orcid.org/0000-0001-5614-3375) |
| Fernando Tiznado-Ramos | Facultad de Ciencias Marinas, UABC, Ensenada, Mexico | [0009-0008-0103-2959](https://orcid.org/0009-0008-0103-2959) |
| Helmut Maske-Rubach | CICESE, Ensenada, Mexico | [0000-0002-8047-6484](https://orcid.org/0000-0002-8047-6484) |
| Ramón Cajal-Medrano | CICESE, Ensenada, Mexico | [0009-0001-7172-3148](https://orcid.org/0009-0001-7172-3148) |
| Maria Marsal-Terés | ICFO – Institut de Ciències Fotòniques, Castelldefels, Spain | [0000-0001-9678-8717](https://orcid.org/0000-0001-9678-8717) |
| Pablo Loza-Álvarez | ICFO – Institut de Ciències Fotòniques, Castelldefels, Spain | [0000-0002-3129-1213](https://orcid.org/0000-0002-3129-1213) |
| **Adán Guerrero**\* | LNMA, Instituto de Biotecnología, UNAM, Cuernavaca, Mexico | [0000-0002-4389-5516](https://orcid.org/0000-0002-4389-5516) |

\* Corresponding authors.

## License

MIT — see [`LICENSE.txt`](LICENSE.txt).
