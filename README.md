# AI Use for Schoolwork, Learning Dispositions and Response Style

Replication code for the article "AI Use for Schoolwork, Learning Dispositions and Response Style: A Person-Centred Study of 15-Year-Olds in Beijing, Shanghai, Jiangsu and Zhejiang (China)".

The code reproduces every estimate reported in the article and its supporting information from the PISA 2025 public-use files. The results are written as CSV files; the tables and figures of the article were formatted from these files. No data are included in this repository.

## Data

The analyses use the PISA 2025 public-use files distributed by the OECD (PISA data and methodology pages, https://www.oecd.org/pisa/data/). Download the following SPSS files:

| File | Content |
|---|---|
| `CY09_MS_STU_PUF.sav` | Student questionnaire, plausible values and weights |
| `CY09_MS_SCH_PUF.sav` | School questionnaire |
| `CY09_MS_STU_TT_PUF.sav` | Student questionnaire timing |
| `CY09_MS_LDW_20260806.sav` | Learning in the digital world: cognitive item responses |

Place them in `data/raw/`, or edit `config.txt` to point to their location. File names that differ in a later release can also be changed in `config.txt`. All analyses are restricted to Beijing, Shanghai, Jiangsu and Zhejiang (China) (`CNT == "QCI"`).

## Software

The analyses were run with R 4.6.1 and Python 3.14.

- R packages: nanoparquet (0.5.2), poLCA (1.6.0.2), mice (3.19.0), multilevLCA (2.1.7). Install them with `Rscript install_r_packages.R`.
- Python packages: listed in `requirements.txt` (`pip install -r requirements.txt`).

The weighted latent class models, the three-step and BCH estimators and the models with direct effects are implemented in the repository (`code/lca_functions.R` and the scripts that use it); poLCA is used only to check the estimator on unweighted data.

## Running the analyses

From the repository root:

```
python run_all.py
```

The script runs all steps in order and writes a log for each to `output/logs/`. Rscript is found on the PATH or through the `RSCRIPT` environment variable. Single steps can be run with `python run_all.py --only 07` or resumed with `python run_all.py --from 07`; every step after 03 needs the class solution written by step 03.

All generated files are written to `output/`:

- `output/derived/`: analysis files built from the PISA data (not to be redistributed)
- `output/results/`: estimates as CSV files

The complete run takes several hours on a desktop computer, most of it for the class enumeration (step 03) and the replicate-weight estimates.

## Scripts and outputs

| Step | Script | Content | Used in |
|---|---|---|---|
| 01 | `01_prepare_data.py` | Analysis file: AI-use items, dispositions, correlates, achievement, weights, timing, LDW takers | all |
| 02 | `02_response_style.py` | Extreme-response and acquiescence indices for the item sets; reliability | Table 3; Section 3.2 |
| 03 | `03_enumerate_classes.R` | Weighted latent class models with 1-8 classes, bivariate residuals, poLCA check, split-half replication | Table 2, Table S2, Figure 2, Figure S1 |
| 04 | `04_careless_responding.R` | Flags for possibly careless respondents; refits without them | Table S7; Section 4.6 |
| 05 | `05_correlates_three_step.R` | Three-step multinomial models (main, careless respondents excluded, response style and age added) | Figure 3, Tables S3, S9 |
| 06 | `06_correlates_imputation.R` | Three-step model with multiply imputed correlates | Section 4.6, Figure S3 |
| 07 | `07_outcomes_bch.R` | BCH class means of dispositions and achievement; Cohen's d | Figure 4, Table S4 |
| 08 | `08_outcomes_composite.R` | Disposition composite and adjusted versions, science, computational problem solving, exclusions by criterion | Figure 5, Table S4, Figure S3 |
| 09 | `09_response_style_bch.R` | Response style by class; share of the disposition difference left after each adjustment | Table 3, Figure 4 |
| 10 | `10_unweighted_multilevel.R` | Unweighted and multilevel latent class models | Section 4.6, Figure S3 |
| 11 | `11_within_school.R` | Class differences within schools (fixed effects) | Section 4.5, Figure S3 |
| 12 | `12_school_composition.R` | Profile composition across schools | Figure S2 |
| 13 | `13_adjacent_solutions.R` | Five-, six- and seven-class solutions compared | Table S8 |
| 14 | `14_continuous_and_ers_correlates.R` | Continuous intensity and breadth models; correlates of the ERS indices | Table S8; Section 4.4 |
| 15 | `15_local_dependence.R` | Six-class models with direct effects for locally dependent item pairs | Table S10; Section 4.1 |
| 16 | `16_descriptives.py` | Weighted descriptive statistics, correlations and the list of response-style items | Tables 1, S1, S6 |

Numbers reported only in the text (for example classification accuracy or class counts) are printed in the corresponding logs.

## Methods in brief

- Latent class models are estimated by weighted pseudo-maximum likelihood with the final student weight rescaled to the sample size; students with at least one valid answer to the four AI-use items are included (missing answers are handled by full-information estimation).
- Standard errors use the 80 BRR-Fay replicate weights (Fay factor 0.5). For each replicate the latent class model, the classification-error matrix and the second-step estimates are recomputed.
- Achievement outcomes are analysed separately for each of the 10 plausible values and combined with Rubin's rules.
- Random starts are seeded in each script, so a complete run reproduces the reported estimates.

## License

The code is released under the MIT License (see `LICENSE`). The PISA data are subject to the OECD's terms of use and are not covered by this license.
