# PREPARE intervention study

Preregistration and power analysis for a randomised crossover study of whether access to AI-generated practice questions improves student performance. The study is embedded in the course Statistik I (BSc Psychology, University of Bern) in the fall semester of 2026. Approximately 400 students are enrolled.

Jamie Cummins, Boris Mayer, Sandra Grinschgl, Alexander Mainetti, Malte Elson, and Michael Schulte-Mecklenbeck.

## Study overview

PREPARE is a pipeline which uses a large language model to generate multiple-choice practice exercises (with feedback and explanations) from lecture transcripts. Students work through the exercises in a web interface and never interact with the model directly. All students are individually randomised once into two groups; in alternating topic-weeks, one group receives the week's exercises before the weekly knowledge assessment and the other group receives them afterwards. Every student therefore receives the complete set of materials every week, and each week has a concurrent randomised control group on the same topic. The primary question is whether early access improves performance on the weekly assessment (intention-to-treat). Secondary questions concern the association between PREPARE usage and performance in the weekly assessments and the final exam.

The canonical preregistration document is `preregistration/preregistration_PREPARE.docx`, and the registration itself is available at (REF).

## Repository structure

```
preregistration/
  preregistration_PREPARE.docx   canonical preregistration document
  power_fit.R                    preregistered RQ1 model and fallback ladder (sourced by the two scripts below)
  power_simulation.R             simulation-based power analysis over the design grid
  power_sensitivity.R            sensitivity checks for assumptions the historical data cannot verify
  power_analysis/                outputs of the two scripts (summaries, raw per-trial results, curves, logs, sessionInfo)
analysis/
  historical_performance.R       per-lesson performance in the HS2025 cohort (ceiling check and calibration)
  historical_performance/        aggregate outputs of that script
```

## Reproducing the power analysis

The scripts require R (4.6 was used) with `lme4`, `lmerTest`, `dplyr`, `tidyr`, and `ggplot2`; the exact package versions are recorded in `preregistration/power_analysis/sessionInfo.txt`. Run the scripts from the repository root:

```bash
Rscript preregistration/power_simulation.R
```

```bash
Rscript preregistration/power_sensitivity.R
```

The full simulation (108 design cells, 500 simulated trials per cell) took about 15 minutes on 10 cores, and the sensitivity checks (39 cells, 400 trials each) about 6 minutes. Both scripts are seeded (`master_seed = 20260905`) and write their outputs to `preregistration/power_analysis/`. Use `--quick` for a one-minute smoke test of the simulation, and `--sims N` or `--cores N` to change the number of trials or cores.

## Data which are not in this repository

The variance components of the simulation were calibrated against the HS2025 weekly assessment exports from Qualtrics (`data/old-reference-data/`, ignored by git). These exports contain self-chosen student pseudonyms and response timestamps and are therefore not shared; the same applies to the per-attempt file `analysis/historical_performance/first_attempts_long.csv`. The aggregate outputs of `analysis/historical_performance.R` (per-lesson summaries, maximum-score diagnostics, and the calibration notes) are included. The ethics application and participant documents are likewise kept out of the repository.
