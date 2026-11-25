# oneshot-shelflife

Code for the article

**Shelf-Life Estimation From Staggered Sensory Tests Through One-Shot Device Methods: Reanalysis of Pezik Pickle and Green Tea Data**
Coşkun Kuş, Department of Statistics, Faculty of Science, Selçuk University, Konya, Türkiye.

Interactive application (runs in the browser, no installation): https://coskunkus.github.io/oneshot-shelflife/

## Method

In a staggered sensory shelf-life study, each evaluation shows only whether the product had already become unacceptable on the day of the test. Such responses are one-shot device test data: an unacceptable response is left-censored and an acceptable one is right-censored. The code implements the non-parametric maximum likelihood estimator of the probability of acceptance (EM algorithm), the semi-parametric proportional hazards model with free hazard ratios or a log-linear link in a covariate such as storage temperature, the expected information matrix, logit, log and log(-log) confidence intervals, the Pearson goodness-of-fit statistic, and prediction at an untested covariate value (Zhu and Balakrishnan, 2022, Reliability Engineering and System Safety 221:108319). The Weibull hazard method is implemented in its original version (Gacula, 1975), in which all responses are ranked, and in the version of Labuza and Schmidl (1988) and Fu and Labuza (1993), followed by Cardelli and Labuza (2001), in which only the unacceptable responses are ranked for comparison.

## Files

| File | Content |
|---|---|
| `code/core_oneshot.R` | Shared core: every other script sources it |
| `code/run_pezik.R` | Pezik pickle data: Tables 1, 3, 4 and Figures 1 to 3 of the article |
| `code/run_greentea.R` | Green tea data: Tables 2, 5, 6 and Figure 4 of the article |
| `code/run_all.R` | Reproduces the examples of Zhu and Balakrishnan (2022) |
| `code/check_formulas.R` | Independent numerical checks of the core (exact and Monte Carlo) and of the two food analyses |
| `code/export_shinylive.R` | Converts the Shiny application into the static site in `docs/` |
| `shiny/app.R` | Shiny application (uses `code/core_oneshot.R`) |
| `shiny/example_greentea.csv` | Example input file for the application |
| `data/pezik.csv`, `data/greentea.csv` | The two data sets of the article (Keklik et al., 2017; Xiao et al., 2022), in the input format of the application |

## Requirements

R (version 4.3 or later) is the only software needed. The analysis scripts use base R only, need no internet connection and no other package, and run on Windows, macOS and Linux. The Shiny application needs the package `shiny`; the export of the browser version needs `shinylive` and `httpuv`.

## Running

All analysis scripts use base R only. Scripts locate the project root themselves, so they can be run from anywhere:

    source("code/check_formulas.R")   # verification
    source("code/run_all.R")          # examples of Zhu and Balakrishnan (2022)
    source("code/run_pezik.R")        # pezik pickle analysis
    source("code/run_greentea.R")     # green tea analysis

Reports are written to `results/` and figures (eps and pdf) to `figures/`. The Monte Carlo layer of `check_formulas.R` uses the seed 20261122.

The Shiny application needs the package `shiny`: open `shiny/app.R` in RStudio and press Run App. The browser version is produced by `code/export_shinylive.R` (packages `shinylive` and `httpuv`).

## Input format for the application

A CSV file with one row per group and test day: `group`, `x` (numeric covariate such as storage temperature; optional), `day`, `N` (number of evaluations) and `r` (number of unacceptable responses).

## Notation shared with the article

`IT` test days, `N` number of evaluations, `r` number of unacceptable responses, `S1` baseline probability of acceptance, `lambda` hazard ratio, `beta` log-linear link coefficient, `m` latent numbers of rejections in the EM algorithm.

## Data

Keklik NM, Işıklı ND, Sur EB (2017). Food Science and Technology 37:125-130. doi:10.1590/1678-457X.33216 (CC BY).
Xiao M, Liu S, Jin H, Xiao M, Wang H, Zhang H, Dai Q (2022). Journal of Chemistry 2022:2955839. doi:10.1155/2022/2955839 (CC BY).
