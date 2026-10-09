# SNPMD replication package

R replication code for the paper "Nonparametric Minimum-Distance Estimation of Simulation Models".

## Contents

| Location | Description |
| --- | --- |
| `MA_MC_Exp/` | MA(1) experiments varying polynomial degree and grid size; `_logP` scripts use 50 logarithmically spaced grid sizes from 21 to 2,001. |
| `Ass_MC_Exp/` | MA(1) diagnostics for objective-function errors, Assumption 2, and Remark 4. |
| `BH_MC_Exp/` | Four-parameter Brock–Hommes experiments: 16 configurations and `N = M = 250, 500, 1000`. |
| `BH_Emp_App/` | S&P 500 and Euro Stoxx 50 estimation, saved results/data, and objective-function plots. |
| `BH_functions.R`, `Mom_functions.R`, `data_helpers.R` | Model, likelihood, moment, and data functions. |

## Requirements and setup

Install R and the required CRAN packages:

```r
install.packages(c("randtoolbox", "tseries", "xtable", "quantmod", "zoo"))
```

The `parallel` package ships with R. After downloading the repository, replace author-specific paths and `"~"`/`"..."` placeholders, check `source()` locations, and create required output directories.

## Replication workflow

1. **MA(1):** run the selected `Est_MA_NPMD_P_K*.R` and `Est_MA_MLE_II_SMM.R`, then the matching `Read_Exp2_MA_SNPMD*.R`. Repeat for `(N,M) = (100,100), (100,1000), (1000,1000)`, updating settings and filenames.
2. **Assumption checks:** run `Est_MA_SNPMD_hetero.R` for those three combinations, then `MC_MA1_Assumption2_Remark4.R`. The default holds the observed series fixed; commented code enables the unconditional experiment.
3. **Brock–Hommes:** run the estimation scripts in `BH_MC_Exp/`, followed by `Read_BH_SNPMD_MC_exp.R`. NPSMLE writes its own summaries.
4. **Empirical application:** run both `Est_BH_*_emp.R` scripts, then `Plot_Est_BH_emp.R`. The supplied `.rds` files and analysis CSVs also support plotting without re-estimation.

The MA experiments default to 10,000 replications; the Brock–Hommes experiments default to 1,000. For a small SNPMD trial, run from the repository root (macOS/Linux):

```bash
BH_QUICK=1 BH_COMBINATIONS=1 BH_CORES=1 \
  Rscript BH_MC_Exp/Est_BH_SNPMD_MC_exp.R
```

Configure Brock–Hommes runs with `BH_REPLICATIONS`, `BH_CORES`, and `BH_OUTPUT_DIR`; use `BH_SAMPLE_SIZES`/`BH_COMBINATIONS` for Monte Carlo experiments and `BH_YEARS` for empirical windows. Keep trial outputs separate.

## Empirical data and outputs

Empirical defaults use nested 250-, 500-, and 1,000-observation windows ending on 31 December 2012, with `M=N`. Prices are detrended using a trailing 61-observation moving average. Supply dated price CSVs through `BH_SP500_FILE` and `BH_EURO_STOXX50_FILE`, or enable Yahoo Finance downloads. Included analysis CSVs contain transformed data.

Outputs include estimates, summaries, R objects, and figures/tables. Empirical `mcse` measures simulation uncertainty in the replication mean.

## Setup

- Set `BH_ROOT_DIR` explicitly to the repository root for empirical NPSMLE. Configure both sections of `Plot_Est_BH_emp.R` to locate the saved results and data.
- MA post-processing also needs definitions for `iR` in the benchmark summary and `vTicks` in the logarithmic-grid reader.

`BH_functions.R` adapts functions from [Sylvain Barde's workshop replication package](https://github.com/Sylvain-Barde/abm-validation-cef-2023); attribution is retained in the source.
