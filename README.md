# ClusterPermutationTests.jl

[![Build Status](https://github.com/lindemann09/ClusterPermutationTests.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lindemann09/ClusterPermutationTests.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Julia](https://img.shields.io/badge/Julia-1.9%2B-blue)](https://julialang.org)

A Julia package for **cluster-based permutation testing** on time-series data.
Cluster permutation tests are a nonparametric approach widely used in neuroscience
(EEG, ERP, fMRI) and psycholinguistics to identify significant temporal effects while
naturally controlling the family-wise error rate across many time points.

## Overview

The test proceeds in three steps:

1. **Fit a statistical model at every time point** — t-test, linear regression, or
   linear mixed model.
2. **Detect clusters** — contiguous windows where the test statistic exceeds a
   threshold.
3. **Build a null-hypothesis distribution via permutation** — shuffle the
   condition labels, refit, and record the maximum cluster mass. Repeat thousands
   of times. The p-value of each observed cluster is its quantile rank against
   this distribution.

Compared to mass-univariate correction methods (Bonferroni, FDR), the cluster
approach exploits the temporal autocorrelation that is intrinsic to physiological
time series and typically has greater power when effects span multiple consecutive
time points.

## Installation

```julia
using Pkg
Pkg.add("ClusterPermutationTests")
```

## Quick Example

```julia
using ClusterPermutationTests, CSV, TypedTables, Downloads

# --- Load data ---
# `epochs` is an (n_epochs × n_timepoints) matrix of observed values.
# `design` is a table with one row per epoch describing the experimental conditions.
epochs = CSV.File(download("…/epochs.dat"), header=false)
design = CSV.read(download("…/design.csv"), Table)

dat = CPData(epochs, design; unit_obs=:subject_id)

# --- Define cluster criterion ---
# Clusters = consecutive time points where |t| ≥ 1.69, spanning at least 50 samples.
cc = ClusterCriterium(threshold=1.69, min_size=50)

# --- Fit the initial time series and detect clusters ---
cpt = fit(CPPairedSampleTTest, @formula(y ~ condition), dat, cc)

# --- Build the null-hypothesis distribution ---
resample!(cpt, 5000; use_threads=true)

# --- Inspect results ---
cluster_table(cpt)
```

```
 cluster │ from  to   size  mass stats  Pr(>|z|)  sign
─────────┼──────────────────────────────────────────────
 1.1     │  210  274    65    -749.6      0.041    *
 1.2     │  350  869   520  -13669.8     <0.001    *
```

## Supported Statistical Models

| Type | Constructor |
|---|---|
| Paired-samples t-test | `CPPairedSampleTTest` |
| Independent-samples t-test (equal variances) | `CPEqualVarianceTTest` |
| Welch's t-test (unequal variances) | `CPUnequalVarianceTTest` |
| OLS linear regression | `CPLinearModel` |
| Linear mixed-effects model | `CPMixedModel` |
| ANOVA F-test from mixed model | `CPAnovaMixedModel` |

All model types share the same interface: `fit` → `resample!` → `cluster_table`.

## Data Input

### `CPData` — the data container

`CPData` holds an `(n_epochs × n_timepoints)` epoch matrix together with a design
table describing the experimental conditions for each epoch.

```julia
# From two files / CSV sources
dat = CPData(epochs_matrix, design_table; unit_obs=:subject_id)

# Restrict to a subset of conditions at construction time
dat = CPData(epochs, design; unit_obs=:subject_id, condition=["A", "B"])

# Subset an existing CPData object
dat2 = select_epochs(dat; condition="A", subject=:all)

# Inspect
epoch_length(dat)   # number of time points
nepochs(dat)        # number of observations (rows)
design_table(dat)   # retrieve the design as a TypedTable
```

### `convert_to_cpdata` — from long-format data

If your data are in long format (one row per observation × time bin), use
`convert_to_cpdata`:

```julia
# data has columns: subject_id, condition, bin (time index), response (DV)
dat = convert_to_cpdata(data; unit_obs=:subject_id, bin=:bin, response=:response)
```

## Study Designs

ClusterPermutationTests automatically classifies your design variables as
**between-subjects**, **within-subjects**, or **covariates** by inspecting which
variables vary across the unit of observation.

```julia
design = study_design(table; unit_obs=:subject_id)

names_between(design)      # between-subject variables
names_within(design)       # within-subject variables
names_covariates(design)   # numeric covariates
unit_observation(design)   # the unit-of-observation column
```

During permutation, within-subject variables are shuffled independently per
participant; between-subject variables are shuffled globally. This preserves the
within-subject correlation structure.

## Fitting Models

All model types are constructed via `fit`:

```julia
# T-tests  (iv must name a binary condition variable)
cpt = fit(CPPairedSampleTTest, @formula(y ~ condition), dat, cc)
cpt = fit(CPEqualVarianceTTest, @formula(y ~ condition), dat, cc)
cpt = fit(CPUnequalVarianceTTest, @formula(y ~ condition), dat, cc)

# OLS regression
cpt = fit(CPLinearModel, @formula(y ~ condition + covariate), dat, cc)

# Linear mixed model
cpt = fit(CPMixedModel,
          @formula(y ~ condition + (1 | subject_id)),
          dat, cc;
          reml=true)

# ANOVA F-statistics from a mixed model
cpt = fit(CPAnovaMixedModel,
          @formula(y ~ condition + (1 | subject_id)),
          dat, cc)
```

The `shuffle_ivs` keyword (regression models) lets you explicitly name the
predictor(s) to permute if the automatic detection is not what you want:

```julia
cpt = fit(CPMixedModel, formula, dat, cc; shuffle_ivs=[:condition])
```

### Cluster mass function

By default the cluster mass is the **sum** of test statistics within a cluster.
Pass `mass_fnc` to use a different aggregation:

```julia
cpt = fit(CPLinearModel, formula, dat, cc; mass_fnc=maximum)
```

## Permutation Resampling

```julia
# Run 5 000 permutations using all available threads
resample!(cpt, 5_000; use_threads=true)

# Limit to a specific number of threads
resample!(cpt, 5_000; use_threads=4)

# Append more permutations to an already-resampled object
resample!(cpt, 5_000)   # total becomes 10 000

npermutations(cpt)      # how many permutations have been accumulated
```

A minimum of **1 000 permutations** is required to compute p-values; at least
**5 000** is recommended for reliable estimates.

## Querying Results

All functions below accept an `effect` argument that can be an integer index, a
`Symbol`, or a `String` matching a coefficient name (see `coefnames`). When a
model has only one effect, `effect` can be omitted.

| Function | Returns |
|---|---|
| `time_series_stats(cpt, effect)` | Per-time-point test statistics (length = epoch_length) |
| `cluster(cpt, effect)` | Detected cluster ranges (`Vector{UnitRange}`) |
| `cluster_mass_stats(cpt, effect)` | Aggregate mass statistic per cluster |
| `cluster_nhd(cpt, effect)` | Null-hypothesis distribution matrix (n_permutations × n_clusters) |
| `cluster_pvalues(cpt, effect)` | Monte Carlo permutation p-values |
| `cluster_table(cpt)` | Summary `CoefTable` for all effects |
| `cluster_table(cpt, effect)` | Summary table for a single effect |
| `time_series_fits(cpt)` | Vector of fitted model objects (one per time point) |

```julia
# Example: retrieve results for a named effect
cluster_pvalues(cpt, :condition)
cluster_mass_stats(cpt, :condition)
cluster(cpt, :condition)

# Inspect the full summary
summary(cpt)
```

## Cluster Specification

### Automatic detection — `ClusterCriterium`

```julia
ClusterCriterium(threshold=2.0, min_size=10, use_absolute=true)
```

| Argument | Description | Default |
|---|---|---|
| `threshold` | Minimum statistic value to enter a cluster | (required) |
| `min_size` | Minimum number of consecutive time points | `10` |
| `use_absolute` | Compare `|statistic|` to threshold | `true` |

A common choice for `threshold` is the critical value of the *t*-distribution at
the desired α-level (e.g., `t_{0.05, df} ≈ 1.69` for large *df*).

### Predefined windows — `ClusterDefinition`

If the time windows of interest are known a priori:

```julia
cc = ClusterDefinition([200:400, 600:900])   # two fixed windows
cc = ClusterDefinition(200:800)              # single window
```

## Visualization

Optional plotting support is available via a Makie extension. Load any Makie
backend (e.g., `GLMakie`, `CairoMakie`) to activate it.

```julia
using CairoMakie, ClusterPermutationTests

fig = Figure(size=(800, 400))
plot_time_series_stats!(Axis(fig[1, 1]), cpt, :condition)   # time series of t-values

fig2 = Figure()
plot_cluster_nhd!(fig2, cpt; effect=:condition, bins=100)   # null distribution histograms
```

`plot_cluster_nhd!` draws one histogram per detected cluster, with the observed
cluster mass overlaid as a vertical line.

## Full Workflow Example

```julia
using ClusterPermutationTests, CSV, TypedTables, Downloads

# 1. Load data
fl_design = "https://raw.githubusercontent.com/lindemann09/JuliaDataSets/refs/heads/main/data/cpt1_design.csv"
fl_epochs = "https://raw.githubusercontent.com/lindemann09/JuliaDataSets/refs/heads/main/data/cpt1_epochs.dat"

epochs = CSV.File(download(fl_epochs), header=false)
design = CSV.read(download(fl_design), Table)
dat    = CPData(epochs, design; unit_obs=:subject_id)

# 2. Define cluster criterion
cc = ClusterCriterium(threshold=1.69, min_size=50)

# 3a. Paired t-test
cpt_t = fit(CPPairedSampleTTest, @formula(y ~ operator_str), dat, cc)
resample!(cpt_t, 5_000; use_threads=true)
cluster_table(cpt_t)

# 3b. Linear mixed model
cpt_mm = fit(CPMixedModel,
             @formula(y ~ operator_str + (1 | subject_id)),
             dat, cc; reml=true)
resample!(cpt_mm, 5_000; use_threads=true)
cluster_table(cpt_mm)

# 3c. ANOVA F-test
cpt_amm = fit(CPAnovaMixedModel,
              @formula(y ~ operator_str + (1 | subject_id)),
              dat, cc)
resample!(cpt_amm, 5_000; use_threads=true)
cluster_table(cpt_amm)
```

## Related Software

- [MixedModels.jl](https://github.com/JuliaStats/MixedModels.jl) — mixed-effects models used internally for `CPMixedModel`
- [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) — t-test implementations used internally
- [clusterperm (R)](https://github.com/dalejbarr/clusterperm) — R reference implementation

## Citation

If you use ClusterPermutationTests.jl in published work, please cite the package
and the methodological reference:

> Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-
> and MEG-data. *Journal of Neuroscience Methods, 164*(1), 177–190.
> https://doi.org/10.1016/j.jneumeth.2007.03.024

## License

MIT License. See [LICENSE](LICENSE) for details.
