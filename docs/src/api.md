# API

## Tests
```@docs
CPPairedSampleTTest
CPEqualVarianceTTest
CPUnequalVarianceTTest
CPRegressionModel
CPMixedModel
CPAnovaMixedModel

CPTTest
CPLinearModel
```

### Configuration
```@docs
CPConfig
```

### Methods

```@docs
cluster
cluster_mass_stats
cluster_pvalues
cluster_table
resample!
permutation_distribution

npermutations
ncoefs

time_series_fits
time_series_stats
```

## Data
```@docs
CPData
```

### Methods

```@docs
convert_to_cpdata
design_table
epoch_length
nepochs
select_epochs
```

## Cluster

```@docs
ClusterCriterium
ClusterPredefined
```

## Study Design

```@docs
BetweenDesign
WithinDesign
MixedDesign
UnitObs
NoUnitObs
```

### Methods

```@docs
study_design
unit_observation
names_between
names_within
names_covariates
is_between
is_within
is_covariate
has_variable
shuffle_variable!
shuffle_variable
```


## Methods from `StatsAPI.jl`, `DataAPI.jl`
```julia
coefnames
fit
ncol
nobs
nrow
@formula
```

## Plotting

An optional plotting extension for `Makie.jl` are available, if the package is loaded.

TODO


