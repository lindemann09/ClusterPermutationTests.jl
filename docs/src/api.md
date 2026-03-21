# ClusterPermutationTests.jl

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

### Methods

```@docs
cluster
cluster_mass_stats
cluster_nhd
cluster_pvalues
cluster_table
npermutations
resample!
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
ClusterDefinition
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

TODO


