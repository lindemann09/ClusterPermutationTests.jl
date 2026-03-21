# ClusterPermutationTests.jl

[![Build Status](https://github.com/lindemann09/ClusterPermutationTests.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lindemann09/ClusterPermutationTests.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-blue)](https://julialang.org)

A Julia package for **cluster-based permutation testing** on time-series data.
Cluster permutation tests are a nonparametric approach widely used in neuroscience
(EEG, ERP, fMRI) and psycholinguistics to identify significant temporal effects while
naturally controlling the family-wise error rate across many time points.

**Documentation**: [https://lindemann09.github.io/ClusterPermutationTests.jl](https://lindemann09.github.io/ClusterPermutationTests.jl)

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

## Related Software

- [HypothesisTests.jl](https://github.com/JuliaStats/HypothesisTests.jl) — t-test implementations used internally
- [MixedModels.jl](https://github.com/JuliaStats/MixedModels.jl) — mixed-effects models used internally for `CPMixedModel`

## Citation

If you use ClusterPermutationTests.jl in published work, please cite the package:

> Lindemann, O. (2024). *ClusterPermutationTests.jl* [Computer software].
> GitHub. https://github.com/lindemann09/ClusterPermutationTests.jl

## License

MIT License. See [LICENSE](LICENSE) for details.
