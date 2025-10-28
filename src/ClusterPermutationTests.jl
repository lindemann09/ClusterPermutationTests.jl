module ClusterPermutationTests

using Reexport: @reexport
using StatsAPI: StatsAPI, fit, pvalue, coef
using StatsBase: coeftable, stderror
using StatsModels
import GLM: lm, LinearModel
using HypothesisTests: HypothesisTests, EqualVarianceTTest, OneSampleTTest,
                        UnequalVarianceTTest, HypothesisTest
import MixedModels: LinearMixedModel, is_randomeffectsterm
using PrettyTables: ft_printf, pretty_table, tf_unicode_rounded
using ProgressMeter: Progress, next!
using UnPack: @unpack
using Logging: with_logger, NullLogger
using Random
using Tables: columntable, getcolumn
using TypedTables: Table, columnnames

include("StudyDesigns/StudyDesigns.jl")
@reexport using .StudyDesigns

export ClusterPermutationTest,
    CPTTest, CPPairedSampleTTest, CPEqualVarianceTTest, CPUnequalVarianceTTest,
    CPRegressionModel, CPLinearModel, CPMixedModel,
    # data,
    CPData,
    select_rows,
    epoch_length,
    nepochs,
    design_table,
    # Cluster
    ClusterCriterium,
    ClusterDefinition,
    cluster_stats,
    cluster_ranges,
    cluster_pvalues,
    cluster_table,
    sample_stats,
    npermutations,
    permutation_stats,
    reset,
    fit,
    initial_fits,
    resample!,
    summary,
    @formula, # reexport
    # plotting
    plot_sample_stats!,
    plot_sample_distribution!

include("StudyDesigns/utilities.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("regressionmodels.jl")



## Makie extensions
_makie_error() = throw(ArgumentError("Have you loaded an appropriate Makie backend?"))
plot_sample_stats!(::Any, ::Any; kwargs...) = _makie_error()
plot_sample_distribution!(::Any, ::Any; kwargs...) = _makie_error()

end;
