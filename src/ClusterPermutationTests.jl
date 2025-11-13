module ClusterPermutationTests

using Reexport: @reexport
using StatsAPI: StatsAPI, fit, coef
using StatsBase: stderror, CoefTable, quantilerank
using StatsModels
import GLM: lm, LinearModel
using HypothesisTests: HypothesisTests, EqualVarianceTTest, OneSampleTTest,
                        UnequalVarianceTTest, HypothesisTest
import MixedModels: LinearMixedModel, is_randomeffectsterm
using PrettyTables: ft_printf, pretty_table, tf_unicode_rounded
using ProgressMeter: Progress, next!
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
    select_epochs,
    epoch_length,
    nepochs,
    design_table,
    # Cluster
    ClusterCriterium,
    ClusterDefinition,
    cluster_mass_stats,
    cluster_ranges,
    cluster_pvalues,
    cluster_table,
    cluster_nhd,
    time_series_stats,
    time_series_fits,
    npermutations,
    reset,
    fit,
    coefnames,
    resample!,
    summary,
    @formula, # reexport
    # plotting
    plot_time_series_stats!,
    plot_cluster_mass_stats_distribution!

include("StudyDesigns/utilities.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("regressionmodels.jl")
include("mixedmodels.jl")



## Makie extensions
_makie_error() = throw(ArgumentError("Have you loaded an appropriate Makie backend?"))
plot_time_series_stats!(::Any, ::Any; kwargs...) = _makie_error()
plot_cluster_mass_stats_distribution!(::Any, ::Any; kwargs...) = _makie_error()

end;
