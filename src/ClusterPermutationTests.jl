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
    cluster_statistics,
    cluster_ranges,
    cluster_pvalues,
    cluster_table,
    cluster_parameter,
    npermutations,
    reset,
    sample_statistics,
    fit,
    initial_fits,
    resample!,
    summary,
    @formula # reexport


include("StudyDesigns/utilities.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("regressionmodels.jl")

end;
