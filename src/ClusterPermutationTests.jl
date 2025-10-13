module ClusterPermutationTests

using Reexport: @reexport
using StatsAPI: StatsAPI, fit, params, pvalue
using StatsModels: RegressionModel, @formula, FormulaTerm, Term, coefnames,
            InteractionTerm, FunctionTerm, ConstantTerm
import GLM: lm, LinearModel
using HypothesisTests: HypothesisTests, EqualVarianceTTest, OneSampleTTest,
                        UnequalVarianceTTest
import MixedModels
using PrettyTables: ft_printf, pretty_table, tf_unicode_rounded
using ProgressMeter: Progress, next!
using UnPack: @unpack
using Random
using TypedTables: Table, columns

include("StudyDesigns/StudyDesigns.jl")
@reexport using .StudyDesigns

export ClusterPermutationTest,
    CPTTest, CPPairedSampleTTest, CPEqualVarianceTTest, CPUnequalVarianceTTest,
    CPRegressionModel, CPLinearModel,
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
    npermutations,
    reset,
    cluster_table,
    params,
    fits,
    fit,
    resample!,
    pvalues,
    summary,
    @formula # reexport

include("StudyDesigns/utilities.jl")
include("utilities.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("lm.jl")

end;
