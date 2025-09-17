module ClusterPermutationTests

import DataAPI
using StatsAPI: StatsAPI, fit, params, pvalue, nobs
using CategoricalArrays: CategoricalValue, CategoricalArray, categorical
using StatsModels: RegressionModel, @formula, FormulaTerm, Term,
            InteractionTerm, FunctionTerm, ConstantTerm
import GLM: lm
using HypothesisTests: HypothesisTests, EqualVarianceTTest, OneSampleTTest,
                        UnequalVarianceTTest
import MixedModels
using PrettyTables: ft_printf, pretty_table, tf_unicode_rounded
using ProgressMeter: Progress, next!
using Random: Random, AbstractRNG, shuffle
import Tables
using TypedTables: TypedTables, Table, columnnames, columns
using UnPack: @unpack

export ClusterPermutationTest,
    CPTTest, CPPairedSampleTTest, CPEqualVarianceTTest, CPUnequalVarianceTTest,
    CPRegressionModel, CPLinearModel,
    # PermuteDesign
    PermutationDesign, BetweenDesign, WithinDesign, MixedDesign,
    UnitObs, NoUnitObs,
    variables_within,
    variables_between,
    has_variable,
    get_variable,
    is_within,
    unit_obs,
    design_table,
    nrow,
    cell_indices,
    get_variable,
    shuffle_variable!,
    shuffle_variable,
    # data,
    CPData,
    select_rows,
    epoch_length,
    nepochs,
    nobs,
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


include("utilities.jl")
include("perm_design/perm_design.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("lm.jl")


DataFrames(::PermutationDesign; kwargs...) = throw(ArgumentError("Have you loaded DataFrames?"))
if !isdefined(Base, :get_extension)
    include("../ext/DataFramesExt.jl")
end

end;
