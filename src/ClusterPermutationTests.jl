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
using Tables: getcolumn, Tables
using TypedTables
using UnPack: @unpack

export ClusterPermutationTest,
    CPTTest, CPPairedSampleTTest, CPEqualVarianceTTest, CPUnequalVarianceTTest,
    CPRegressionModel, CPLinearModel,
    # PermuteDesign
    PermutationDesign, BetweenDesign, WithinDesign, MixedDesign,
    UnitObs, NoUnitObs,
    names,
    names_within,
    names_between,
    names_covariates,
    has_variable,
    is_covariate,
    is_within,
    is_between,
    unit_observation,
    cell_indices,
    shuffle_variable!,
    shuffle_variable,
    getcolumn,
    # data,
    CPData,
    select_rows,
    epoch_length,
    nepochs,
    nobs,
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


include("utilities.jl")
include("perm_design/perm_design.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("lm.jl")

end;
