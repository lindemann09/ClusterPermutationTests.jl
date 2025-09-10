module ClusterPermutationTests

using DataAPI
using DataFrames
using StatsModels
using Random
using UnPack
using HypothesisTests
using ProgressMeter
using PrettyTables
using StatsAPI

import HypothesisTests: nobs
import StatsAPI: fit, summary, params

export ClusterPermutationTest,
    ClusterPermutationTTest,
    # PermuteDesign
    PermutationDesign,
    design_table,
    nrow,
    cell_indices,
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

include("perm_design/perm_design.jl")
include("cpdata.jl")
include("cluster.jl")
include("cptest/cptest.jl")
include("sampling.jl")

include("ttest.jl")

end;
