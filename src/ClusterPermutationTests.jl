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

export ClusterPermutationTestDefinition,
    ClusterPermutationTest,
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
    ClusterCriteria,
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

include("cpdesign.jl")
include("cpdata.jl")
include("cpdefinition.jl")
include("cluster.jl")
include("cpcollection.jl")
include("clusterpermutationtest.jl")

include("general_test.jl")
include("ttest.jl")

end;
