module ClusterPermutationTests

using DataAPI
using TypedTables
using StatsModels
using Random
using UnPack
using HypothesisTests
using ProgressMeter
using PrettyTables

import StatsAPI: fit, summary

export ClusterPermutationTestDefinition,
    ClusterPermutationTest,
    ClusterPermutationTTest,
#    CPData,
    ClusterCriteria,
    ClusterDefinition,
 #   PermuteDesign,
    cluster_statistics,
    cluster_ranges,
    unit_obs,
    unit_obs_name,
    iv_names,
    data_matrix,
    nepoch_samples,
    nepochs,
    npermutations,
    reset,
    cluster_table,
    sample_statistics,
    fits,
    fit,
    resample!,
    pvalues,
    summary,
    @formula # reexport

include("permutation.jl")
include("cpdata.jl")
include("cpdefinition.jl")
include("cluster.jl")
include("cpcollection.jl")
include("clusterpermutationtest.jl")

include("general_test.jl")
include("ttest.jl")

end;
