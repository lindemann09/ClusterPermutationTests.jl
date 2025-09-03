module ClusterPermutationTests

using DataAPI
using Tables
using TypedTables
using StatsModels
using Random
using UnPack
using HypothesisTests
using ProgressMeter
using PrettyTables
using StatsAPI

import StatsAPI: fit, summary, params
import Random: randperm!

export ClusterPermutationTestDefinition,
    ClusterPermutationTest,
    ClusterPermutationTTest,
    PermutationDesign,
#    CPData,
    ClusterCriteria,
    ClusterDefinition,
    cluster_statistics,
    cluster_ranges,
    unit_obs,
    unit_obs_name,
    design_table,
    data_matrix,
    epoch_length,
    nepochs,
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
