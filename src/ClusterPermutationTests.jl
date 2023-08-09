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
import Random: randperm!

export ClusterPermutationTestDefinition,
    ClusterPermutationTest,
    ClusterPermutationTTest,
#    CPData, CPDesign,
    ClusterCriteria,
    ClusterDefinition,
    cluster_statistics,
    cluster_ranges,
    unit_obs,
    unit_obs_name,
    design_table,
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

include("cpdesign.jl")
include("cpdata.jl")
include("cpdefinition.jl")
include("cluster.jl")
include("cpcollection.jl")
include("clusterpermutationtest.jl")

include("general_test.jl")
include("ttest.jl")

end;
