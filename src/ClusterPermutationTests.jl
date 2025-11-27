module ClusterPermutationTests

using Reexport: @reexport
using SplitApplyCombine: groupfind
using Random
using Tables: columntable, getcolumn, istable, matrix
using TypedTables: Table, columnnames
using ProgressMeter: Progress, next!, finish!
using Logging: with_logger, NullLogger, SimpleLogger, AbstractLogger, global_logger
using CSV

import DataAPI: nrow, ncol
import StatsAPI: StatsAPI, fit, nobs
using StatsBase: coef, stderror, CoefTable, quantilerank
using StatsModels
import GLM: lm, LinearModel
using HypothesisTests: HypothesisTests, EqualVarianceTTest, OneSampleTTest,
	UnequalVarianceTTest, HypothesisTest
import MixedModels: LinearMixedModel, is_randomeffectsterm, refit!
using AnovaMixedModels: anova, teststat


include("StudyDesigns/StudyDesigns.jl")
@reexport using .StudyDesigns

export ClusterPermutationTest,
	CPTTest, CPPairedSampleTTest, CPEqualVarianceTTest, CPUnequalVarianceTTest,
	CPRegressionModel, CPLinearModel, CPMixedModel, CPAnovaMixedModel,
	# DataAPI, StatsAPI
	nrow, ncol,
	nobs, fit, coefnames, @formula,
	# data,
	CPData,
	select_epochs,
	epoch_length,
	nepochs,
	design_table,
	convert_to_cpdata,
	# Cluster
	ClusterCriterium,
	ClusterDefinition,
	cluster_mass_stats,
	cluster,
	cluster_pvalues,
	cluster_table,
	cluster_nhd,
	time_series_stats,
	time_series_fits,
	npermutations,
	resample!,
	# plotting
	plot_time_series_stats!,
	plot_cluster_nhd!


include("StudyDesigns/utilities.jl")
include("cluster.jl")
include("cpdata.jl")
include("cptype.jl")
include("sampling.jl")

include("ttest.jl")
include("regressionmodels.jl")
include("mixedmodels.jl")
include("anovamixedmodels.jl")


## Makie extensions
_makie_error() = throw(ArgumentError("Have you loaded an appropriate Makie backend?"))
plot_time_series_stats!(::Any, ::Any; kwargs...) = _makie_error()
plot_cluster_nhd!(::Any, ::Any; kwargs...) = _makie_error()

end;
