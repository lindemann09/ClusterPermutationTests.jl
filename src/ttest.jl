abstract type CPTTest <: ClusterPermutationTest end
abstract type CPTwoSampleTTest <: CPTTest end

struct CPPairedSampleTTest <: CPTTest
	cpc::CPCollection{OneSampleTTest}
	dat::CPData
	compare::Tuple
end;

struct CPEqualVarianceTTest <: CPTwoSampleTTest
	cpc::CPCollection{EqualVarianceTTest}
	dat::CPData
	compare::Tuple
end;

struct CPUnequalVarianceTTest <: CPTwoSampleTTest
	cpc::CPCollection{UnequalVarianceTTest}
	dat::CPData
	compare::Tuple
end;

n_threads_default(::CPTTest) = Threads.nthreads()

####
#### Test info
####
function test_info(x::CPTTest)
	return "$(typeof(x)), (compare=$(string.(x.compare)), $(x.cpc.mass_fnc))"
end

####
#### Fit TTests
####
function StatsAPI.fit(T::Type{<:CPTTest}, # TODO: two value comparison only, needs to be more general
	iv::SymbolOString,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum)

	iv = Symbol(iv)
	paired = is_within(dat.design, iv)
	if T == CPTTest
		# choose test based on design
		T = paired ? CPPairedSampleTTest : CPEqualVarianceTTest
	end

	if paired
		T == CPPairedSampleTTest || throw(ArgumentError(
			"$(iv) is a within-subject variable. A t-test for independent samples ($(T)) is not possible."))
		compare = unique(getproperty(dat.design.within, iv))
	else
		T == CPPairedSampleTTest && throw(ArgumentError(
			"$(iv) is a between-subject variable. $(T) is not possible."))
		compare = unique(getproperty(dat.design.between, iv))
	end

	length(compare) == 2 || throw(ArgumentError(
		"'$iv' comprises $(length(compare)) categories; two required."))

	if T == CPPairedSampleTTest
		M = OneSampleTTest
	elseif T != CPEqualVarianceTTest
		M = EqualVarianceTTest
	elseif T == CPUnequalVarianceTTest
		M = UnequalVarianceTTest
	else
		throw(ArgumentError("Test $(T) not supported for t.test."))
	end
	cpc = CPCollection{M}(iv, mass_fnc, cluster_criterium)

	rtn = T(cpc, dat, (compare[1], compare[2]))

	fit_initial_time_series!(rtn)
	return rtn
end;

# formula interface
function StatsAPI.fit(T::Type{<:CPTTest},
	f::FormulaTerm,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	kwargs...)

	(f.lhs isa Term && f.rhs isa Term) || throw(
		ArgumentError("Incorrect t.test formula: '$f'"))
	return fit(T, Symbol(f.rhs), dat, cluster_criterium; kwargs...)
end

####
#### Time series statistics
####
time_series_stats(cpt::CPTTest) = get_parameter(cpt.cpc.ts, 1)

####
#### Parameter_estimates
####
@inline function parameter_estimates(cpt::CPTTest,
	design::StudyDesign,
	cl_ranges::Vector{TClusterRange};
	store_fits::Bool = false)::TParameterVector
	# Estimate parameters for a specific cluster (range)
	epochs, design_tbl = _prepare_data(cpt, cpt.dat.epochs, design) # TODO view?
	param = TParameterVector()
	for cr in cl_ranges
		for t in cr
			tt = _estimate(cpt, view(epochs, :, t), design_tbl)
			if store_fits
				push!(cpt.cpc.m, tt)
				push!(cpt.cpc.ts, TPStats(t, tt.t))
			else
				push!(param, TPStats(t, tt.t))
			end
		end
	end
	return param
end

####
#### CPPairedSampleTTest
####
@inline function _prepare_data(cpt::CPPairedSampleTTest,
	epochs::Matrix{<:Real},
	permutation::StudyDesign)::Tuple{Matrix{eltype(epochs)}, Table}

	iv = getcolumn(permutation, cpt.cpc.iv)
	tbl = Table((; cpt.cpc.iv => iv))
	a = @view epochs[iv .== cpt.compare[1], :]
	b = @view epochs[iv .== cpt.compare[2], :]
	return b - a, tbl # equal size required
end

@inline function _estimate(::CPPairedSampleTTest, values::SubArray{<:Real}, ::Table)
	return OneSampleTTest(values)
end

####
#### CPTwoSampleTTest
####

@inline function _prepare_data(cpt::CPTwoSampleTTest,
	epochs::Matrix{<:Real},
	permutation::StudyDesign)::Tuple{Matrix{eltype(epochs)}, Table}

	return epochs, Table((; cpt.cpc.iv => getproperty(permutation, cpt.cpc.iv)))
end

@inline function _estimate(cpt::CPEqualVarianceTTest,
	values::SubArray{<:Real},
	permutation::Table)

	(dat_a, dat_b) = _ttest_get_data(cpt, values, permutation)
	return EqualVarianceTTest(dat_a, dat_b)
end

@inline function _estimate(cpt::CPUnequalVarianceTTest,
	values::SubArray{<:Real},
	permutation::Table)

	(dat_a, dat_b) = _ttest_get_data(cpt, values, permutation)
	return UnequalVarianceTTest(dat_a, dat_b)
end

@inline function _ttest_get_data(cpt::CPTTest, values::SubArray{<:Real}, permutation::Table)
	# perform sequential ttests -> parameter
	iv = getproperty(permutation, cpt.cpc.iv)
	dat_a = @view values[iv .== cpt.compare[1]]
	dat_b = @view values[iv .== cpt.compare[2]]
	return dat_a, dat_b
end

##
## permutation stats
##
permutation_stats(cpt::CPTTest) = _permutation_stats(cpt.cpc, cluster_ranges(cpt), 1)

##
## Cluster Functions
##
cluster_ranges(cpt::CPTTest) = _cluster_ranges(time_series_stats(cpt), cpt.cpc.cc)

function cluster_mass(cpt::CPTTest)
	ts = time_series_stats(cpt)
	cl_ranges = _cluster_ranges(ts, cpt.cpc.cc)
	return _cluster_mass(cpt.cpc.mass_fnc, ts, cl_ranges)
end

function cluster_pvalues(cpt::CPTTest; inhibit_warning::Bool = false)
	return _cluster_pvalues(permutation_stats(cpt), cluster_mass(cpt), inhibit_warning)
end

cluster_table(cpt::CPTTest) = _cluster_table(time_series_stats(cpt),cluster_ranges(cpt), cluster_pvalues(cpt))
