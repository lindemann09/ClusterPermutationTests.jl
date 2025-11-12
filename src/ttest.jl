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
#### coefnames
####
StatsAPI.coefnames(cpt::CPTTest) = [string(cpt.cpc.iv)]


####
#### Parameter_estimates
####
"""returns vector (time) of vector (parameters)"""
@inline function parameter_estimates(cpt::CPTTest,
	design::StudyDesign;
	fit_cluster_only::Bool = true,
	store_fits::Bool = false)::TVecTimeXParameter

	# Estimate parameters for a specific cluster (range)
	epochs, design_tbl = _prepare_data(cpt, cpt.dat.epochs, design) # TODO view?
	param = TVecTimeXParameter()
	if fit_cluster_only
		time_points = cpt.cpc.tp
	else
		time_points = Int32(1):Int32(epoch_length(cpt.dat))
	end
	for t in time_points
		tt = _estimate(cpt, view(epochs, :, t), design_tbl)
		push!(param, [tt.t])
		if store_fits
			push!(cpt.cpc.m, tt)
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


####
#### T-test have just a single coefficient
####
time_series_coefs(cpt::CPTTest) = time_series_coefs(cpt, 1)
cluster_mass_permutations(cpt::CPTTest) = cluster_mass_permutations(cpt, 1)
cluster_ranges(cpt::CPTTest) = cluster_ranges(cpt, 1)
cluster_mass(cpt::CPTTest) = cluster_mass(cpt, 1)
cluster_pvalues(cpt::CPTTest; inhibit_warning::Bool = false) = cluster_pvalues(cpt, 1; inhibit_warning)
cluster_table(cpt::CPTTest) = cluster_table(cpt, 1)
