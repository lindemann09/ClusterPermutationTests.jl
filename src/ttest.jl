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
function StatsAPI.fit(T::Type{<:CPTTest},
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
	cpc = CPCollection{M}([iv], mass_fnc, cluster_criterium)

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
StatsAPI.coefnames(::CPTTest) = ["contrast"]


####
#### Parameter_estimates
####
"""returns vector (time) of vector (parameters)"""
@inline function parameter_estimates(cpt::CPTTest,
	design::AbstractStudyDesign,
	time_points::Vector{Int32};
	store_fits::Bool = false)::TVecTimeXParameter

	# Estimate parameters for a specific cluster (range)
	T = typeof(cpt)
	epochs, design_tbl = _prepare_data(T, cpt.dat.epochs, design,
			first(cpt.cpc.shuffle_ivs), cpt.compare)
	iv = first(cpt.cpc.shuffle_ivs)
	param = TVecTimeXParameter()
	for t in time_points
		tt = _estimate(T, view(epochs, :, t), design_tbl, iv, cpt.compare)
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
@inline function _prepare_data(::Type{CPPairedSampleTTest},
	epochs::Matrix{<:Real},
	permutation::AbstractStudyDesign,
	iv::Symbol,
	compare::Tuple)::Tuple{Matrix{eltype(epochs)}, Table}

	iv_dat = getcolumn(permutation, iv)
	tbl = Table((; iv => iv_dat))
	a = @view epochs[iv_dat .== compare[1], :]
	b = @view epochs[iv_dat .== compare[2], :]
	return b - a, tbl # equal size required
end

@inline function _estimate(::Type{CPPairedSampleTTest}, values::SubArray{<:Real}, ::Table, ::Symbol, ::Tuple)
	return OneSampleTTest(values)
end

####
#### CPTwoSampleTTest
####

@inline function _prepare_data(::Type{<:CPTwoSampleTTest},
	epochs::Matrix{<:Real},
	permutation::AbstractStudyDesign,
	iv::Symbol,
	::Tuple)::Tuple{Matrix{eltype(epochs)}, Table}

	return epochs, Table((; iv => getproperty(permutation, iv)))
end

@inline function _estimate(::Type{CPEqualVarianceTTest},
	values::SubArray{<:Real},
	permutation::Table,
	iv::Symbol,
	compare::Tuple)

	iv_dat = getproperty(permutation, iv)
	dat_a = @view values[iv_dat .== compare[1]]
	dat_b = @view values[iv_dat .== compare[2]]
	return EqualVarianceTTest(dat_a, dat_b)
end

@inline function _estimate(::Type{CPUnequalVarianceTTest},
	values::SubArray{<:Real},
	permutation::Table,
	iv::Symbol,
	compare::Tuple)

	iv_dat = getproperty(permutation, iv)
	dat_a = @view values[iv_dat .== compare[1]]
	dat_b = @view values[iv_dat .== compare[2]]
	return UnequalVarianceTTest(dat_a, dat_b)
end


####
#### T-test have just a single coefficient
####
time_series_stats(cpt::CPTTest) = time_series_stats(cpt, 1)
cluster_nhd(cpt::CPTTest) = cluster_nhd(cpt, 1)
cluster_ranges(cpt::CPTTest) = cluster_ranges(cpt, 1)
cluster_mass_stats(cpt::CPTTest) = cluster_mass_stats(cpt, 1)
cluster_pvalues(cpt::CPTTest; inhibit_warning::Bool = false) = cluster_pvalues(cpt, 1; inhibit_warning)
#cluster_table(cpt::CPTTest) = cluster_table(cpt, 1)
