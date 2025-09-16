abstract type CPTTest <: ClusterPermutationTest end
abstract type CPTwoSampleTTest <: CPTTest end

struct CPPairedSampleTTest <: CPTTest
	cpc::CPCollection
	dat::CPData
	iv::Symbol
	compare::Tuple{String, String}
end;

struct CPEqualVarianceTTest <: CPTwoSampleTTest
	cpc::CPCollection
	dat::CPData
	iv::Symbol
	compare::Tuple{String, String}
end;

struct CPUnequalVarianceTTest <: CPTwoSampleTTest
	cpc::CPCollection
	dat::CPData
	iv::Symbol
	compare::Tuple{String, String}
end;

function StatsAPI.fit(T::Type{<:CPTTest}, # TODO: two value comparison only, needs to be more general
	iv::SymbolOString,
	dat::CPData,
	cluster_criterium::ClusterCritODef;
	mass_fnc::Function = sum)

	iv = String(iv)
	paired = is_within(iv, dat.design)
	if T == CPTTest
		# choose test based on design
		T = paired ? CPPairedSampleTTest : CPEqualVarianceTTest
	end

	if paired
		T == CPPairedSampleTTest || throw(ArgumentError(
			"$(iv) is a within-subject variable. A t-test for independent samples ($(T)) is not possible."))
		compare = unique(dat.design.within[:, iv])
	else
		T == CPPairedSampleTTest && throw(ArgumentError(
			"$(iv) is a between-subject variable. $(T) is not possible."))
		compare = unique(dat.design.between[:, iv])
	end

	length(compare) == 2 || throw(ArgumentError(
		"'$iv' comprises $(length(compare)) categories; two required."))

	cpc = CPCollection(cluster_criterium, mass_fnc)

	rtn = T(cpc, dat, Symbol(iv), (compare[1], compare[2]) )
	initial_fit!(rtn)
	return rtn
end;

# formula interface
function StatsAPI.fit(T::Type{<:CPTTest},
	f::FormulaTerm,
	data::CPData,
	cluster_criterium::ClusterCritODef;
	kwargs...)

	(f.lhs isa Term && f.rhs isa Term) || throw(
		ArgumentError("Incorrect t.test formula: '$f'"))
	return fit(T, Symbol(f.rhs), data, cluster_criterium; kwargs...)
end

####
#### CPPairedSampleTTest
####

function prepare_data(cpt::CPPairedSampleTTest,
	mtx::Matrix{<:Real},
	permutation::PermutationDesign)::Tuple{Matrix{eltype(mtx)}, Table}

	iv = get_variable(permutation, cpt.iv)
	tbl = Table((; cpt.iv => iv))
	a = @view mtx[iv .== cpt.compare[1], :]
	b = @view mtx[iv .== cpt.compare[2], :]
	return b - a, tbl # equal size required
end

function estimate(::CPPairedSampleTTest, samples::SubArray{<:Real}, ::Table)::Float64
	tt = OneSampleTTest(samples)
	return tt.t
end


####
#### CPTwoSampleTTest
####

function prepare_data(cpt::CPTwoSampleTTest,
	mtx::Matrix{<:Real},
	permutation::PermutationDesign)::Tuple{Matrix{eltype(mtx)}, Table}

	return mtx, Table((; cpt.iv => get_variable(permutation, cpt.iv)))
end

function estimate(cpt::CPEqualVarianceTTest,
	samples::SubArray{<:Real},
	permutation::Table)::Float64

	(dat_a, dat_b) = _ttest_get_data(cpt, samples, permutation)
	tt = EqualVarianceTTest(dat_a, dat_b)
	return tt.t
end

function estimate(cpt::CPUnequalVarianceTTest,
	samples::SubArray{<:Real},
	permutation::Table)::Float64

	(dat_a, dat_b) = _ttest_get_data(cpt, samples, permutation)
	tt = UnequalVarianceTTest(dat_a, dat_b)
	return tt.t
end

@inline function _ttest_get_data(cpt::CPTTest, samples::SubArray{<:Real}, permutation::Table)
	# perform sequential ttests -> parameter
	iv = getproperty(permutation, cpt.iv)
	dat_a = @view samples[iv .== cpt.compare[1]] # FIXME check
	dat_b = @view samples[iv .== cpt.compare[2]]
	return dat_a, dat_b
end

## TODO improve show for different tests