abstract type CPTTest <: ClusterPermutationTest end

struct CPPairedSampleTTest <: CPTTest
	cpc::CPCollection
	dat::CPData
	specs::NamedTuple
end;

struct CPEqualVarianceTTest <: CPTTest
	cpc::CPCollection
	dat::CPData
	specs::NamedTuple
end;

struct CPUnequalVarianceTTest <: CPTTest
	cpc::CPCollection
	dat::CPData
	specs::NamedTuple
end;

function StatsAPI.fit(T::Type{<:CPTTest}, # TODO: two value comparison only, needs to be more general
	iv::SymbolOString,
	dat::CPData;
	cluster_criterium::ClusterCritODef,
	mass_fnc::Function = sum)

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
	specs = (; iv = Symbol(iv), compare)

	rtn = T(cpc, dat, specs)
	initial_fit!(rtn)
	return rtn
end;

# formula interface
function StatsAPI.fit(T::Type{<:CPTTest},
	f::FormulaTerm,
	data::CPData;
	kwargs...)

	(f.lhs isa StatsModels.Term && f.rhs isa StatsModels.Term) || throw(
		ArgumentError("Incorrect t.test formula: '$a'"))
	return fit(T, Symbol(f.rhs), data; kwargs...)
end

####
#### CPPairedSampleTTest
####

function prepare_data(cpt::CPPairedSampleTTest,
	mtx::Matrix{<:Real},
	permutation::PermutationDesign)::Tuple{Matrix{eltype(mtx)}, Table}

	@unpack specs = cpt
	iv = get_variable(permutation, specs.iv)
	tbl = Table((; specs.iv => iv))
	a = @view mtx[iv .== specs.compare[1], :]
	b = @view mtx[iv .== specs.compare[2], :]
	return b - a, tbl # equal size required
end

function estimate(::CPPairedSampleTTest, samples::SubArray{<:Real}, ::Table)::Float64
	tt = OneSampleTTest(samples)
	return tt.t
end


####
#### CPEqualVarianceTTest, CPUnequalVarianceTTest
####

function prepare_data(cpt::Union{CPEqualVarianceTTest, CPUnequalVarianceTTest},
	mtx::Matrix{<:Real},
	permutation::PermutationDesign)::Tuple{Matrix{eltype(mtx)}, Table}

	return mtx, Table((; cpt.specs.iv => get_variable(permutation, cpt.specs.iv)))
end

function estimate(cpt::CPEqualVarianceTTest,
	samples::SubArray{<:Real},
	permutation::Table)::Float64

	(dat_a, dat_b) = _ttest_get_data(cpt.specs, samples, permutation)
	tt = EqualVarianceTTest(dat_a, dat_b)
	return tt.t
end

function estimate(cpt::CPUnequalVarianceTTest,
	samples::SubArray{<:Real},
	permutation::Table)::Float64

	(dat_a, dat_b) = _ttest_get_data(cpt.specs, samples, permutation)
	tt = UnequalVarianceTTest(dat_a, dat_b)
	return tt.t
end

@inline function _ttest_get_data(specs::NamedTuple, samples::SubArray{<:Real}, permutation::Table)
	# perform sequential ttests -> parameter
	iv = getproperty(permutation, specs.iv)
	dat_a = @view samples[iv .== specs.compare[1]] # FIXME check
	dat_b = @view samples[iv .== specs.compare[2]]
	return dat_a, dat_b
end