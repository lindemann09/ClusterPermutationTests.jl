struct ClusterPermutationTTest <: ClusterPermutationTest
	cpc::CPCollection
	data::CPData
end;

function StatsAPI.fit(::Type{ClusterPermutationTTest}, # TODO: two value comparison only, needs to be more general
	iv::SymbolOString,
	dat::CPData;
	equal_variance = true,
	cluster_criterium::ClusterCritODef)

	paired = is_within(iv, dat.design)

	if !paired
        compare = unique(dat.design.between[:, iv])
		throw(ArgumentError("Independent samples t-test not implemented yet. Use paired=true."))
	else
        compare = unique(dat.design.within[:, iv])
	end

	length(compare) == 2 || throw(ArgumentError("'$iv' comprises $(length(compare)) categories; two required."))

	tdef = CPTestDefinition(;
				estimate_fnc=ttest,
				preprocess_fnc=ttest_preprocess,
				mass_fnc=sum,
				# parameter for specs
				paired, iv=Symbol(iv), compare, equal_variance)
	cpc = CPCollection(tdef, cluster_criterium)
	initial_fit!(cpc, dat)
	return ClusterPermutationTTest(cpc, dat)
end;

# formula interface
function StatsAPI.fit(::Type{ClusterPermutationTTest},
	f::FormulaTerm,
	data::CPData;
	kwargs...)

	(f.lhs isa StatsModels.Term && f.rhs isa StatsModels.Term) || throw(
		ArgumentError("Incorrect t.test formula: '$a'"))
	return fit(ClusterPermutationTTest, Symbol(f.rhs), data; kwargs...)
end

function ttest_preprocess(
	mtx::Matrix{<:Real}, design::PermutationDesign, specs::NamedTuple, ### FIXME use CPData?
)::Matrix
	iv = design_table(design)[:, specs.iv]
	if specs[:paired]
		@unpack compare = specs
		a = mtx[iv .== compare[1], :]
		b = mtx[iv .== compare[2], :]
		return b - a # equal size required
	else
		return mtx
	end
end

function ttest(dat::Vector{<:Real}, design::PermutationDesign, specs::NamedTuple)::Float64 ## FIXME use CPData?
	# perform sequential ttests -> parameter
	if specs[:paired]
		tt = OneSampleTTest(dat)
	else
		iv = design_table(design)[:, specs.iv]
		@unpack compare = specs
		dat_a = dat[iv .== compare[1]] # FIXME check
		dat_b = dat[iv .== compare[2]]
		if specs[:equal_variance]
			tt = EqualVarianceTTest(dat_a, dat_b)
		else
			tt = UnequalVarianceTTest(dat_a, dat_b)
		end
	end
	return tt.t
end
