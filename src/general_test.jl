struct ClusterPermutationTestGeneral <: ClusterPermutationTest
	def::ClusterPermutationTestDefinition
	cpc::ClusterPermutationCollection
	data::CPData
end

function StatsAPI.fit(def::ClusterPermutationTestDefinition,
	data_mtx::AbstractMatrix{<:Real},
	design::Any; ##FIXME should use permutationDesign
	cluster_criteria::ClusterDef,
	unit_obs::OptSymbolOString,
	ftype::Type{<:Float64} = Float64,
	kwargs...)::ClusterPermutationTestGeneral

    specs = (; kwargs...)
	T = ftype
	cpc = ClusterPermutationCollection(specs, T[], cluster_criteria, Vector{T}[])
	data = CPData(data_mtx, design; unit_obs)
	initial_fit!(cpc; cpt_def = def, data)
	return ClusterPermutationTestGeneral(def, cpc, data)
end
