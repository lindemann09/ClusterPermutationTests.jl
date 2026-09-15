##
## Definitions of types used internally for sampling and test definitions
##
const TParameterVector = Vector{Float64}
const T2DParameterVector = Vector{TParameterVector}
const TParameterMatrix = Matrix{Float64}

struct ClusterMasses # results of a single permutation
	max_mass::TParameterVector # max mass for each effect
	cluster_mass::T2DParameterVector # cluster-wise masses for each effect (effect X cluster)
end
ClusterMasses(cluster_masses::T2DParameterVector) = ClusterMasses(Float64[], cluster_masses)


###
### CPCollection
###

mutable struct CPCollection{M}
	shuffle_ivs::Vector{Symbol} # name of the to be shuffled independent variable
	Md::Vector{M} # fitted models of initial fit
	coefs::TParameterMatrix # (time X effect) time series statistics of the initial fit
	X::Vector{ClusterMasses} # cluster masses for each permutation (permutation)
end;

function CPCollection{M}(shuffle_ivs::Vector{Symbol}) where {M}
	return CPCollection{M}(shuffle_ivs, M[], zeros(Float64, 0, 0), ClusterMasses[])
end

