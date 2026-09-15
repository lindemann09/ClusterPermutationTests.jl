"""
    CPMixedModel

Cluster permutation test using a linear mixed-effects model (via MixedModels.jl).

Use `fit(CPMixedModel, formula, dat, cluster_criterium)` to construct.
"""
struct CPMixedModel <: CPRegressionModel
	config::CPConfig
	cpc::CPCollection{LinearMixedModel}
	dat::CPData

	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
	reml::Bool # use REML estimation
end;

n_threads_default(::CPMixedModel) = 2

"""
    fit(::Type{<:CPMixedModel}, f::FormulaTerm, shuffle_ivs, dat::CPData, cluster_criterium;
        mass_fnc=sum, contrasts=Dict(), logger=NullLogger(), reml=false)

Fit a cluster permutation test using a linear mixed-effects model.

`shuffle_ivs` explicitly specifies which predictor variables are shuffled during permutation.
`reml`: if `true`, use REML estimation (default: `false`, i.e. ML).
`mass_fnc` is the cluster mass function (default: `sum`).
"""
function StatsAPI.fit(::Type{<:CPMixedModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	config::CPConfig;
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}(),
	reml::Bool = false) ::CPMixedModel

	data, shuffle_ivs = _prepare_regression_data(f, dat, shuffle_ivs)
	cpc = CPCollection{LinearMixedModel}(shuffle_ivs)
	rtn = CPMixedModel(config, cpc, data, f, contrasts, reml)
	fit_initial_time_series!(rtn)
	return rtn
end


####
#### Parameter estimates
####
@inline function parameter_estimates(cpt::CPMixedModel,
	design::AbstractStudyDesign,
	time_points::Vector{<:Integer};
	store_model_fits::Bool = false)::T2DParameterVector

	design = columntable(design)
	param = T2DParameterVector()

	md = LinearMixedModel(cpt.f, design; contrasts = cpt.contrasts)
	for t in time_points
		md = refit!(md, view(cpt.dat.epochs, :, t); progress = false, REML = cpt.reml)
		z = coef(md) ./ stderror(md) # parameter: t-value of effect
		push!(param, z[2:end])
		if store_model_fits
			push!(cpt.cpc.Md, md)
		end
	end
	return param
end