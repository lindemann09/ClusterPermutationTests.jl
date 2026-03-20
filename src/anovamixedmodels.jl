"""
    CPAnovaMixedModel

Cluster permutation test using ANOVA F-statistics from a mixed-effects model
(via AnovaMixedModels.jl).

Use `fit(CPAnovaMixedModel, formula, dat, cluster_criterium)` to construct.
"""
struct CPAnovaMixedModel <: CPRegressionModel
	cpc::CPCollection{LinearMixedModel}
	dat::CPData

	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts}
	type::Int
end;

n_threads_default(::CPAnovaMixedModel) = 2

"""
    fit(::Type{<:CPAnovaMixedModel}, f::FormulaTerm, shuffle_ivs, dat::CPData, cluster_criterium;
        mass_fnc=sum, contrasts=Dict(), logger=NullLogger(), type=3)

Fit a cluster permutation test using ANOVA F-statistics from a linear mixed-effects model.

`shuffle_ivs` explicitly specifies which predictor variables are shuffled during permutation.
`type`: ANOVA type (1, 2, or 3; default: 3).
`mass_fnc` is the cluster mass function (default: `sum`).
"""
function StatsAPI.fit(::Type{<:CPAnovaMixedModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}(),
	logger::Union{AbstractLogger, Nothing} = NullLogger(),
	type::Int = 3)

	data, shuffle_ivs = _prepare_regression_data(f, dat, shuffle_ivs)
	cpc = CPCollection{LinearMixedModel}(shuffle_ivs, mass_fnc, cluster_criterium)
	rtn = CPAnovaMixedModel(cpc, data, f, contrasts, type)
	fit_initial_time_series!(rtn; logger)
	return rtn
end

####
#### Parameter estimates
####
function parameter_estimates(cpt::CPAnovaMixedModel,
	design::AbstractStudyDesign,
	time_points::Vector{Int32};
	store_fits::Bool = false)::T2DParamVector

	design = columntable(design)
	param = T2DParamVector()
	md = LinearMixedModel(cpt.f, design; contrasts = cpt.contrasts)
	for t in time_points
		md = refit!(md, view(cpt.dat.epochs, :, t); progress = false, REML = false)
		f = teststat(anova(md, type = cpt.type))
		push!(param, collect(f[2:end]))
		if store_fits
			push!(cpt.cpc.M, md)
		end
	end
	return param
end

time_series_fits(x::CPAnovaMixedModel) = anova.(x.cpc.M, type = x.type)

function StatsAPI.coefnames(cpt::CPAnovaMixedModel)
	rtn = anovatable(anova(first(cpt.cpc.M), type = cpt.type))
	return rtn.rownms[2:end] # remove Intercept
end
