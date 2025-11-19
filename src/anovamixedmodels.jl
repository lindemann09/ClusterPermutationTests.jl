struct CPAnovaMixedModel <: CPRegressionModel
	cpc::CPCollection{AnovaResult}
	dat::CPData

	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
	reml::Bool # use REML estimation
	type::Int
end;

n_threads_default(::CPAnovaMixedModel) = n_threads_default(::CPMixedModel)

function StatsAPI.fit(::Type{<:CPAnovaMixedModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}(),
	logger::Union{AbstractLogger, Nothing} = NullLogger(),
	reml::Bool = false,
	type::Int = 3)

	data, shuffle_ivs = _prepare_regression_data(f, dat, shuffle_ivs)
	cpc = CPCollection{AnovaResult}(shuffle_ivs, mass_fnc, cluster_criterium)
	rtn = CPAnovaMixedModel(cpc, data, f, contrasts, reml, type)
	fit_initial_time_series!(rtn; logger)
	return rtn
end


####
#### Parameter estimates
####
@inline function parameter_estimates(cpt::CPAnovaMixedModel,
	design::AbstractStudyDesign,
	time_points::Vector{Int32};
	store_fits::Bool = false)::T2DParamVector

	design = columntable(design)
	param = T2DParamVector()

	md = LinearMixedModel(cpt.f, design; contrasts = cpt.contrasts)
	for t in time_points
		md = refit!(md, view(cpt.dat.epochs, :, t); progress = false, REML = cpt.reml)
		aov = anova(md, type = cpt.type)
		f = teststat(aov)
		push!(param, f[2:end])
		if store_fits
			push!(cpt.cpc.m, aov)
		end
	end
	return param
end