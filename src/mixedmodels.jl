struct CPMixedModel <: CPRegressionModel
	cpc::CPCollection{LinearMixedModel}
	dat::CPData

	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
	log_file::Union{IOStream, Nothing} # iostream to log fitting
end;

n_threads_default(::CPMixedModel) = 2

function StatsAPI.fit(::Type{<:CPMixedModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}(),
	log_file::Union{IOStream, Nothing} = nothing)

	tbl = _prepare_design_table(f, dat.design; dv_dtype = eltype(dat.epochs))
	data = CPData(dat.epochs, tbl; unit_obs = unit_observation(dat.design))

	shuffle_ivs = StudyDesigns.to_symbol_vector(shuffle_ivs)
	for v in shuffle_ivs
		has_variable(data.design, v) || throw(ArgumentError("Variable '$(v)' not found in formula or design table!"))
	end

	cpc = CPCollection{LinearMixedModel}(shuffle_ivs, mass_fnc, cluster_criterium)
	rtn = CPMixedModel(cpc, data, f, contrasts, log_file)
	fit_initial_time_series!(rtn)
	return rtn
end


####
#### Parameter estimates
####
@inline function parameter_estimates(cpt::CPMixedModel,
	design::AbstractStudyDesign,
	time_points::Vector{Int32};
	store_fits::Bool = false)::T2DParamVector

	design = columntable(design)
	param = T2DParamVector()
	if cpt.log_file isa IOStream
		logger = SimpleLogger(cpt.log_file)
	else
		logger = NullLogger()
	end
	local md

	md = LinearMixedModel(cpt.f, design; contrasts = cpt.contrasts)
	for t in time_points
		with_logger(logger) do   # TODO improve logging
			md = refit!(md, view(cpt.dat.epochs, :, t); progress = false, REML = true)
		end # logger
		z = coef(md) ./ stderror(md) # parameter: t-value of effect
		push!(param, z[2:end])
		if store_fits
			push!(cpt.cpc.m, md)
		end
	end
	return param
end