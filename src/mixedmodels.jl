struct CPMixedModel <: CPRegressionModel
	cpc::CPCollection{LinearMixedModel}
	dat::CPData

	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

n_threads_default(::CPMixedModel) = 1 #floor(Int64, Threads.nthreads()/4)

function StatsAPI.fit(::Type{<:CPMixedModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	tbl = _prepare_design_table(f, dat.design, dv_dtype = eltype(dat.epochs))
	data = CPData(dat.epochs, tbl; unit_obs = unit_observation(dat.design))

	shuffle_ivs = StudyDesigns.to_symbol_vector(shuffle_ivs)
	for v in shuffle_ivs
		has_variable(data.design, v) || throw(ArgumentError("Variable '$(v)' not found in formula or design table!"))
	end

	cpc = CPCollection{LinearMixedModel}(shuffle_ivs, mass_fnc, cluster_criterium)
	rtn = CPMixedModel(cpc, data, f, contrasts)
	fit_initial_time_series!(rtn)
	return rtn
end


####
#### Parameter estimates
####
@inline function parameter_estimates(cpt::CPMixedModel,
	design::AbstractStudyDesign;
	fit_cluster_only::Bool = true,
	store_fits::Bool = false)::TVecTimeXParameter

	design = columntable(design)
	param = TVecTimeXParameter()
	if fit_cluster_only
		time_points = cpt.cpc.tp
	else
		time_points = Int32(1):Int32(epoch_length(cpt.dat))
	end

	dv_data = getproperty(design, cpt.f.lhs.sym)
	logger = NullLogger()
	local md
	for t in time_points
		dv_data[:] = cpt.dat.epochs[:, t] # update dependent variable FIXME view?
		with_logger(logger) do   # TODO improve logging
			md = fit(LinearMixedModel, cpt.f, design; contrasts = cpt.contrasts,
				progress = false, REML = true) ## fit model!
		end # logger
		z = coef(md) ./ stderror(md) # parameter: t-value of effect
		push!(param, z[2:end])
		if store_fits
			push!(cpt.cpc.m, md)
		end
	end
	return param
end