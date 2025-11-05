struct CPMixedModel <: CPRegressionModel
	cpc::CPCollection{LinearMixedModel}
	dat::CPData

	f::FormulaTerm
	effect::String # name of effect to be tested
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

n_threads_default(::CPMixedModel) = floor(Int64, Threads.nthreads()/4)

function StatsAPI.fit(::Type{<:CPMixedModel},
	f::FormulaTerm,
	effect::SymbolOString,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	effect = String(effect)
	iv = first(split(effect, ": "))
	tbl = _prepare_design_table(f, dat.design, dv_dtype = eltype(dat.epochs))
	cpc = CPCollection{LinearMixedModel}(iv, mass_fnc, cluster_criterium)
	rtn = CPMixedModel(cpc,
			CPData(dat.epochs, tbl; unit_obs = unit_observation(dat.design)),
			f, effect, contrasts)
	initial_fit!(rtn)
	_check_effects(rtn.cpc.m[1], effect)
	return rtn
end


####
#### Parameter estimates
####

@inline function parameter_estimates(cpt::CPMixedModel, dat::CPData; store_fits::Bool = false)::TParameterVector
	design = columntable(dat.design)
	param = TParameterVector() # TODO would be vector preallocation faster?
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, cpt.f.lhs.sym)
	logger = NullLogger()
	local md
	for dv in eachcol(dat.epochs)
		dv_data[:] = dv
		with_logger(logger) do # FIXME improve logging
			md = fit(LinearMixedModel, cpt.f, design; contrasts = cpt.contrasts,
					progress=false, REML=true) ## fit model!
		end # logger

		if store_fits
			push!(cpt.cpc.m, md)
		else
			#calc z for effect
			if isnothing(i)
				i = _get_coefficient_row(md, cpt.effect)
			end
			z = coef(md)[i] / stderror(md)[i] # parameter: t-value of effect
			push!(param, z)
		end
	end
	return param
end