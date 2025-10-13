abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection
	dat::CPData
	iv::Symbol
	f::FormulaTerm
end;

function StatsAPI.fit(::Type{<:CPLinearModel}, # TODO: two value comparison only, needs to be more general
	f::FormulaTerm,
	iv::Symbol,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum)

	if is_mixedmodel(f)
		throw(ArgumentError("Mixed models are not yet supported."))
	end
	# TODO only pure between design: check is unit obs must be random effect -> mixedModel

	tbl = prepare_design_table(f, dat.design, dv_dtype = eltype(dat.mtx))

	unit_obs = unit_observation(dat.design)
	if unit_obs ∉ columnnames(tbl)
		unit_obs = nothing
	end
	prepared_data = CPData(dat.mtx, tbl; unit_obs)

	rtn = CPLinearModel(CPCollection(cluster_criterium, mass_fnc),
					prepared_data, iv, f)
	return
	initial_fit!(rtn)
	return rtn
end

####
#### definition of parameter_estimates
####
parameter_estimates(cpt::ClusterPermutationTest, dat::CPData) = 	parameter_estimates(cpt, dat.design, ZERO_RANGE)

"""estimates for a specific section in the time series (cluster) for a given permutation"""
function parameter_estimates(cpt::CPLinearModel, design::StudyDesign, range::TClusterRange)::TParameterVector
	dv_name = Symbol(cpt.f.lhs.sym)
	design = Table(design)
	if length(range) == 0 # take entire time series if zero_range
		dat = cpt.dat.mtx
	else
		dat = @view cpt.dat.mtx[:, range]
	end
	rtn = TParameterVector() # TODO would be vector preallocation faster?
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, dv_name)
	for dv in eachcol(dat)
		dv_data[:] = dv
		md = fit(LinearModel, cpt.f, design) ## fit model!
		if isnothing(i)
			i = findfirst(x->x == cpt.iv, coefnames(md))
		end
		push!(rtn, coef(md)[i])
	end
	return rtn
end



###
### Utilities for regression design tables
###

"""select columns from formula and add empty column for dependent variable"""
function prepare_design_table(f::FormulaTerm, design::StudyDesign;
	 dv_dtype::Type=Float64)::Table

	 # dv
	dv = Vector{dv_dtype}(undef, length(design))
	dv_name = Symbol(f.lhs.sym)

	pred = predictors(f)
	for v in pred
		has_variable(design, v) || throw(ArgumentError("Variable '$(v)' not found in design table!"))
	end

	perm_design = select_columns(columns(design), pred)  # select required variables

	return Table(perm_design, (; dv_name => dv)) # add dependent variable column
end
