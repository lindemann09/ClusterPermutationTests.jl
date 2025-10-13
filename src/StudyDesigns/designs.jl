
const EMPTYTABLE = Table((; _ = Vector{Int64}()))

### unit opbs

function UnitObs(name::Symbol, values::CategoricalArray)
	ids_uo = [values .== u for u in unique(values)]
	X = reduce(hcat, ids_uo) # vecvec to matrix, convert to matrix of Bool
	return UnitObs(name, values, X, findfirst.(eachrow(X)))
end

unit_observation(uo::UnitObs) = uo.name
unit_observation(::NoUnitObs) = nothing

### StudyDesign

function StudyDesign(
	design::Table;
	unit_obs::OptSymbolOString = nothing,
	covariates::OptMultiSymbolOString = nothing)

	design_vars = columnnames(design)
	within_names = Symbol[]
	between_names = Symbol[]
	covariate_names = Symbol[]

	# check variables
	if !isnothing(unit_obs)
		unit_obs = Symbol(unit_obs)
		unit_obs in design_vars || _err_not_in_design(unit_obs)
	end
	# covariates
	if !isnothing(covariates)
		append!(covariate_names, to_symbol_vector(covariates))
		for x in covariate_names
			x in design_vars || _err_not_in_design(x)
		end
	end
	# classify variables
	for v in design_vars
		if v == unit_obs || v in covariate_names
			continue
		elseif !_guess_is_categorical(design, v)
			push!(covariate_names, v)
		elseif !isnothing(unit_obs) && _guess_is_within(design, v, unit_obs) # also checks if v in design
			push!(within_names, v)
		else
			push!(between_names, v)
		end
	end

	isempty(between_names) && isempty(within_names) && throw(ArgumentError("No categorical variables found."))

	return make_design(design, unit_obs; between_names, within_names, covariate_names)
end

StudyDesign(design::Any; kwargs...) = StudyDesign(ensure_table(design); kwargs...)

# Helper function to create a `StudyDesign` object from a given experimental design `Table`, specifying which variables are between-subject and which are within-subject.
# converts design variable to categorical
function make_design(design::Table, unit_obs::Symbol;
	between_names::Vector{Symbol},
	within_names::Vector{Symbol},
	covariate_names::Vector{Symbol})

	# unit obs is defined
	uo = UnitObs(unit_obs, categorical(getproperty(design, unit_obs)))

	if !isempty(within_names)
		within = select_col(design, within_names, convert_categorical = true)
	else
		within = nothing
	end

	if !isempty(between_names)
		if unit_obs ∉ between_names
			between_names = vcat(uo.name, between_names)
		end
		between = select_col(design, between_names, convert_categorical = true)
		between = Table(unique(between))
	else
		between = nothing
	end
	if !isempty(within_names)
		within = select_col(design, within_names, convert_categorical = true)
	else
		within = nothing
	end

	co_var_tbl = isempty(covariate_names) ? EMPTYTABLE : select_col(design, covariate_names)

	if !isnothing(between) && !isnothing(within)
		return MixedDesign(between, within, co_var_tbl, uo)
	elseif !isnothing(between)
		return BetweenDesign(between, co_var_tbl, uo)
	elseif !isnothing(within)
		return WithinDesign(within, co_var_tbl, uo)
	else
		throw(ArgumentError("No independent variable found."))
	end
end

function make_design(design::Table,	unit_obs::Nothing;
			between_names::Vector{Symbol},
			covariate_names::Vector{Symbol},
			kwargs...)
	# unit_obs is not defined: It has to be a pure between design
	# each row is a unit of observation: get cell indices of each unique combination of between variables
	between = select_col(design, between_names, convert_categorical = true)
	ids_uo, between = cell_indices(between, between_names)
	X = reduce(hcat, ids_uo) # vecvec to matrix
	i = findfirst.(eachrow(X))
	co_var_tbl = isempty(covariate_names) ? EMPTYTABLE : select_col(design, covariate_names)
	return BetweenDesign(between, co_var_tbl, NoUnitObs(i))
end

unit_observation(d::StudyDesign) = unit_observation(d.uo)
names_between(d::BetweenDesign) = collect(columnnames(d.between))
names_between(d::MixedDesign) = setdiff(columnnames(d.between), [unit_observation(d)])
names_between(::WithinDesign) = Symbol[]
names_within(d::Union{WithinDesign, MixedDesign}) = collect(columnnames(d.within))
names_within(::BetweenDesign) = Symbol[]
names_covariates(d::StudyDesign) = !isempty(d.covariates) ? collect(columnnames(d.covariates)) : Symbol[]
function Base.names(d::StudyDesign)
	uo = unit_observation(d)
	rtn = vcat(names_between(d), names_within(d), names_covariates(d))
	if isnothing(uo)
		return rtn
	else
		return vcat(unit_observation(d), rtn)
	end
end

Base.length(x::StudyDesign) = length(x.uo.i)
Base.copy(x::BetweenDesign) = BetweenDesign(copy(x.between), x.covariates, x.uo) # does not make a copy of covariates and uo, since they will no be modified
Base.copy(x::WithinDesign) = WithinDesign(copy(x.within), x.covariates, x.uo)
Base.copy(x::MixedDesign) = MixedDesign(copy(x.between), copy(x.within), x.covariates, x.uo)

function StatsAPI.nobs(x::StudyDesign)
	ids, combis = cell_indices(x)
	n = (; nobs = [sum(i) for i in ids])
	return Table(combis, n)
end

function Base.show(io::IO, mime::MIME"text/plain", x::StudyDesign)
	println(io, "$(typeof(x)) with $(nrow(x)) rows")
	tmp = unit_observation(x)
	!isnothing(tmp) && println(io, "  unit obs: $(tmp)")
	tmp = names_between(x)
	!isempty(tmp) && println(io, "  between: $(join(tmp, ", "))")
	tmp = names_within(x)
	!isempty(tmp) && println(io, "  within: $(join(tmp, ", "))")
	tmp = names_covariates(x)
	!isempty(tmp) && println(io, "  covariates: $(join(tmp, ", "))")
	return nothing
end

is_covariate(d::StudyDesign, var::Symbol) = var in names_covariates(d)
is_within(::BetweenDesign, var::Symbol) = false
is_between(::WithinDesign, var::Symbol) = false
is_within(d::Union{WithinDesign, MixedDesign}, var::Symbol) = var in names_within(d)
is_between(d::Union{BetweenDesign, MixedDesign}, var::Symbol) = var in names_between(d)
has_variable(d::StudyDesign, var::Symbol) =
	is_between(d, var) || is_within(d, var) || is_covariate(d, var) || var == unit_observation(d)


# utilities / helper

function cell_indices(x::StudyDesign; variables::OptMultiSymbolOString = nothing)
	d = Table(x)
	if isnothing(variables)
		return cell_indices(d, names(x))
	else
		return cell_indices(d, to_symbol_vector(variables))
	end
end


"""get the indices (bool vector) of all combinations of the columns"""
function cell_indices(dat::Table, columns::SymbolVecOrTuple)
	# returns cell indices (bool vectors) and unique combinations (dataframe)

	# find indices (bool vectors) of each unique value in each column and write to dict
	# this intermediate index dict vectors avoids redundant searches in dataframe and speed up code
	unique_vals = [unique(getproperty(dat, col)) for col in columns]
	index_dict = Dict{Symbol, Dict{Any, BitVector}}() # dict[column][val] = [...indices...]
	for (values, col) in zip(unique_vals, columns)
		index_dict[col] = Dict{Any, BitVector}()
		dat_column = getproperty(dat, col)
		for val in values
			index_dict[col][val] = dat_column .== val
		end
	end

	# all combinations of column values => logically combine bool vector of indices
    ids = BitVector[]
	existing_combis = NamedTuple[]
    for uval in Iterators.product(unique_vals...) # loop each possible combi
		x = true
        row = [n=>c for (n, c) in zip(columns, uval)] # naming unique combis,  Vector{Pair{Symbol, String}}
		for (col, val) in row
			x = x .& index_dict[col][val]
		end
		if any(x)
			push!(ids, x)
            push!(existing_combis, NamedTuple(row))
		end
	end

	return ids, Table(existing_combis)
end


function _guess_is_within(dat::Table, var::Symbol, unit_obs::Symbol)
	values = getproperty(dat, var)
	uobs = getproperty(dat, unit_obs)
	for unit in unique(uobs)
		i = uobs .== unit
		length(unique(skipmissing(values[i]))) > 1 && return true
	end
	return false
end

function _guess_is_categorical(dat::Table, var::Symbol)
	T = eltype(getproperty(dat, var))
	return T === Any || T <: Union{CategoricalValue, Missing, Nothing, AbstractString, Bool, Symbol}
end

_err_not_in_design(var::Symbol) = throw(ArgumentError("Variable '$(var)' is not in the design table."))

