const EMPTYTABLE = Table((; _ = Vector{Int64}()))

###
### Unit of Observation struct
###
struct UnitObs
	name::Symbol
	values::CategoricalArray
	# X: bit matrix for selection of unit of observation to speeds up later processing
	# each column corresponds to bitvector for one unit_obs
	X::BitMatrix
	# indices to reconstruct rows in between design (which only stores unique combinations of between variables)
	i::Vector{Int}
end

struct NoUnitObs
	i::Vector{Int}
end

function UnitObs(name::Symbol, values::CategoricalArray)
	ids_uo = [values .== u for u in unique(values)]
	X = reduce(hcat, ids_uo) # vecvec to matrix, convert to matrix of Bool
	return UnitObs(name, values, X, findfirst.(eachrow(X)))
end

unit_observation(uo::UnitObs) = uo.name
unit_observation(::NoUnitObs) = nothing


###
### StudyDesign struct
###

#abstract type StudyDesign <: Tables.AbstractColumns end
abstract type StudyDesign end

struct BetweenDesign{U <: Union{UnitObs, NoUnitObs}} <: StudyDesign
	between::Table # unique combinations of between variables
	covariates::Table # covariates (never used for permutations)
	uo::U
end

struct WithinDesign <: StudyDesign
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

struct MixedDesign <: StudyDesign
	between::Table # unique combinations of between variables
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

include("cell_indices.jl")
include("shuffle_variables.jl")

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
		within = select_columns(design, within_names, convert_categorical = true)
	else
		within = nothing
	end

	if !isempty(between_names)
		if unit_obs ∉ between_names
			between_names = vcat(uo.name, between_names)
		end
		between = select_columns(design, between_names, convert_categorical = true)
		between = Table(unique(between))
	else
		between = nothing
	end
	if !isempty(within_names)
		within = select_columns(design, within_names, convert_categorical = true)
	else
		within = nothing
	end

	co_var_tbl = isempty(covariate_names) ? EMPTYTABLE : select_columns(design, covariate_names)

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
	between = select_columns(design, between_names, convert_categorical = true)
	ids_uo, between = cell_indices(between, between_names)
	X = reduce(hcat, ids_uo) # vecvec to matrix
	i = findfirst.(eachrow(X))
	co_var_tbl = isempty(covariate_names) ? EMPTYTABLE : select_columns(design, covariate_names)
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

function HypothesisTests.nobs(x::StudyDesign) ## TODO TEST
	ids, combis = cell_indices(x)
	n = (; nobs = [sum(i) for i in ids])
	return Table(combis, n)
end

function Base.show(io::IO, mime::MIME"text/plain", x::StudyDesign)
	println(io, "$(typeof(x)) with $(length(x)) rows")
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


# required Tables.AbstractColumns object methods
# ## Tables interface

function Tables.getcolumn(d::StudyDesign, var::Symbol) # get a single variable
	if is_within(d, var)
		return getproperty(d.within, var)
	elseif is_between(d, var)
		btw_var = getproperty(d.between, var)
		return btw_var[d.uo.i]
	elseif is_covariate(d, var)
		return getproperty(d.covariates, var)
	else
		_err_not_in_design(var)
	end
end

Tables.istable(::Type{<:StudyDesign}) = true
Tables.columnaccess(::Type{<:StudyDesign}) = true
function Tables.columns(d::StudyDesign)::NamedTuple
	cov = isempty(d.covariates) ? (;) : columns(d.covariates)
	if d isa BetweenDesign
		return merge(_expand_between(d), cov)
	elseif d isa MixedDesign
		return merge(_expand_between(d), columns(d.within), cov)
	else
		# within design
		return merge((; d.uo.name => d.uo.values), columns(d.within), cov)
	end
end
Tables.getcolumn(d::StudyDesign, ::Type{T}, col::Int, var::Symbol) where {T} = getcolumn(d, var)
Tables.getcolumn(d::StudyDesign, i::Int) = getcolumn(d, names(d)[i])
Tables.columnnames(d::StudyDesign) = names(d)


# utilities
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

@inline function _expand_between(d::StudyDesign)::NamedTuple
	if d isa WithinDesign
		return (;)
	else
		# expand between design to full length with unit of observations
		return select_rows(columns(d.between), d.uo.i)
	end
end

