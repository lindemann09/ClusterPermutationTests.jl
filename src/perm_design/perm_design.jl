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
### PermutationDesign struct
###

#abstract type PermutationDesign <: Tables.AbstractColumns end
abstract type PermutationDesign end

struct BetweenDesign{U <: Union{UnitObs, NoUnitObs}} <: PermutationDesign
	between::Table # unique combinations of between variables
	covariates::Table # covariates (never used for permutations)
	uo::U
end

struct WithinDesign <: PermutationDesign
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

struct MixedDesign <: PermutationDesign
	between::Table # unique combinations of between variables
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

include("cell_indices.jl")
include("shuffle_variables.jl")

function PermutationDesign(design::Table;
	unit_obs::OptSymbolOString = nothing,
	covariates::OptMultiSymbolOString = nothing)

	design_vars = columnnames(design)
	within_names = Symbol[]
	between_names = Symbol[]
	covariates_names = Symbol[]

	# check variables
	if !isnothing(unit_obs)
		unit_obs = Symbol(unit_obs)
		unit_obs in design_vars || _err_not_in_design(unit_obs)
	end
	# covariates
	if !isnothing(covariates)
		append!(covariates_names, to_symbol_vector(covariates))
		for x in covariates_names
			x in design_vars || _err_not_in_design(x)
		end
	end
	# classify variables
	for v in design_vars
		if v == unit_obs || v in covariates_names
			continue
		elseif !_guess_is_categorical(design, v)
			push!(covariates_names, v)
		elseif !isnothing(unit_obs) && _guess_is_within(design, v, unit_obs) # also checks if v in design
			push!(within_names, v)
		else
			push!(between_names, v)
		end
	end

	isempty(between_names) && isempty(within_names) && throw(ArgumentError("No categorical variables found."))

	if isnothing(unit_obs)
		# no unit of observation: pure between design with no repeated measures/unit_obs
		return make_design(design, between_names, covariates_names)
	else
		return make_design(design, between_names, within_names, covariates_names, unit_obs)
	end
end

PermutationDesign(design::Any; kwargs...) = PermutationDesign(ensure_table(design); kwargs...)

# Helper function to create a `PermutationDesign` object from a given experimental design `Table`, specifying which variables are between-subject and which are within-subject.
# converts design variable to categorical
function make_design(
	design::Table,
	between_vars::Vector{Symbol},
	within_vars::Vector{Symbol},
	covariates::Vector{Symbol},
	unit_obs::Symbol)

	# unit obs is defined
	uo = UnitObs(unit_obs, categorical(getproperty(design, unit_obs)))

	if !isempty(within_vars)
		within = select_columns(design, within_vars, convert_categorical = true)
	else
		within = nothing
	end

	if !isempty(between_vars)
		if unit_obs ∉ between_vars
			between_vars = vcat(uo.name, between_vars)
		end
		between = select_columns(design, between_vars, convert_categorical = true)
		between = Table(unique(between))
	else
		between = nothing
	end
	if !isempty(within_vars)
		within = select_columns(design, within_vars, convert_categorical = true)
	else
		within = nothing
	end

	co_var_tbl = isempty(covariates) ? EMPTYTABLE : select_columns(design, covariates)

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

function make_design(
	design::Table,
	between_vars::Vector{Symbol},
	covariates::Vector{Symbol})
	# unit_obs is not defined in a PURE between design
	# each row is a unit of observation: get cell indices of each unique combination of between variables
	between = select_columns(design, between_vars, convert_categorical = true)
	ids_uo, between = cell_indices(between, between_vars)
	X = reduce(hcat, ids_uo) # vecvec to matrix
	i = findfirst.(eachrow(X))

	co_var_tbl = isempty(covariates) ? EMPTYTABLE : select_columns(design, covariates)

	return BetweenDesign(between, co_var_tbl, NoUnitObs(i))
end

unit_observation(pd::PermutationDesign) = unit_observation(pd.uo)
names_between(pd::BetweenDesign) = collect(columnnames(pd.between))
names_between(pd::MixedDesign) = setdiff(columnnames(pd.between), [unit_observation(pd)])
names_between(::WithinDesign) = Symbol[]
names_within(pd::Union{WithinDesign, MixedDesign}) = collect(columnnames(pd.within))
names_within(::BetweenDesign) = Symbol[]
names_covariates(pd::PermutationDesign) = !isempty(pd.covariates) ? collect(columnnames(pd.covariates)) : Symbol[]

Base.length(x::PermutationDesign) = length(x.uo.i)
Base.copy(x::BetweenDesign) = BetweenDesign(copy(x.between), x.covariates, x.uo) # does not make a copy of covariates and uo, since they will no be modified
Base.copy(x::WithinDesign) = WithinDesign(copy(x.within), x.covariates, x.uo)
Base.copy(x::MixedDesign) = MixedDesign(copy(x.between), copy(x.within), x.covariates, x.uo)

function HypothesisTests.nobs(x::PermutationDesign) ## TODO TEST
	ids, combis = cell_indices(x)
	n = (; nobs = [sum(i) for i in ids])
	return Table(columns(combis), n)
end

function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
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

is_covariate(pd::PermutationDesign, var::Symbol) = var in names_covariates(pd)
is_within(::BetweenDesign, var::Symbol) = false
is_between(::WithinDesign, var::Symbol) = false
is_within(pd::Union{WithinDesign, MixedDesign}, var::Symbol) = var in names_within(pd)
is_between(pd::Union{BetweenDesign, MixedDesign}, var::Symbol) = var in names_between(pd)
has_variable(pd::PermutationDesign, var::Symbol) =
	is_between(pd, var) || is_within(pd, var) || is_covariate(pd, var) || var == unit_observation(pd)

TypedTables.Table(pd::PermutationDesign) = design_table(pd)
function design_table(pd::PermutationDesign)::TypedTables.Table
	cov = isempty(pd.covariates) ? (;) : pd.covariates
	between = _expand_between(pd)
	if pd isa BetweenDesign
		return TypedTables.Table(between, cov)
	elseif pd isa MixedDesign
		return TypedTables.Table(between, pd.within, cov)
	else
		# within design
		return TypedTables.Table((; pd.uo.name => pd.uo.values), pd.within, cov)
	end
end


function design_table(pd::PermutationDesign, only_columns::AbstractVector{Symbol})::TypedTables.Table
	if  !isempty(pd.covariates)
		cov = select_columns(columns(pd.covariates), only_columns)
	else
		cov = (;)
	end
	between = select_columns(_expand_between(pd), only_columns)
	if pd isa BetweenDesign
		return TypedTables.Table(between, cov)
	else
		within = select_columns(columns(pd.within), only_columns)
		if pd isa MixedDesign
			return TypedTables.Table(between, within, cov)
		else
			if pd.uo.name in only_columns
				# within design with selected unit obs
				return TypedTables.Table((; pd.uo.name => pd.uo.values), within, cov)
			else
				# within design without selected unit obs
				return TypedTables.Table(within, cov)
			end
		end
	end
end


# ## Tables interface
# Tables.istable(::Type{<:PermutationDesign}) = true
# Tables.columnaccess(::Type{<:PermutationDesign}) = true
# Tables.columns(pd::PermutationDesign) = design_table(pd)
# # required Tables.AbstractColumns object methods
function Tables.getcolumn(pd::PermutationDesign, var::Symbol) # get a single variable
	if is_within(pd, var)
		return getproperty(pd.within, var)
	elseif is_between(pd, var)
		return getproperty(pd.between, var)
	elseif is_covariate(pd, var)
		return getproperty(pd.covariates, var)
	else
		_err_not_in_design(var)
	end
end
Tables.getcolumn(pd::PermutationDesign, var::String) = getproperty(pd, Symbol(var))
Tables.getcolumn(pd::PermutationDesign, ::Type{T}, col::Int, var::Symbol) where {T} = getcolumn(pd, col)
# FIXME Tables.getcolumn(pd::PermutationDesign, i::Int) = design_table(pd)[:, i]
Tables.columnnames(pd::PermutationDesign) = vcat(names_between(pd), names_within(pd), names_covariates(pd))


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

function _expand_between(pd::PermutationDesign)::NamedTuple
	if pd isa WithinDesign
		return (;)
	else
		# expand between design to full length with unit of observations
		return select_rows(columns(pd.between), pd.uo.i)
	end
end

