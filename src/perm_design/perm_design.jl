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

unit_obs(uo::UnitObs) = uo.name
unit_obs(::NoUnitObs) = nothing


###
### PermutationDesign struct
###

abstract type PermutationDesign end

struct BetweenDesign{U <: Union{UnitObs, NoUnitObs}} <: PermutationDesign
	between::Table # unique combinations of between variables
	uo::U
end

struct WithinDesign <: PermutationDesign
	within::Table # unique combinations of within variables
	uo::UnitObs
end

struct MixedDesign <: PermutationDesign
	between::Table # unique combinations of between variables
	within::Table # unique combinations of within variables
	uo::UnitObs
end

include("cell_indices.jl")
include("shuffle_variables.jl")

function PermutationDesign(design::Table;
		unit_obs::OptSymbolOString = nothing)

	if isnothing(unit_obs)
		# no unit of observation: pure between design with no repeated measures/unit_obs
		isempty(design) && throw(ArgumentError("No independent variable found."))
		return make_permutation_design(design, collect(columnnames(design)))
	else
		# unit of observation is defined: check which variables are within and between
		unit_obs =  Symbol(unit_obs)
		within_vars = Symbol[]
		between_vars = Symbol[]
		for v in columnnames(design)
			if v == unit_obs
				continue
			elseif !isnothing(unit_obs) && is_within(v, design, unit_obs)
				push!(within_vars, v)
			else
				push!(between_vars, v)
			end
		end
		isnothing(unit_obs) && !isempty(within_vars) && throw(ArgumentError(
			"A 'unit of observation' variable must be specified if 'within' variables are specified."))
		return make_permutation_design(design, between_vars, within_vars, unit_obs)
	end
end

function PermutationDesign(design::Any; kwargs...)
	Tables.istable(design) || throw(ArgumentError("Design must be a Tables.jl compatible table (e.g., DataFrame or TypedTable)."))
	return PermutationDesign(Table(design); kwargs...)
end


# Helper function to create a `PermutationDesign` object from a given experimental design `Table`, specifying which variables are between-subject and which are within-subject.
function make_permutation_design(
		design::Table,
		between_vars::Vector{Symbol},
		within_vars::Vector{Symbol},
		unit_obs::Symbol)

	# unit obs is defined
	uo = UnitObs(unit_obs, categorical(getproperty(design, unit_obs)))

	if !isempty(within_vars)
		within = select_columns(design, within_vars, convert_categorical=true)
	else
		within = nothing
	end

	if !isempty(between_vars)
		if unit_obs ∉ between_vars
			between_vars = vcat(uo.name, between_vars)
		end
		between = select_columns(design, between_vars, convert_categorical=true)
		between = Table(unique(between))
	else
		between = nothing
	end

	if !isnothing(between) && !isnothing(within)
		return MixedDesign(between, within, uo)
	elseif !isnothing(between)
		return BetweenDesign(between, uo)
	elseif !isnothing(within)
		return WithinDesign(within, uo)
	else
		throw(ArgumentError("No independent variable found."))
	end
end

function make_permutation_design(design::Table, between_vars::Vector{Symbol})
	# unit_obs is not defined in a PURE between design
	# each row is a unit of observation: get cell indices of each unique combination of between variables
	between = select_columns(design, between_vars, convert_categorical=true)
	ids_uo, between = cell_indices(between, between_vars)
	X = reduce(hcat, ids_uo) # vecvec to matrix
	i = findfirst.(eachrow(X))
	return BetweenDesign(between, NoUnitObs(i))
end

unit_obs(perm_design::PermutationDesign) = unit_obs(perm_design.uo)
variables_between(perm_design::Union{BetweenDesign, MixedDesign}) =
	setdiff(columnnames(perm_design.between), [unit_obs(perm_design)])
variables_between(::WithinDesign) = Symbol[]
variables_within(perm_design::Union{WithinDesign, MixedDesign}) =
			collect(columnnames(perm_design.within))
variables_within(::BetweenDesign)  = Symbol[]

Base.length(x::PermutationDesign) = length(x.uo.i)
Base.copy(x::BetweenDesign) = BetweenDesign(copy(x.between), x.uo)
Base.copy(x::WithinDesign) = WithinDesign(copy(x.within), x.uo)
Base.copy(x::MixedDesign) = MixedDesign(copy(x.between), copy(x.within), x.uo)

function HypothesisTests.nobs(x::PermutationDesign) ## TODO TEST
	ids, combis = cell_indices(x)
	n = (; nobs = [sum(i) for i in ids])
	return Table(columns(combis), n)
end

function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
	println(io, "$(typeof(x)) with $(length(x)) rows")
	tmp = unit_obs(x)
	uo_str = isnothing(tmp) ? " -- " :  tmp
	tmp = variables_between(x)
	between_str = isempty(tmp) ? " -- " : join(tmp, ", ")
	tmp = variables_within(x)
	within_str = isempty(tmp) ? " -- " : join(tmp, ", ")
	println(io, "  unit obs: $(uo_str)")
	println(io, "  between: $(between_str)")
	println(io, "  within: $(within_str)")
	return nothing
end

function get_variable(x::PermutationDesign, var::Symbol) # get a single variable
	if var in variables_within(x)
		return getproperty(x.within, var)
	else
		return getproperty(x.between, var)
	end
end
get_variable(x::PermutationDesign, var::String) = get_variable(x, Symbol(var))

TypedTables.Table(x::PermutationDesign) = design_table(x)
function design_table(x::PermutationDesign)::Table
	if x isa MixedDesign
		return Table(_expand_between(x), columns(x.within))
	elseif x isa BetweenDesign
		return Table(_expand_between(x))
	else
		return Table((; x.uo.name => x.uo.values), x.within)
	end
end

function is_within(var::Symbol, dat::Table, unit_obs::Symbol)
	## raises error if variable is not in design
    values = getproperty(dat, var)
    uobs = getproperty(dat, unit_obs)
	for unit in unique(uobs)
		i = uobs .== unit
		length(unique(skipmissing(values[i]))) > 1 && return true
	end
	return false
end

function is_within(var::Symbol, perm_design::BetweenDesign)
	## raises error if variable is not a between or within design variable
	if var in variables_between(perm_design)
		return false
	else
		throw(ArgumentError("$(var)' is not a between or within variable."))
	end
end

function is_within(var::Symbol, perm_design::WithinDesign)
	## raises error if variable is not in design
	if var in variables_within(perm_design)
		return true
	else
		throw(ArgumentError("$(var)' is not a between or within variable."))
	end
end

function is_within(var::Symbol, perm_design::MixedDesign)
	## raises error if variable is not in design
	if var in variables_within(perm_design)
		return true
	elseif var in variables_between(perm_design)
		return false
	else
		throw(ArgumentError("$(var)' is not a between or within variable."))
	end
end

function select_rows(nt::NamedTuple,
	idx::Union{BitVector, AbstractVector{Bool}, AbstractVector{Int}})::NamedTuple
	# utility to select rows of a NamedTuple representation df tabular data
	rtn = NamedTuple()
	for (col, val) in pairs(nt)
		rtn = merge(rtn, (; col=> val[idx]))
	end
	return rtn
end

function select_columns(nt::NamedTuple, cols::AbstractVector{Symbol};
    convert_categorical::Bool=false)::NamedTuple
	# utility to select rows of a NamedTuple representation df tabular data
    rtn = NamedTuple()
    for (col, val) in pairs(nt)
        if col in cols
            val = convert_categorical ? categorical(val) : val
            rtn = merge(rtn, (; col => val))
        end
    end
    return rtn
end

select_columns(tbl::Table, cols::AbstractVector{Symbol}; convert_categorical::Bool=false) =
	 Table(select_columns(columns(tbl), cols; convert_categorical))

function _expand_between(perm_design::PermutationDesign)::NamedTuple
	# expand between design to full length with unit of observations
    return select_rows(columns(perm_design.between), perm_design.uo.i)
end

