###
### Unit of Observation struct
###
struct UnitObs
	name::String
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

function UnitObs(name::String, values::CategoricalArray)
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
	between::DataFrame # unique combinations of between variables
	uo::U
end

struct WithinDesign <: PermutationDesign
	within::DataFrame # unique combinations of within variables
	uo::UnitObs
end

struct MixedDesign <: PermutationDesign
	between::DataFrame # unique combinations of between variables
	within::DataFrame # unique combinations of within variables
	uo::UnitObs
end

include("cell_indices.jl")
include("shuffle_variables.jl")

function PermutationDesign(design::DataFrame;
		unit_obs::OptSymbolOString = nothing)

	if isnothing(unit_obs)
		# no unit of observation: pure between design with no repeated measures/unit_obs
		isempty(design) && throw(ArgumentError("No independent variable found."))
		return make_permutation_design(design, names(design))
	else
		# unit of observation is defined: check which variables are within and between
		unit_obs =  String(unit_obs)
		within_vars = String[]
		between_vars = String[]
		for v in names(design)
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

# Helper function to create a `PermutationDesign` object from a given experimental design `DataFrame`, specifying which variables are between-subject and which are within-subject.
function make_permutation_design(
		design::DataFrame,
		between_vars::Vector{String},
		within_vars::Vector{String},
		unit_obs::String)

	# unit obs is defined
	uo = UnitObs(unit_obs, categorical(design[:, unit_obs]))

	if !isempty(within_vars)
		within = design[:, within_vars]
		_convert_to_categorical!(within)
	else
		within = nothing
	end

	if !isempty(between_vars)
		if unit_obs ∉ between_vars
			between_vars = vcat(uo.name, between_vars)
		end
		between = unique(design[:, between_vars])
		_convert_to_categorical!(between)
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

function make_permutation_design(design::DataFrame, between_vars::Vector{String})
	# unit_obs is not defined in a PURE between design
	# each row is a unit of observation: get cell indices of each unique combination of between variables
	ids_uo, between = cell_indices(design[:, between_vars], between_vars)
	X = reduce(hcat, ids_uo) # vecvec to matrix
	i = findfirst.(eachrow(X))
	_convert_to_categorical!(between)
	return BetweenDesign(between, NoUnitObs(i))
end

unit_obs(perm_design::PermutationDesign) = unit_obs(perm_design.uo)
variables_between(perm_design::Union{BetweenDesign, MixedDesign}) =
	setdiff(names(perm_design.between), [unit_obs(perm_design)])
variables_between(::WithinDesign) = String[]
variables_within(perm_design::Union{WithinDesign, MixedDesign}) =
			names(perm_design.within)
variables_within(::BetweenDesign)  = String[]

DataFrames.nrow(x::PermutationDesign) = length(x.uo.i)

Base.copy(x::BetweenDesign) = BetweenDesign(copy(x.between), x.uo)
Base.copy(x::WithinDesign) = WithinDesign(copy(x.within), x.uo)
Base.copy(x::MixedDesign) = MixedDesign(copy(x.between), copy(x.within), x.uo)

function HypothesisTests.nobs(x::PermutationDesign) ## TODO TEST
	ids, combis = cell_indices(x)
	combis[:, "nobs"] = [sum(i) for i in ids]
	return combis
end

function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
	println(io, "$(typeof(x)) with $(nrow(x)) rows")
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

function get_variable(x::PermutationDesign, var::String) # get a single variable
	if var in names(x.within)
		return x.within[:, var]
	else
		return x.between[x.uo.i, var]
	end
end
get_variable(x::PermutationDesign, var::Symbol) = get_variable(x, String(var))

DataFrames.DataFrame(x::PermutationDesign) = design_table(x::PermutationDesign)
function design_table(x::PermutationDesign)::DataFrame

	if x isa MixedDesign
		return hcat(x.between[x.uo.i, :], x.within)
	elseif x isa BetweenDesign
		return x.between[x.uo.i, :]
	else
		uodf = DataFrame(x.uo.name => x.uo.values)
		rtn = hcat(uodf, x.within)
	end
end


function is_within(var::SymbolOString, dat::DataFrame, unit_obs::SymbolOString)
	## raises error if variable is not in design
	for uo in unique(dat[:, unit_obs])
		i = dat[:, unit_obs] .== uo
		length(unique(skipmissing(dat[i, var]))) > 1 && return true
	end
	return false
end

function is_within(var::String, perm_design::BetweenDesign)
	## raises error if variable is not a between or within design variable
	if var in variables_between(perm_design)
		return false
	else
		throw(ArgumentError("$(var)' is not a between or within variable."))
	end
end

function is_within(var::String, perm_design::WithinDesign)
	## raises error if variable is not in design
	if var in variables_within(perm_design)
		return true
	else
		throw(ArgumentError("$(var)' is not a between or within variable."))
	end
end

function is_within(var::String, perm_design::MixedDesign)
	## raises error if variable is not in design
	if var in variables_within(perm_design)
		return true
	elseif var in variables_between(perm_design)
		return false
	else
		throw(ArgumentError("$(var)' is not a between or within variable."))
	end
end


function _convert_to_categorical!(df::DataFrame)
	#helper
    # converts strings, symbols and bool variables to categorical
	for v in names(df)
		if getproperty(df, v) isa CategoricalArray
			continue
		end
		transform!(df, v => categorical => v)
	end
	return df
end

function _expand_between(perm_design::PermutationDesign)
    i = perm_design.uo.i
    rtn = NamedTuple()
    for (n, c) in pairs(columns(perm_design.between))
        rtn = merge(rtn, (; n=> c[i]))
    end
    Table(rtn)
end

