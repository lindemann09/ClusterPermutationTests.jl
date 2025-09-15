const StringSymbolOReal = Union{AbstractString, Real, Symbol}
const SymbolOString = Union{Symbol, AbstractString}
const OptSymbolOString = Union{SymbolOString, Nothing}
const OptMultiSymbolOString = Union{SymbolOString, Base.AbstractVecOrTuple{SymbolOString}, Nothing}

struct PermutationDesign{S<:Union{Nothing, String}}
	between::DataFrame
	within::DataFrame
	unit_obs::S
	X::BitMatrix # id matrix for the unit of observation; each column corresponds to one value of unit_obs
end

include("cell_indices.jl")
include("shuffle_variables.jl")

function PermutationDesign(design::DataFrame;
		unit_obs::OptSymbolOString = nothing,
		convert_categorical::Bool = true)

	if isnothing(unit_obs)
		# no unit of observation: pure between design (with no repeated measures)
		isempty(design) && throw(ArgumentError("No independent variable found."))
		return make_permutation_design(design, names(design), convert_categorical)
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
		isempty(between_vars) && isempty(within_vars) && throw(ArgumentError("No independent variable found."))
		isnothing(unit_obs) &&	!isempty(within_vars) && throw(ArgumentError(
			"A 'unit of observation' variable must be specified if 'within' variables are specified."))
		return make_permutation_design(design, between_vars, within_vars, unit_obs, convert_categorical)
	end
end

"""
	make_permutation_design(design::DataFrame, between_vars::Vector{String}, within_vars::Vector{String}, unit_obs::SymbolOString)

Create a `PermutationDesign` object from a given experimental design `DataFrame`, specifying which variables are between-subject and which are within-subject.

This function is for internal use only.
"""
function make_permutation_design(design::DataFrame, between_vars::Vector{String}, within_vars::Vector{String},
		unit_obs::SymbolOString, convert_categorical::Bool)

	# unit obs is defined
	unit_obs = String(unit_obs)
	design = transform(design, unit_obs => categorical => unit_obs) # make copy & ensure unit_obs is categorical
	if convert_categorical
		_convert_to_categorical!(design)
	end

	# add to between and within vars, if required
	if !isempty(within_vars)
		within_vars = vcat(unit_obs, within_vars)
	end
	if !isempty(between_vars)
		between_vars = vcat(unit_obs, between_vars)
	end

	### make X
	# X: matrix with indices for unit of observation to speeds up later processing
	# each column corresponds to one unique value of unit_obs
	unit_obs_values = getproperty(design, unit_obs)
	ids_uo = [unit_obs_values .== u for u in unique(unit_obs_values)]
	X = reduce(hcat, ids_uo) # vecvec to matrix, convert to matrix of Bool
	between = unique(design[:, between_vars])
	within = design[:, within_vars]
	return PermutationDesign(between, within, unit_obs, X)
end

function make_permutation_design(design::DataFrame, between_vars::Vector{String}, convert_categorical::Bool)
	# unit_obs is not defined in a PURE between design
	# each row is a unit of observation: get cell indices of each unique combination of between variables
	ids_uo, between = cell_indices(design[:, between_vars], between_vars)
	if convert_categorical
		_convert_to_categorical!(between)
	end
	X = BitMatrix(reduce(hcat, ids_uo)) # vecvec to matrix
	return PermutationDesign(between, DataFrame(), nothing, X)
end


Base.propertynames(::PermutationDesign) = (:between, :within, :unit_obs, :X,
	:between_variables, :within_variables, :type)
function Base.getproperty(x::PermutationDesign, s::Symbol)
	if s === :type
		return _design_type(x)
	elseif s === :between_variables
		n = names(x.between)
		return isnothing(x.unit_obs) ? n : setdiff(n, [x.unit_obs])
	elseif s === :within_variables
		n = names(x.within)
		return isnothing(x.unit_obs) ? n : setdiff(n, [x.unit_obs])
	else
		return getfield(x, s)
	end
end

DataFrames.nrow(x::PermutationDesign) = size(x.X, 1)

function Base.copy(x::PermutationDesign)::PermutationDesign
	return PermutationDesign(copy(x.between), copy(x.within), x.unit_obs, copy(x.X))
end

function HypothesisTests.nobs(x::PermutationDesign)
	ids, combis = cell_indices(x)
	combis[:, "nobs"] = [sum(i) for i in ids]
	return combis
end


function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
	println(io, "PermutationDesign ($(x.type)) with $(nrow(x)) rows")
	isempty(x.between) ? between = " -- " : (between = join(x.between_variables, ", "))
	isnothing(x.unit_obs) ? unit_obs = " -- " : (unit_obs = x.unit_obs)
	isempty(x.within) ? within = " -- " : (within = join(x.within_variables, ", "))
	println(io, "  between: $(between)")
	println(io, "  within: $(within)")
	println(io, "  unit obs: $(unit_obs)")
	return nothing
end

function get_variable(x::PermutationDesign, var::String) # get a single variable
	if var in names(x.within)
		return x.within[:, var]
	else
		i = findfirst.(eachrow(x.X))
		return x.between[i, var]
	end
end
get_variable(x::PermutationDesign, var::Symbol) = get_variable(x, String(var))

DataFrames.DataFrame(x::PermutationDesign) = design_table(x::PermutationDesign)
function design_table(x::PermutationDesign)
	tp = x.type
	if tp == :mixed
		return hcat(_between_design_table(x), x.within[:, Not(x.unit_obs)])
		# alternatively, but slower and probably less robust:
		# return rightjoin(x.between, x.within, on = x.unit_obs, makeunique = true)
	elseif tp == :between
		return _between_design_table(x)
	else
		return x.within
	end
end

function _between_design_table(x::PermutationDesign)
	# multiple measures per unit of observation possible
	# i: for each row (in X) in which column (-> row in between) is the "1"
	# thus, which row of between corresponds to the rows of the design
	i = findfirst.(eachrow(x.X))
	return x.between[i, :]
end

is_within(variable::SymbolOString, perm_design::PermutationDesign) =
	is_within(String(variable), perm_design.between_variables, perm_design.within_variables)

function is_within(variable::SymbolOString, dat::DataFrame, unit_obs::SymbolOString)
	for uo in unique(dat[:, unit_obs])
		i = dat[:, unit_obs] .== uo
		length(unique(skipmissing(dat[i, variable]))) > 1 && return true
	end
	return false
end

function is_within(variable::String, between_cols::Vector{String}, within_cols::Vector{String})
	# helper function, checks also if variable is in design at all
	if variable in between_cols
		return false
	elseif variable in within_cols
		return true
	else
		throw(ArgumentError("$(variable)' is not a between or within design variable."))
	end
end


function _convert_to_categorical!(df::DataFrame)
	#helper
    # converts strings, symbols and bool variables to categorical
	for v in names(df)
		v_type = eltype(getproperty(df, v))
        if v_type <: Union{Missing, AbstractString} || v_type <: Union{Missing, Symbol} || v_type <: Union{Missing, Bool}
            transform!(df, v => categorical => v)
        end
	end
	return df
end


function _design_type(x::PermutationDesign)
	b = !isempty(x.between)
	w = !isempty(x.within)
	if b && w
		return :mixed
	elseif b
		return :between
	elseif w
		return :within
	else
		return :none
	end
end
