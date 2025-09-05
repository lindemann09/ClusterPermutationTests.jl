const StringSymbolOReal = Union{AbstractString, Real, Symbol}
const SymbolOString = Union{Symbol, AbstractString}
const OptSymbolOString = Union{SymbolOString, Nothing}
const OptMultiSymbolOString = Union{SymbolOString, Base.AbstractVecOrTuple{SymbolOString}, Nothing}

struct PermutationDesign
	between::DataFrame
	within::DataFrame
	unit_obs::Union{String, Nothing}
	ids_uo::Matrix{Bool}
end

function PermutationDesign(between::DataFrame, within::DataFrame, unit_obs::Union{String, Nothing})

	# check variables
	nrow(between)==0 && nrow(within)==0 && throw(ArgumentError("No independent variable specified/found."))

	if !isnothing(unit_obs)
		if (!isempty(between) && unit_obs ∉ names(between)) ||
		   (!isempty(within) && unit_obs ∉ names(within))
			throw(ArgumentError("'$unit_obs' is not a valid 'unit of observation' variable. " *
								"It needs to be specified in between and within dataframe as an independent variable."))
		end
	end

	#  matrix with indices for unit of observation, speeds up later processing
	if isempty(within)
		ids_uo = fill(true, nrow(between), 1) # all true
	else
        # within or mixed design
        if isnothing(unit_obs)
            throw(ArgumentError("A 'unit of observation' variable must be specified when 'within' variables are specified."))
        end
		uo_values = getproperty(within, unit_obs)
		# id matrix
		ids_uo = [getproperty(within, unit_obs) .== u for u in unique(uo_values)]
		ids_uo = reduce(hcat, ids_uo) # vecvec to matrix
	end
	PermutationDesign(between, within, unit_obs, ids_uo)
end

function PermutationDesign(design::DataFrame,
	between_vars::Vector{String},
	within_vars::Vector{String};
	unit_obs::OptSymbolOString = nothing)

	if !isnothing(unit_obs)
		unit_obs = String(unit_obs)
		if !isempty(within_vars)
			within_vars = vcat([unit_obs], within_vars)
		end
		if !isempty(between_vars)
			between_vars = vcat([unit_obs], between_vars)
		end
	end
	# check variables
	for v in vcat(between_vars, within_vars)
		check_variable(design, v)
	end
    between = unique(design[:, between_vars])
	PermutationDesign(between, design[:, within_vars], unit_obs)
end

function PermutationDesign(design::DataFrame; unit_obs::OptSymbolOString = nothing)

	unit_obs = isnothing(unit_obs) ? nothing : String(unit_obs)
	within_vars = String[]
	between_vars = String[]
	for v in names(design)
		if v == unit_obs
			continue
		elseif !isnothing(unit_obs) && _is_within(design, v, unit_obs)
			push!(within_vars, v)
		else
			push!(between_vars, v)
		end
	end
	return PermutationDesign(design, between_vars, within_vars; unit_obs)
end


Base.propertynames(::PermutationDesign) = (:between, :within, :unit_obs, :ids_uo,
	:between_variables, :within_variables, :design_type)
function Base.getproperty(x::PermutationDesign, s::Symbol)
	if s === :design_type
		return __design_type(x.between, x.within)
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

DataFrames.nrow(x::PermutationDesign) = maximum((nrow(x.between), nrow(x.within)))

function design_table(x::PermutationDesign)
	tp = x.design_type
	if tp == :mixed
        # i: in which column is the "1" for each row (in the within design), and
        # thus, which row of the between design is need and matched to within design
        i = findfirst.(eachrow(x.ids_uo))
        return hcat(x.between[i, :], x.within[:, Not(x.unit_obs)])
        # alternatively, but slower and probably less robust:
        #return rightjoin(x.between, x.within, on = x.unit_obs, makeunique = true)
	elseif tp == :between
		return x.between
	else
		return x.within
	end
end


function Base.copy(x::PermutationDesign)::PermutationDesign
	return PermutationDesign(copy(x.between), copy(x.within), x.unit_obs, copy(x.ids_uo))
end

function cell_indices(x::PermutationDesign; factors::OptMultiSymbolOString = nothing)
	d = design_table(x)
	if !isnothing(factors)
		d = d[:, _to_string_vector(factors)]
	end
	return _ids_column_combinations(d)
end

function HypothesisTests.nobs(x::PermutationDesign)
	ids, combis = cell_indices(x)
	combis[:, "nobs"] = [sum(i) for i in ids]
	return combis
end

function Random.randperm!(rng::AbstractRNG,
	perm_design::PermutationDesign;
	permute_independent = true)

	if perm_design.design_type != :within
		throw(ArgumentError("Permutation design is not within"))
		return nothing
	end

	design = perm_design.within
	ivnames = perm_design.within_variables ## FIXME needs to treat within and between differently
	if permute_independent || length(ivnames) == 1
		for iv in ivnames
			for i in eachcol(perm_design.ids_uo)
				# shuffle inside each nit of observation
				design[i, iv] = shuffle(rng, design[i, iv])
			end
		end
	else
		# permute multiple ivs together (in the same fashion) via ids
		shuffled_ids = collect(1:nrow(design))
		for i in eachcol(perm_design.ids_uo)
			# shuffle ids per unit of observation
			shuffled_ids[i] = shuffle(rng, shuffled_ids[i])
		end
		for iv in ivnames
			design[:, iv] = design[shuffled_ids, iv]
		end
	end
	return nothing
end

Random.randperm!(perm_design::PermutationDesign; kwargs...) = randperm!(Random.GLOBAL_RNG, perm_design; kwargs...)
function Random.randperm(rng::AbstractRNG, perm_design::PermutationDesign; kwargs...)
	rtn = copy(perm_design)
	randperm!(rng, rtn; kwargs...)
	return rtn
end
Random.randperm(perm_design::PermutationDesign; kwargs...) = randperm(Random.GLOBAL_RNG, perm_design; kwargs...)

function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
	println(io, "PermutationDesign ($(x.design_type)) with $(nrow(x)) rows")
	isempty(x.between) ? between = " -- " : (between = join(x.between_variables, ", "))
	isnothing(x.unit_obs) ?	unit_obs = " -- " : (unit_obs = x.unit_obs)
	isempty(x.within) ? within = " -- " : (within = join(x.within_variables, ", "))
	println(io, "  between: $(between)")
	println(io, "  within: $(within)")
    println(io, "  unit obs: $(unit_obs)")
	return nothing
end

#
# utilities
#

function __design_type(between::DataFrame, within::DataFrame)
	if !isempty(between) && !isempty(within)
		return :mixed
	elseif !isempty(between)
		return :between
	elseif !isempty(within)
		return :within
	else
		return :none
	end
end

function check_variable(dat::DataFrame, var::SymbolOString)
	return hasproperty(dat, var) || throw(
		ArgumentError("Variable '$var' is not in design table"))
end

function _is_within(dat::DataFrame, column::SymbolOString, unit_obs::SymbolOString)
	for uo in unique(dat[:, unit_obs])
		i = dat[:, unit_obs] .== uo
		length(unique(skipmissing(dat[i, column]))) > 1 && return true
	end
	return false
end

_to_string_vector(x::SymbolOString) = return [String(x)]
_to_string_vector(x::Base.AbstractVecOrTuple{SymbolOString}) = return [String(v) for v in x]

"""get the indices (bool vector) of all combinations of the columns"""
function _ids_column_combinations(dat::DataFrame, columns::Vector{Symbol})
	# find indices (bool vectors) of each unique value in each column and write to dict
	# this intermediate index dict vectors avoids redundant searches in dataframe and speed up code
	unique_vals = [unique(getproperty(dat, col)) for col in columns]
	index_dict = Dict{Symbol, Dict{Any, Vector{Bool}}}() # dict[column][val] = [...indices...]
	for (values, col) in zip(unique_vals, columns)
		index_dict[col] = Dict{Any, Vector{Bool}}()
		for val in values
			index_dict[col][val] = dat[:, col] .== val
		end
	end

	# all combinations of column values => logically combine bool vector of indices
	all_combis = sort!(DataFrame(Iterators.product(unique_vals...), columns))
	ids = Vector{Bool}[]
	rm_combis = [] # delete this combination, because they don't exists
	for row in eachrow(all_combis)
		x = true
		for (col, val) in zip(columns, row)
			x = x .&& index_dict[col][val]
		end
		if any(x)
			push!(ids, x)
		else
			# combination does not exist
			push!(rm_combis, getfield(row, :rownumber))
		end
	end

	return ids, delete!(all_combis, rm_combis)
end
_ids_column_combinations(dat::DataFrame, columns::Vector{<:AbstractString}) =
	_ids_column_combinations(dat, Symbol.(columns))
_ids_column_combinations(dat::DataFrame) = _ids_column_combinations(dat, names(dat))
