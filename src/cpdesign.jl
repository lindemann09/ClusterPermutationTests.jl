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
		_check_variable(design, v)
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
		elseif !isnothing(unit_obs) && _is_within(v, design, unit_obs)
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
		return _design_type(x.between, x.within)
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
DataFrames.DataFrame(x::PermutationDesign) = design_table(x::PermutationDesign)

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

function cell_indices(x::PermutationDesign; variables::OptMultiSymbolOString = nothing)
	d = design_table(x)
	if isnothing(variables)
		return _cell_indices(d, names(d))
	else
		return _cell_indices(d, _to_string_vector(variables))
	end
end

function HypothesisTests.nobs(x::PermutationDesign)
	ids, combis = cell_indices(x)
	combis[:, "nobs"] = [sum(i) for i in ids]
	return combis
end

shuffle_variable!(perm_design::PermutationDesign, iv::Union{Symbol, String}; kwargs...) =
	shuffle_variable!(Random.GLOBAL_RNG, perm_design, iv; kwargs...)
shuffle_variable!(rng::AbstractRNG, perm_design::PermutationDesign, iv::Symbol;
	kwargs...) = shuffle_variable!(rng, perm_design, String(iv); kwargs...)

function shuffle_variable!(rng::AbstractRNG,
	perm_design::PermutationDesign,
	iv::String;
	synchronize::OptMultiSymbolOString = String[])

	within_vars = perm_design.within_variables
	between_vars = perm_design.between_variables
	iv_is_within = _is_within(iv, between_vars, within_vars) # to be shuffled variable is within (also checks if in design at all)

	design_df = iv_is_within ? perm_design.within : perm_design.between

	# check variables and find all relevant sync variables
	sync = String[] # needed variables
	for s in _to_string_vector(synchronize)
		sync_var_is_within = _is_within(s, between_vars, within_vars)
		if !sync_var_is_within && iv_is_within
			@warn "'$(s)' is a between variable. " *
				  "Between variables arn't affected by the shuffling of a variable ('$(iv)') " *
				  "within the unit of observations."
		elseif sync_var_is_within && !iv_is_within
			@warn "'$(s)' is a within variable. " *
				  "Within variables arn't affected by the shuffling of a property ('$(iv)') " *
				  "of the unit of observations."
		else
			push!(sync, s)
		end
	end
	if isempty(sync)
		# no sync variables: one vector with all true
		shuffle_group_ids = [fill(true, nrow(design_df))]
	else
		shuffle_group_ids, _ = _cell_indices(design_df, sync)
	end

	if !iv_is_within
		# between: shuffle inside cells of synchronized variables
		for i in shuffle_group_ids
			design_df[i, iv] = shuffle(rng, design_df[i, iv])
		end
	else
		# within: consider additionally unit of observations and shuffle inside cells
		for uo in eachcol(perm_design.ids_uo)
			for x in shuffle_group_ids
				i = x .&& uo
				design_df[i, iv] = shuffle(rng, design_df[i, iv])
			end
		end
	end
	return perm_design
end

shuffle_variable(perm_design::PermutationDesign, iv::Union{Symbol, String}; kwargs...) =
	shuffle_variable(Random.GLOBAL_RNG, perm_design, iv; kwargs...)
shuffle_variable(rng::AbstractRNG, perm_design::PermutationDesign, iv::Symbol;
	kwargs...) = shuffle_variable(rng, perm_design, String(iv); kwargs...)

function shuffle_variable(rng::AbstractRNG, perm_design::PermutationDesign, iv::String;
	synchronize::OptMultiSymbolOString = String[])

	pd = copy(perm_design)
	shuffle_variable!(rng, pd, iv; synchronize)
	return pd
end




function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
	println(io, "PermutationDesign ($(x.design_type)) with $(nrow(x)) rows")
	isempty(x.between) ? between = " -- " : (between = join(x.between_variables, ", "))
	isnothing(x.unit_obs) ? unit_obs = " -- " : (unit_obs = x.unit_obs)
	isempty(x.within) ? within = " -- " : (within = join(x.within_variables, ", "))
	println(io, "  between: $(between)")
	println(io, "  within: $(within)")
	println(io, "  unit obs: $(unit_obs)")
	return nothing
end

#
# utilities
#

function _design_type(between::DataFrame, within::DataFrame)
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

function _check_variable(dat::DataFrame, var::SymbolOString)
	return hasproperty(dat, var) || throw(
		ArgumentError("Variable '$var' is not in design table"))
end

function _is_within(variable::SymbolOString, dat::DataFrame, unit_obs::SymbolOString)
	for uo in unique(dat[:, unit_obs])
		i = dat[:, unit_obs] .== uo
		length(unique(skipmissing(dat[i, variable]))) > 1 && return true
	end
	return false
end

function _is_within(variable::String, between_cols::Vector{String}, within_cols::Vector{String})
	# helper function, checks also if variable is in design at all
	if variable in between_cols
		return false
	elseif variable in within_cols
		return true
	else
		throw(ArgumentError("$(variable)' is not a between or within design variable."))
	end
end


_to_string_vector(::Nothing) = String[]
_to_string_vector(x::SymbolOString) = [String(x)]
_to_string_vector(x::Base.AbstractVecOrTuple{Symbol}) = [String(v) for v in x]
_to_string_vector(x::Base.AbstractVecOrTuple{String}) = x

"""get the indices (bool vector) of all combinations of the columns"""
function _cell_indices(dat::DataFrame, columns::Vector{String})
	# find indices (bool vectors) of each unique value in each column and write to dict
	# this intermediate index dict vectors avoids redundant searches in dataframe and speed up code
	unique_vals = [unique(getproperty(dat, col)) for col in columns]
	index_dict = Dict{String, Dict{Any, Vector{Bool}}}() # dict[column][val] = [...indices...]
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
