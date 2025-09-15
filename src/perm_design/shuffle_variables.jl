function shuffle_variable!(rng::AbstractRNG,
	perm_design::PermutationDesign,
	iv::String;
	synchronize::OptMultiSymbolOString = nothing)

	_is_categorical(perm_design, iv) || throw(ArgumentError(
		"Can't shuffle design; '$iv' is not a categorical variable."))
	within_vars = perm_design.within_variables
	between_vars = perm_design.between_variables
	iv_is_within = is_within(iv, between_vars, within_vars) # to be shuffled variable is within (also checks if in design at all)

	# prepare shuffle id groups (to improve performance)
	sync_vars = String[]
	if !isnothing(synchronize)
		# check and find all relevant sync variables
		for s in _to_string_vector(synchronize)
			sync_var_is_within = is_within(s, between_vars, within_vars)
			if !sync_var_is_within && iv_is_within
				@warn "'$(s)' is a between variable. " *
					"Between variables can't affect the shuffling of a variable ('$(iv)') " *
					"within the unit of observations."
			elseif sync_var_is_within && !iv_is_within
				@warn "'$(s)' is a within variable. " *
					"Within variables can't affect the shuffling of a property ('$(iv)') " *
					"of the unit of observations."
			else
				_is_categorical(perm_design, s) || throw(ArgumentError(
					"Can't shuffle design; sync variable '$s' is not a categorical variable."))
				push!(sync_vars, s)
			end
		end
	end

	tmp_df = iv_is_within ? perm_design.within : perm_design.between
	if isempty(sync_vars)
		# no sync variables: one vector with all true
		shuffle_group_ids = [fill(true, nrow(tmp_df))]
	else
		shuffle_group_ids, _ = cell_indices(tmp_df, sync_vars)
	end

	# shuffle
	if iv_is_within
		# shuffle within, consider unit of observations and shuffle inside cells
		for uo in eachcol(perm_design.X)
			for x in shuffle_group_ids
				i = x .&& uo
				perm_design.within[i, iv] = shuffle(rng, perm_design.within[i, iv])
			end
		end
	else
		# shuffle between inside cells of synchronized variables (ignore unit of observations)
		for i in shuffle_group_ids
			perm_design.between[i, iv] = shuffle(rng, perm_design.between[i, iv])
		end
	end
	return perm_design
end

shuffle_variable!(perm_design::PermutationDesign, iv::Union{Symbol, String}; synchronize::OptMultiSymbolOString = nothing) =
	shuffle_variable!(Random.GLOBAL_RNG, perm_design, iv; synchronize)
shuffle_variable!(rng::AbstractRNG, perm_design::PermutationDesign, iv::Symbol;
	synchronize::OptMultiSymbolOString = nothing) = shuffle_variable!(rng, perm_design, String(iv); synchronize)

function shuffle_variable(rng::AbstractRNG, perm_design::PermutationDesign, iv::String;
	synchronize::OptMultiSymbolOString = nothing)

	pd = copy(perm_design)
	shuffle_variable!(rng, pd, iv; synchronize)
	return pd
end

shuffle_variable(perm_design::PermutationDesign, iv::Union{Symbol, String}; synchronize::OptMultiSymbolOString = nothing) =
	shuffle_variable(Random.GLOBAL_RNG, perm_design, iv; synchronize)
shuffle_variable(rng::AbstractRNG, perm_design::PermutationDesign, iv::Symbol;
	synchronize::OptMultiSymbolOString = nothing) = shuffle_variable(rng, perm_design, String(iv); synchronize)

## utilities
@inline function _is_categorical(x::PermutationDesign, var::String)
	if var in names(x.within)
		v = getproperty(x.within, var)
	elseif var in names(x.between)
		v = getproperty(x.between, var)
	else
		throw(ArgumentError("Variable '$var' not found in design!"))
	end

	return eltype(v) <: Union{Missing, CategoricalValue}
end