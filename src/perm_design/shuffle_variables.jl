function shuffle_variable!(rng::AbstractRNG,
	perm_design::PermutationDesign,
	iv::String;
	synchronize::OptMultiSymbolOString = nothing)

	within_vars = perm_design.within_variables
	between_vars = perm_design.between_variables
	iv_is_within = is_within(iv, between_vars, within_vars) # to be shuffled variable is within (also checks if in design at all)

	# check variables and find all relevant sync variables
	sync_vars = String[] # needed variables
	for s in _to_string_vector(synchronize)
		sync_var_is_within = is_within(s, between_vars, within_vars)
		if !sync_var_is_within && iv_is_within
			@warn "'$(s)' is a between variable. " *
				  "Between variables arn't affected by the shuffling of a variable ('$(iv)') " *
				  "within the unit of observations."
		elseif sync_var_is_within && !iv_is_within
			@warn "'$(s)' is a within variable. " *
				  "Within variables arn't affected by the shuffling of a property ('$(iv)') " *
				  "of the unit of observations."
		else
			push!(sync_vars, s)
		end
	end
	# prepare shuffle id groups (to improve performance)
	tmp = iv_is_within ? perm_design.within : perm_design.between
	if isempty(sync_vars)
		# no sync variables: one vector with all true
		shuffle_group_ids = [fill(true, nrow(tmp))]
	else
		shuffle_group_ids, _ = cell_indices(tmp, sync_vars)
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
