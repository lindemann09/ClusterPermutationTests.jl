function shuffle_variable!(rng::AbstractRNG,
	perm_design::StudyDesign,
	iv::Symbol;
	synchronize::OptMultiSymbolOString = nothing)

	iv_is_within = is_within(perm_design, iv) # to be shuffled variable is within (also checks if in design at all)

	# prepare shuffle id groups (to improve performance)
	sync_vars = Symbol[]
	if !isnothing(synchronize)
		# check and find all relevant sync variables
		for s in to_symbol_vector(synchronize)
			is_covariate(perm_design, s) && throw(ArgumentError("Covariates can't be used to synchronize shuffling."))
			sync_var_is_within = is_within(perm_design, s)
			if unit_observation(perm_design) == s
				@warn "Unit of observation variable ('$(s)') is always synchronized during shuffling."
			elseif !sync_var_is_within && iv_is_within
				@warn "'$(s)' is a between variable. " *
					  "Between variables can't affect the shuffling of a variable ('$(iv)') " *
					  "within the unit of observations."
			elseif sync_var_is_within && !iv_is_within
				@warn "'$(s)' is a within variable. " *
					  "Within variables can't affect the shuffling of a property ('$(iv)') " *
					  "of the unit of observations."
			else
				push!(sync_vars, s)
			end
		end
	end

	tmp_df = iv_is_within ? perm_design.within : perm_design.between
	if isempty(sync_vars)
		# no sync variables: one vector with all true
		shuffle_group_ids = [fill(true, length(tmp_df))]
	else
		shuffle_group_ids, _ = cell_indices(tmp_df, sync_vars)
	end

	# shuffle
	if iv_is_within
		dat = getproperty(perm_design.within, iv)
		# shuffle within, consider unit of observations and shuffle inside cells
		for uo in eachcol(perm_design.uo.X)
			for x in shuffle_group_ids
				i = x .&& uo
				dat[i] = shuffle(rng, dat[i])
			end
		end
	else
		is_covariate(perm_design, iv) && throw(ArgumentError("Covariates like '$(iv)' can't be shuffled."))
		dat = getproperty(perm_design.between, iv)
		# shuffle between inside cells of synchronized variables (ignore unit of observations)
		for i in shuffle_group_ids
			dat[i] = shuffle(rng, dat[i])
		end
	end
	return perm_design
end

function shuffle_variable!(rng::AbstractRNG, perm_design::StudyDesign, iv::Vector{Symbol}; kwargs...)
	for v in iv
		shuffle_variable!(rng, perm_design, v; kwargs...)
	end
	return perm_design
end

shuffle_variable!(rng::AbstractRNG, perm_design::StudyDesign, iv::String; kwargs...) =
	shuffle_variable!(rng, perm_design, Symbol(iv); kwargs...)
shuffle_variable!(rng::AbstractRNG, perm_design::StudyDesign, iv::Vector{String}; kwargs...) =
	shuffle_variable!(rng, perm_design, Symbol.(iv); kwargs...)

shuffle_variable!(perm_design::StudyDesign, iv::Union{Symbol, Vector{Symbol}, String, Vector{String}}; kwargs...) =
	shuffle_variable!(Random.GLOBAL_RNG, perm_design, iv; kwargs...)

shuffle_variable(rng::AbstractRNG, perm_design::StudyDesign, iv::Union{Symbol, Vector{Symbol}, String, Vector{String}}; kwargs...) =
	shuffle_variable!(rng, copy(perm_design), iv; kwargs...)
shuffle_variable(perm_design::StudyDesign, iv::Union{Symbol, Vector{Symbol}, String, Vector{String}}; kwargs...) =
	shuffle_variable!(Random.GLOBAL_RNG, copy(perm_design), iv; kwargs...)

