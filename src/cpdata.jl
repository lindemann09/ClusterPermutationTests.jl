"""
	?

Data for a cluster permutation analysis
"""
struct CPData{T <: Real}
	mtx::Matrix{T}
	design::PermutationDesign

	function CPData(mtx::Matrix{T}, design::PermutationDesign; kwargs...) where {T <: Real}
		size(mtx, 1) == nrow(design) || throw(
			DimensionMismatch(
				"Matrix and design table must have the same number of rows!"),
		)
		rtn = new{T}(mtx, design)
		# select subset, if conditions specified
		return length(kwargs) > 0 ? select_rows(rtn; kwargs...) : rtn
	end
end

"""
	CPData(data_mtx::AbstractMatrix{<:Real},
	...) TODO

Data for a cluster permutation analysis
"""
function CPData(data_mtx::AbstractMatrix{<:Real},
	design::DataFrame;
	unit_obs::OptSymbolOString,
	kwargs...)

	vars = keys(kwargs)
	if length(vars) == 0
		# take all, except unit_obs
		return CPData(data_mtx, PermutationDesign(design; unit_obs))
	else
		if !isnothing(unit_obs)
			unit_obs = Symbol(unit_obs)
			unit_obs ∈ vars && throw(ArgumentError("unit_obs variable '$unit_obs' also specified in conditions!"))
			vars = vcat(unit_obs, vars...)
		end
		perm_design = PermutationDesign(design[:, vars]; unit_obs)
		return CPData(data_mtx, perm_design; kwargs...)
	end
end

function select_rows(dat::CPData; kwargs...)

	ivs = String.(keys(kwargs))
	length(ivs) > 0 || throw(ArgumentError("No variables and conditions specified!"))

	dsgn = dat.design
	df = design_table(dsgn)

	# select subset with specified conditions (ivs)
	idx = fill(true, nrow(df))
	for (val, col) in zip(values(kwargs), ivs)
		if !(val === :all || val === "all") # select, if not all rows
			idx = idx .&& in.(df[:, col], Ref(val))
		end
	end

	perm_design = PermutationDesign(df[idx, :], dsgn.between_variables, dsgn.within_variables;
					unit_obs=dsgn.unit_obs)
	return CPData(dat.mtx[idx, :], perm_design)
end


design_table(x::CPData) = design_table(x.design)
epoch_length(x::CPData) = size(x.mtx, 2)
nepochs(x::CPData) = size(x.mtx, 1)


