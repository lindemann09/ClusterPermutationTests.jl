"""
	?

Data for a cluster permutation analysis
"""
struct CPData{T <: Real}
	mtx::Matrix{T}
	design::PermutationDesign
end

"""
	CPData(; data_mtx::AbstractMatrix{<:Real},
	...) TODO

Data for a cluster permutation analysis
"""
function CPData(; data_mtx::AbstractMatrix{<:Real},
	design::DataFrame,
	unit_obs::OptSymbolOString,
	kwargs...)

	design_tbl = design
	nrows = nrow(design_tbl)

	size(data_mtx, 1) == nrows|| throw(
		DimensionMismatch(
			"Matrix and design table must have the same number of rows!"),
	)

	ivs = keys(kwargs)
	if length(ivs) == 0
		# take all, except unit_obs
		ivs = collect(names(design_tbl))
		filter!(!=(unit_obs), ivs)
	end

	# select subset with specified conditions
	ids = fill(true, nrows)
	for (values, var) in zip(values(kwargs), ivs)
		check_variable(design_tbl, var)
		if !(values isa DataAPI.All) # select all columns
			ids = ids .&& in.(design_tbl[:, var], Ref(values))
		end
	end
	return CPData(data_mtx[ids, :],
		PermutationDesign(design_tbl[ids, [unit_obs, ivs...]]; unit_obs))
end

design_table(x::CPData) = design_table(x.design)
epoch_length(x::CPData) = size(x.mtx, 2)
nepochs(x::CPData) = size(x.mtx, 1)
