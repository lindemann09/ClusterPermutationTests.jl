"""
	?

Data for a cluster permutation analysis
"""
struct CPData{T <: Real}
	mtx::Matrix{T}
	design::CPDesign
end

"""
	CPData(; data_mtx::AbstractMatrix{<:Real},
	...) TODO

Data for a cluster permutation analysis
"""
function CPData(; data_mtx::AbstractMatrix{<:Real},
	design::Any, ## TODO ? anything that can be converted to a table
	unit_obs::OptSymbolOString,
	kwargs...)

	Tables.istable(design) || throw(ArgumentError(
		"design is not a table, but $(typeof(design)).",))

	design_tbl = TypedTables.Table(design)
	nrows = length(design_tbl)

	size(data_mtx, 1) == nrows|| throw(
		DimensionMismatch(
			"Matrix and design table must have the same number of rows!"),
	)

	ivs = keys(kwargs)
	if length(ivs) == 0
		# take all, excpt unit_obs
		ivs = collect(columnnames(design_tbl))
		filter!(!=(unit_obs), ivs)
	end

	# select subset with specified conditions
	ids = ones(Bool, nrows)
	for (values, var) in zip(values(kwargs), ivs)
		if !(values isa DataAPI.All) # select all columns
			ids = ids .& has_values(design_tbl, var, values)
		end
	end

	return CPData(data_mtx[ids, :],
		CPDesign(design_tbl[ids]; ivs, uo = unit_obs))
end

unit_obs(x::CPData) = unit_obs(x.design)
design_table(x::CPData) = design_table(x.design)
data_matrix(x::CPData) = x.mtx
epoch_length(x::CPData) = size(x.mtx, 2)
nepochs(x::CPData) = size(x.mtx, 1)

## utilities
function has_values(tbl::Table, col::Symbol, values::Base.AbstractVecOrTuple)
	check_variable(tbl, col)
	dat = getproperty(tbl, col)
	return [i in values for i in dat]
end
