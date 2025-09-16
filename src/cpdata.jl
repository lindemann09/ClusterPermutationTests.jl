"""
	?

Data for a cluster permutation analysis
"""
struct CPData{T <: Real}
	mtx::Matrix{T}
	design::PermutationDesign

	function CPData(mtx::Matrix{T}, design::PermutationDesign) where {T <: Real}
		size(mtx, 1) == length(design) || throw(
			DimensionMismatch(
				"Matrix and design table must have the same number of rows!"),
		)
		return new{T}(mtx, design)
	end
end



"""
	CPData(data_mtx::AbstractMatrix{<:Real},
	...) TODO

Data for a cluster permutation analysis
"""
function CPData(data_mtx::AbstractMatrix{<:Real},
	design::Table;
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
		tbl = select_columns(design, vars)
		perm_design = PermutationDesign(tbl; unit_obs) # select variables
		return select_rows(CPData(data_mtx, perm_design); kwargs...)
	end
end
function CPData(cpdat::CPData; kwargs...)
	if :unit_obs in keys(kwargs)
		unit_obs = kwargs[:unit_obs]
	else
		unit_obs = unit_obs(cpdat.design.uo)
	end
	CPData(cpdat.mtx, design_table(cpdat.design); unit_obs, kwargs...) # copy with selection
end

function CPData(data_mtx::AbstractMatrix{<:Real}, design::Any;
	unit_obs::OptSymbolOString,	kwargs...)
	Tables.istable(design) || throw(ArgumentError("Design must be a Tables.jl compatible table (e.g., DataFrame or TypedTable)."))
	return CPData(data_mtx, Table(design); unit_obs, kwargs...)
end

function select_rows(dat::CPData; kwargs...)

	ivs = keys(kwargs)
	length(ivs) > 0 || throw(ArgumentError("No variables and conditions specified!"))

	dsgn = dat.design
	d_columns = columns(design_table(dsgn))

	# select subset with specified conditions (ivs)
	idx = true
	for (val, col) in zip(values(kwargs), ivs)
		if !(val === :all || val === "all") # select, if not all rows
			idx = idx .& in.(d_columns[col], Ref(val))
		end
	end

	df = Table(select_rows(d_columns, idx)) # selected design table
	if dsgn.uo isa NoUnitObs # no within
		perm_design = make_permutation_design(df, variables_between(dsgn))
	else
		perm_design = make_permutation_design(df, variables_between(dsgn),
			variables_within(dsgn), unit_obs(dsgn.uo))
	end
	return CPData(dat.mtx[idx, :], perm_design)
end

design_table(x::CPData) = design_table(x.design)
epoch_length(x::CPData) = size(x.mtx, 2)
nepochs(x::CPData) = size(x.mtx, 1)

function Base.show(io::IO, mime::MIME"text/plain", x::CPData)
	println(io, "CPData: matrix $(size(x.mtx))")
	show(io, mime, x.design)
	return nothing
end

