"""
	?

Data for a cluster permutation analysis
"""
struct CPData{T <: Real}
	epochs::Matrix{T}
	design::StudyDesign

	function CPData(epochs::Matrix{T}, design::StudyDesign) where {T <: Real}
		size(epochs, 1) == length(design) || throw(
			DimensionMismatch(
				"Matrix and design table must have the same number of rows!"),
		)
		return new{T}(epochs, design)
	end
end



"""
	CPData(epochs::AbstractMatrix{<:Real},
	...) TODO

Data for a cluster permutation analysis
"""
function CPData(epochs::AbstractMatrix{<:Real},
	design::Table;
	unit_obs::OptSymbolOString,
	kwargs...)

	vars = keys(kwargs)
	if length(vars) == 0
		# take all, except unit_obs
		return CPData(epochs, StudyDesign(design; unit_obs))
	else
		if !isnothing(unit_obs)
			unit_obs = Symbol(unit_obs)
			if unit_obs ∉ vars
				vars = vcat(unit_obs, vars...)
			end
		end
		tbl = select_col(design, vars)
		perm_design = StudyDesign(tbl; unit_obs) # select variables
		return select_epochs(CPData(epochs, perm_design); kwargs...)
	end
end

function CPData(cpdat::CPData; kwargs...)
	if :unit_obs in keys(kwargs)
		unit_obs = kwargs[:unit_obs]
	else
		unit_obs = unit_obs(cpdat.design.uo)
	end
	CPData(cpdat.epochs, columntable(cpdat.design); unit_obs, kwargs...) # copy with selection
end

CPData(epochs::AbstractMatrix{<:Real}, design::Any; unit_obs::OptSymbolOString, kwargs...) =
	CPData(epochs, StudyDesigns.ensure_table(design); unit_obs, kwargs...)


function select_epochs(dat::CPData; kwargs...)

	ivs = keys(kwargs)
	length(ivs) > 0 || throw(ArgumentError("No variables and conditions specified!"))

	dsgn = dat.design
	d_columns = columntable(dsgn)

	# select subset with specified conditions (ivs)
	idx = true
	for (val, col) in zip(values(kwargs), ivs)
		if !(val === :all || val === "all") # select, if not all rows
			idx = idx .& in.(d_columns[col], Ref(val))
		end
	end

	df = Table(select_rows(d_columns, idx)) # selected design table
	perm_design = StudyDesigns.make_design(df, unit_observation(dsgn.uo);
				between_names = names_between(dsgn),
				within_names = names_within(dsgn),
				covariate_names = names_covariates(dsgn))
	return CPData(dat.epochs[idx, :], perm_design)
end

design_table(x::CPData) = Table(x.design)
epoch_length(x::CPData) = size(x.epochs, 2)
nepochs(x::CPData) = size(x.epochs, 1)
StudyDesigns.unit_observation(x::CPData) = unit_observation(x.design.uo)

function Base.show(io::IO, mime::MIME"text/plain", x::CPData)
	println(io, "CPData: matrix $(size(x.epochs))")
	show(io, mime, x.design)
	return nothing
end

