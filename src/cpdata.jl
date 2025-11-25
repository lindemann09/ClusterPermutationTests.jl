"""
	?

Data for a cluster permutation analysis
"""
struct CPData{T <: Real}
	epochs::Matrix{T}
	design::AbstractStudyDesign

	function CPData(epochs::Matrix{T}, design::AbstractStudyDesign) where {T <: Real}
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
		return CPData(epochs, study_design(design; unit_obs))
	else
		if !isnothing(unit_obs)
			unit_obs = Symbol(unit_obs)
			if unit_obs ∉ vars
				vars = vcat(unit_obs, vars...)
			end
		end
		tbl = select_cols(design, vars)
		perm_design = study_design(tbl; unit_obs) # select variables
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

function CPData(epochs::Any, design::Any; unit_obs::OptSymbolOString, kwargs...)

	if design isa AbstractString # try to read file
		design = CSV.File(design, header = true)
	end

	if epochs isa AbstractString # try to read file
		epochs = matrix(CSV.File(epochs, header = false))
	end
	if istable(epochs)
		epochs = matrix(epochs)
	end

	istable(design) ||
		throw(ArgumentError("Design must be a Tables.jl compatible table (e.g., DataFrame or TypedTable)."))
	epochs isa AbstractMatrix{<:Real} ||
		throw(ArgumentError("Epochs must be a matrix or Tables.jl compatible table of real values."))
	CPData(epochs, Table(design); unit_obs, kwargs...)
end


"""
	convert_to_cpdata(data::Any; unit_obs::Union{Symbol, AbstractString},
		bin::Union{Symbol, AbstractString}, response::Union{Symbol, AbstractString})

Convert a Tables.jl compatible table (e.g., DataFrame, TypedTable or Arrow.Table) to a CPData object.

The `data` must contain the columns for the design variables, a `bin` variable indicating the
time bin for each observation, and a `response` variable containing the observed values.

The `unit_obs` specifies the column representing the unit of observation.
"""
function convert_to_cpdata(data::Any; unit_obs::SymbolOString, bin::SymbolOString, response::SymbolOString)::CPData

	if data isa AbstractString # try to read file
		data = CSV.File(data, header = true)
	end
	istable(data) || throw(ArgumentError("Data must be a Tables.jl compatible table (e.g., DataFrame or TypedTable)."))
	data = columntable(data)

	design_vars = [x for x in keys(data) if x ∉ [bin, response]] # design variables
	#tbl_data = Table(select_cols(data, data_vars))
	tbl_design = Table(select_cols(data, design_vars))

	resp = getproperty(data, response)
	unique_bins = sort(unique(getproperty(data, bin)))
	bin_ids = groupfind(getproperty(data, bin))

	resp_mtx = Vector{Vector{Float64}}()
	unique_row_ids = Int[]
	@info "Found $(length(unique_bins)) bins: $(join(unique_bins, ", "))"
	for g_ids in groupfind(tbl_design)
		push!(unique_row_ids, first(g_ids))
		row = Float64[]
		for b in unique_bins
			i = intersect(g_ids, bin_ids[b]) # find all indices for bin, b, in these group indices
			if length(i) > 1
				row_dat = join(tbl_design[first(g_ids)], ", ")
				throw(ArgumentError("Two time bin $(b) for condition ($(row_dat))"))
			end
			val = isempty(i) ? NaN : resp[first(i)]
			push!(row, val)
		end
		push!(resp_mtx, row)
	end

	CPData(stack(resp_mtx, dims = 1), tbl_design[unique_row_ids]; unit_obs)
end


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

