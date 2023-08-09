const StringSymbolOReal = Union{AbstractString,Real,Symbol}
const SymbolOString = Union{Symbol,AbstractString}
const OptSymbolOString = Union{SymbolOString, Nothing}
const MultiSymbolOString = Union{SymbolOString,Base.AbstractVecOrTuple{SymbolOString}}

struct CPDesign
	ivs::Table
	uo::Table
	uo_ids::Matrix{Bool}
end

function CPDesign(design::Table;
	ivs::MultiSymbolOString,
	uo::OptSymbolOString=nothing)

    # ivs
    if ivs isa SymbolOString
		ivs = [Symbol(ivs)]
	else
		ivs = [Symbol(v) for v in ivs]
	end
    for var in ivs
        check_variable(design, var)
    end

	values = [getproperty(design, v) for v in ivs]
    ivs_tbl = NamedTuple(zip(ivs, values))

    # uo
	if isnothing(uo)
		uo_tbl = (; none=[])
		uo_ids = Matrix{Bool}(undef, 0, 0)
	else
		uo = Symbol(uo)
		check_variable(design, uo)
        uo in ivs && throw(ArgumentError(
			"'$uo' can't be a variable and the unit of observation.")
		)
        uo_values = getproperty(design, uo)
        uo_tbl = NamedTuple{(uo,)}((uo_values,))
        # id matrix
		uo_ids = [getproperty(design, uo) .== u for u in unique(uo_values)]
        uo_ids = reduce(hcat, uo_ids) # vecvec to matrix
	end

	return CPDesign(Table(ivs_tbl), Table(uo_tbl), uo_ids)
end

Base.length(x::CPDesign) = length(x.ivs)

function Base.copy(x::CPDesign)::CPDesign
    return CPDesign(copy(x.ivs), copy(x.uo), copy(x.uo_ids))
end

function design_table(x::CPDesign)
    length(x.uo) == 0 ? x.ivs : Table(x.uo, x.ivs)
end

function unit_obs(x::CPDesign)
    length(x.uo) == 0 ? nothing : first(columnnames(x.uo))
end

function randperm!(rng::AbstractRNG, perm_design::CPDesign;
    permute_independent=true)

    ivnames = columnnames(perm_design.ivs)
    if permute_independent || length(ivnames) == 1
        for n in ivnames
            iv = getproperty(perm_design.ivs, n)
            for i in eachcol(perm_design.uo_ids)
                iv[i] = shuffle(rng, iv[i]) # change inside table
            end
        end
    else
        # permute multiple ivs together via ids
        shuffled_ids = collect(1:length(perm_design))
        for i in eachcol(perm_design.uo_ids)
            shuffled_ids[i] = shuffle(rng, shuffled_ids[i])
        end
        for n in ivnames
            iv = getproperty(perm_design.ivs, n)
            iv[:] .= iv[shuffled_ids]
        end
    end
    return nothing
end

# utilities
function check_variable(design::Table, var::Symbol)
    return hasproperty(design, var) || throw(
        ArgumentError("Variable '$var' is not in design table"))
end
