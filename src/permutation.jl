const StringSymbolOReal = Union{AbstractString,Real,Symbol}
const ColumnIndex = Union{Symbol,AbstractString}
const OptColumnIndex = Union{ColumnIndex,Nothing}
const MultiColumnIndex = Union{ColumnIndex,Base.AbstractVecOrTuple{ColumnIndex}}

struct PermuteDesign #n=nvars, c=n units
    ivs::NamedTuple
    uo::Vector # unique units of observation
    uo_name::Symbol # unit of observation
    unit_ids::Vector{BitVector} # rename uo_ids
end

function PermuteDesign(design::Table;
    ivs::MultiColumnIndex,
    unit_obs::OptColumnIndex)
    if ivs isa ColumnIndex
        names = [Symbol(ivs)]
    else
        names = [Symbol(v) for v in ivs]
    end

    if isnothing(unit_obs)
        uo_name = :none
        uo = []
        uo_ids = [repeat(BitVector([true]), length(design))]
    else
        uo_name = Symbol(unit_obs)
        uo_name in names && throw(
            ArgumentError(
                "'$uo_name' can't be a variable and the unit of observation.")
        )
        uo = unique(getproperty(design, uo_name))
        uo_ids = [getproperty(design, uo_name) .== u for u in uo]
    end

    sel_var = [getproperty(design, v) for v in names]
    return PermuteDesign(NamedTuple(zip(names, sel_var)), uo, uo_name, uo_ids)
end

Base.length(x::PermuteDesign) = length(x.ivs[1])
iv_names(x::PermuteDesign) = keys(x.ivs)

function deepcopy_ivs(design::PermuteDesign)::PermuteDesign
    copy_ivs = NamedTuple(zip(keys(design.ivs), deepcopy(values(design.ivs))))
    return PermuteDesign(copy_ivs, design.uo, design.uo_name, design.unit_ids)
end

function as_table(x::PermuteDesign)::Table
    d = Dict()
    for (n, v) in pairs(x.ivs)
        d[n] = v
    end
    if length(x.uo) > 0
        # rebuild unit_of_observations
        unit_obs = Vector(undef, length(x))
        for (v, id) in zip(x.uo, x.unit_ids)
            unit_obs[id] .= v
        end
        d[x.uo_name] = unit_obs
    end
    return Table(d)
end

function randperm!(rng::AbstractRNG, permvar::PermuteDesign;
    permute_independent=true)
    if permute_independent || length(permvar.ivs) == 1
        for iv in permvar.ivs
            for i in permvar.unit_ids
                iv[i] = shuffle(rng, iv[i])
            end
        end
    else
        # permute together via ids
        shuffled_ids = collect(1:length(permvar))
        for i in permvar.unit_ids
            shuffled_ids[i] = shuffle(rng, shuffled_ids[i])
        end
        for iv in permvar.ivs
            iv[:] .= iv[shuffled_ids]
        end
    end
    return nothing
end
