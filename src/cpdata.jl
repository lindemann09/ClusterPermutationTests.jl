"""
    ?

Data for a cluster permutation analysis
"""
struct CPData{T<:Real}
    mtx::Matrix{T}
    design::PermuteDesign
end

"""
    CPData(; data_mtx::AbstractMatrix{<:Real},
    ...) TODO

Data for a cluster permutation analysis
"""
function CPData(; data_mtx::AbstractMatrix{<:Real},
    design::Any, ## TODO ? anything that can be converted to a table
    unit_obs::OptColumnIndex,
    kwargs...)

    design_tbl = Table(design)
    # check nrows
    size(data_mtx, 1) == length(design_tbl) || throw(
        DimensionMismatch(
            "Matrix and design table must have the same number of rows!"),
    )

    # check variables
    if !(isnothing(unit_obs))
        unit_obs = Symbol(unit_obs)
        hasproperty(design_tbl, unit_obs) || throw(
            ArgumentError(
                "Unit of observation ('$unit_obs') is not in the design tablee."),
        )
    end
    selected_var = []
    if !isnothing(unit_obs)
        check_variable(design_tbl, unit_obs)
        push!(selected_var, unit_obs)
    end
    ivs = keys(kwargs)
    if length(ivs) == 0
        # take all, excpt unit_obs
        ivs = collect(columnnames(design_tbl))
        filter!(!=(unit_obs), ivs)
    end
    for var in ivs
        check_variable(design_tbl, var)
        push!(selected_var, var)
    end

    # select subset with specified conditions
    ids = ones(Bool, length(design_tbl))
    for (values, var) in zip(values(kwargs), ivs)
        if !(values isa DataAPI.All) # select all columns
            type_check_conditons(values)
            ids = ids .& has_values(design_tbl, var, values)
        end
    end
    # select rows in data (copy data)
    dat = data_mtx[ids, :]
    #select!(design, Not(:row))
    return CPData(dat, PermuteDesign(design_tbl[ids]; ivs, unit_obs))
end

unit_obs(x::CPData) = x.design.uo
unit_obs_name(x::CPData) = x.design.uo_name

iv_names(x::CPData) = iv_names(x.design)
data_matrix(x::CPData) = x.mtx
nepoch_samples(x::CPData) = size(x.mtx, 2)
nepochs(x::CPData) = size(x.mtx, 1)

## utilities
function has_values(tbl::Table, col::Symbol, values::Base.AbstractVecOrTuple)
    dat = getproperty(tbl, col)
    return [i in values for i in dat]
end

function check_variable(design::Table, var::Symbol)
    return hasproperty(design, var) || throw(
        ArgumentError("Variable '$var' is not in design table"))
end

function check_variable(::Table, ::Any)
    throw(
        ArgumentError("Variables have to be a Symbol, not " * string(typeof(var)))
    )
end

type_check_conditons(::StringSymbolOReal) = println("nothing")
function type_check_conditons(values::Union{Tuple,AbstractVector})
    for x in values
        typeof(x) <: StringSymbolOReal || throw(
            ArgumentError(
                "Conditions have to be $(string(StringSymbolOReal)) not $(typeof(x))"),
        )
    end
end

function type_check_conditons(x::Any)
    throw(
        ArgumentError(
            "Conditions have to be $(string(StringSymbolOReal)) not $(typeof(x))")
    )
end
