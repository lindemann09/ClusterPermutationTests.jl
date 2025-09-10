struct CPTestDefinition
    # estim
    #   signature: function(dat::Vector{<:Real}, design::PermutationDesign, specs::NamedTuple ## FIXME use CPData?
    #   return:    ftype::Type{<:Float64},
    # process
    #   signature: function signature(dat::Matrix{<:Real}, design::PermutationDesign, specs::NamedTuple)
    #   return:    Matrix
    estimate_fnc::Function
    preprocess_fnc::Function
    mass_fnc::Function # required signature function(dat::Vector{<:Real})::Real
    specs::NamedTuple # all parameter passed (together with data) to estimate_fnc
end

function CPTestDefinition(;
    estimate_fnc::Function,
    preprocess_fnc::Union{Nothing,Function},
    mass_fnc::Function=sum,
    kwargs...) # default mass_fnc
    if isnothing(preprocess_fnc)
        preprocess_fnc = _nopreprocessing
    end

    # check return type
    ArgsPreprocessFnc = (Matrix{<:Real}, PermutationDesign, NamedTuple)
    ArgsEstimFnc = (Vector{<:Real}, PermutationDesign, NamedTuple)
    _check_function(estimate_fnc, ArgsEstimFnc, Real, "Parameter estimation")
    _check_function(preprocess_fnc, ArgsPreprocessFnc, Matrix, "Data preprocess")
    specs = (; kwargs...)
    return CPTestDefinition(estimate_fnc, preprocess_fnc, mass_fnc, specs)
end

function test_info(x::CPTestDefinition)
    specs_str = join(["$k=$(v)" for (k,v) in pairs(x.specs)], ", ")
    return "mass_fnc = $(x.mass_fnc), $specs_str"
end;


# utilities

_nopreprocessing(data::Matrix{<:Real}, ::PermutationDesign, ::NamedTuple) = data

function _check_function(fnc::Function, arg_types, rtn_type, func_label::String)
    x = methods(fnc, arg_types)
    if length(x) < 1
        x = methods(fnc)
        throw(
            ArgumentError(
                "$func_label function has to have the " *
                "following positional arguments $arg_types. It is defined like this:\n$x",
            ),
        )
    end

    rtn_types = Base.return_types(fnc)
    for rt in rtn_types
        if !(rt <: rtn_type)
            throw(
                ArgumentError(
                    "$func_label function has to return always a real number " *
                    "and not $rtn_types"),
            )
        end
    end
end
