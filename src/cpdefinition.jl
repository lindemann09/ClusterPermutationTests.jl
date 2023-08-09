struct ClusterPermutationTestDefinition
    # estim
    #   signature: function(dat::Vector{<:Real}, design::CPDesign, specs::NamedTuple
    #   return:    ftype::Type{<:Float64},
    # process
    #   signature: function signature(dat::Matrix{<:Real}, design::CPDesign, specs::NamedTuple)
    #   return:    Matrix
    estimate_fnc::Function
    preprocess_fnc::Function
    mass_fnc::Function # required signature function(dat::Vector{<:Real})::Real
end

function ClusterPermutationTestDefinition(;
    estimate_fnc::Function,
    preprocess_fnc::Union{Nothing,Function},
    mass_fnc::Function=sum) # default mass_fnc
    if isnothing(preprocess_fnc)
        preprocess_fnc = _nopreprocessing
    end
    # check return type
    ArgsPreprocessFnc = (Matrix{<:Real}, CPDesign, NamedTuple)
    ArgsEstimFnc = (Vector{<:Real}, CPDesign, NamedTuple)
    _check_function(estimate_fnc, ArgsEstimFnc, Real, "Parameter estimation")
    _check_function(preprocess_fnc, ArgsPreprocessFnc, Matrix, "Data preprocess")
    return ClusterPermutationTestDefinition(estimate_fnc, preprocess_fnc, mass_fnc)
end

function repr_functions(x::ClusterPermutationTestDefinition)
    return "$(x.estimate_fnc), $(x.preprocess_fnc), $(x.mass_fnc)"
end;

# helper
_nopreprocessing(data::Matrix{<:Real}, ::CPDesign, ::NamedTuple) = data

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
