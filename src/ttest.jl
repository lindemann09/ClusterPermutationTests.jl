struct ClusterPermutationTTest <: ClusterPermutationTest
    def::ClusterPermutationTestDefinition
    cpc::ClusterPermutationCollection
    data::CPData
end;

function fit(::Type{ClusterPermutationTTest},
    iv::ColumnIndex,
    data_mtx::AbstractMatrix{<:Real},
    design::Any;
    unit_obs::OptColumnIndex,
    paired::Bool,
    compare::Union{Nothing,Base.AbstractVecOrTuple{<:StringSymbolOReal}}=nothing,
    equal_variance=true,
    cluster_criteria::ClusterDef)
    if isnothing(compare)
        compare = unique(design[:, iv])
        length(compare) == 2 || throw(
            ArgumentError(
                "'$iv' comprises not two values: '$compare'. Maybe use compare-parameter to specify conditions.",
            ),
        )
    end

    length(compare) == 2 || throw(ArgumentError(
        "compare has to be tuple of two conditions/groups"))
    cpt_def = ClusterPermutationTestDefinition(ttest,
        ttest_preprocess, sum)
    var = Dict(Symbol(iv) => compare)
    data = CPData(; data_mtx, design, unit_obs, var...) # check &  prepare data
    cpc = ClusterPermutationCollection(;
        cluster_criteria,
        # parameter for specs
        paired, iv, compare, equal_variance)
    initial_fit!(cpc; cpt_def, data)
    return ClusterPermutationTTest(cpt_def, cpc, data)
end;

function fit(::Type{ClusterPermutationTTest},
    f::FormulaTerm,
    data_mtx::AbstractMatrix{<:Real},
    design::Any;
    unit_obs::OptColumnIndex,
    paired::Bool,
    kwargs...)
    (f.lhs isa StatsModels.Term && f.rhs isa StatsModels.Term) || throw(
        ArgumentError("Incorrect t.test formula: '$a'"))
    return fit(ClusterPermutationTTest, Symbol(f.rhs), data_mtx, design;
        unit_obs, paired, kwargs...)
end

function ttest_preprocess(
    mtx::Matrix{<:Real}, design::PermuteDesign, specs::NamedTuple
)::Matrix
    iv = design.ivs[1]
    if specs[:paired]
        @unpack compare = specs
        a = mtx[iv .== compare[1], :]
        b = mtx[iv .== compare[2], :]
        return b - a # equal size required
    else
        return mtx
    end
end

function ttest(dat::Vector{<:Real}, design::PermuteDesign, specs::NamedTuple)::Float64
    # perform sequential ttests -> parameter
    iv = design.ivs[1]
    if specs[:paired]
        tt = OneSampleTTest(dat)
    else
        @unpack compare = specs
        dat_a = dat[iv .== compare[1]]
        dat_b = dat[iv .== compare[2]]
        if specs[:equal_variance]
            tt = EqualVarianceTTest(dat_a, dat_b)
        else
            tt = UnequalVarianceTTest(dat_a, dat_b)
        end
    end
    return tt.t
end
