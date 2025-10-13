module StudyDesigns

using DataAPI: DataAPI, nrow, ncol
using StatsAPI
using CategoricalArrays: CategoricalValue, CategoricalArray, categorical
using Tables
using Tables: getcolumn, columnnames
using TypedTables: Table
using Random: Random, AbstractRNG, shuffle

export StudyDesign, BetweenDesign, WithinDesign, MixedDesign,
    UnitObs, NoUnitObs,
    names,
    names_within,
    names_between,
    names_covariates,
    has_variable,
    is_covariate,
    is_within,
    is_between,
    unit_observation,
    shuffle_variable!,
    shuffle_variable,
    nobs,
	nrows,
	ncols


###
### Unit of Observation struct
###
struct UnitObs
	name::Symbol
	values::CategoricalArray
	# X: bit matrix for selection of unit of observation to speeds up later processing
	# each column corresponds to bitvector for one unit_obs
	X::BitMatrix
	# indices to reconstruct rows in between design (which only stores unique combinations of between variables)
	i::Vector{Int}
end

struct NoUnitObs
	i::Vector{Int}
end

###
### StudyDesign
###

#abstract type StudyDesign <: Tables.AbstractColumns end
abstract type StudyDesign end

struct BetweenDesign{U <: Union{UnitObs, NoUnitObs}} <: StudyDesign
	between::Table # unique combinations of between variables
	covariates::Table # covariates (never used for permutations)
	uo::U
end

struct WithinDesign <: StudyDesign
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

struct MixedDesign <: StudyDesign
	between::Table # unique combinations of between variables
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

include("utilities.jl")
include("designs.jl")
include("tables.jl")
include("shuffle_variables.jl")

end
