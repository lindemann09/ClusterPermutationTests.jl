module StudyDesigns

using DataAPI: DataAPI, nrow, ncol
using StatsAPI: StatsAPI, fit, nobs
using Tables
using Tables: getcolumn, columnnames
using CategoricalArrays: CategoricalValue, CategoricalArray, categorical
using TypedTables: Table
using Random: Random, AbstractRNG, shuffle

export AbstractStudyDesign,
	BetweenDesign, WithinDesign, MixedDesign,
	UnitObs, NoUnitObs,
	study_design,
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
	nrow, ncol, nobs

###
### Unit of Observation struct
###

"""
    UnitObs

Stores the unit-of-observation variable name, its values as a `CategoricalArray`, and
auxiliary indexing structures used to speed up permutation.
"""
struct UnitObs
	name::Symbol
	values::CategoricalArray
	# X: bit matrix for selection of unit of observation to speeds up later processing
	# each column corresponds to bitvector for one unit_obs
	X::BitMatrix
	# indices to reconstruct rows in between design (which only stores unique combinations of between variables)
	i::Vector{Int}
end

"""
    NoUnitObs

Placeholder for designs without an explicit unit-of-observation variable.
Each row of the epoch matrix is treated as a separate unit of observation.
"""
struct NoUnitObs
	i::Vector{Int}
end

###
### Study Design
###

"""
Abstract base type for experimental study designs used in `CPData`.

Concrete subtypes are `BetweenDesign`, `WithinDesign`, and `MixedDesign`.
Each subtype implements the Tables.jl interface and can be converted to a `DataFrame` via
`using DataFrames; DataFrame(design)`.
"""
#abstract type AbstractStudyDesign <: Tables.AbstractColumns end
abstract type AbstractStudyDesign end

"""
    BetweenDesign

Study design with only between-subject variables.

Each unit of observation has a single level for every independent variable.
"""
struct BetweenDesign{U <: Union{UnitObs, NoUnitObs}} <: AbstractStudyDesign
	between::Table # unique combinations of between variables
	covariates::Table # covariates (never used for permutations)
	uo::U
end

"""
    WithinDesign

Study design with only within-subject variables.

Each unit of observation has multiple levels (one per condition) for each independent variable.
"""
struct WithinDesign <: AbstractStudyDesign
	within::Table # unique combinations of within variables
	covariates::Table # covariates (never used for permutations)
	uo::UnitObs
end

"""
    MixedDesign

Study design with both between- and within-subject variables.
"""
struct MixedDesign <: AbstractStudyDesign
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
