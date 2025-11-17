using Aqua
using ClusterPermutationTests
using Test

#Aqua.test_all(ClusterPermutationTests; ambiguities = false, deps_compat = false)

@testset "ClusterPermutationTests.jl" begin end

@testset "StudyDesigns" begin
	using MixedModelsDatasets: dataset
	dat = dataset("kb07")
	d = study_design(dat; unit_obs = :subj)
	@test d isa WithinDesign
	@test length(names(d)) == 7
	@test length(names_within(d)) == 4
	@test length(names_between(d)) == 0
	@test length(names_covariates(d)) == 2
	@test unit_observation(d) == :subj
	@test nrow(d) == 1789
	@test ncol(d) == 7
	@test has_variable(d, :subj)
	@test has_variable(d, :not_existing) == false

	d = study_design(dataset("verbagg"); unit_obs = :subj)
	@test d isa MixedDesign
	@test length(names_between(d)) == 1
	@test length(names_within(d)) == 6
	@test length(names_covariates(d)) == 1
	shuffle_variable(d, :item)
	shuffle_variable!(d, :mode)
	shuffle_variable!(d, [:item, :situ, :r2])

    d = study_design(dataset("d3"); exclude_columns=[:g,:h])
	@test d isa BetweenDesign
    @test length(names_between(d)) == 1
	@test length(names_covariates(d)) == 2
end
