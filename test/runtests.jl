using Aqua
using ClusterPermutationTests
using Test

using MixedModelsDatasets: dataset
using CSV
using TypedTables
using Downloads: download

Aqua.test_all(ClusterPermutationTests; ambiguities = false, deps_compat = true)

@testset "StudyDesigns" begin
	tbl = (A = ["A1", "A2", "A3", "A1", "A2", "A3", "A1", "A2", "A3"],
		B = ["B1", "B1", "B1", "B1", "B1", "B1", "B1", "B1", "B1"],
		C = ["C1", "C1", "C1", "C1", "C1", "C1", "C2", "C2", "C2"])
	d = study_design(tbl)
	@test d isa BetweenDesign
	@test length(names_between(d)) == 3
	@test length(names_within(d)) == 0
	@test unit_observation(d) === nothing

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

	d = study_design(dataset("d3"); exclude_columns = [:g, :h])
	@test length(names_between(d)) == 1
	@test length(names_covariates(d)) == 2
end



@testset "ClusterPermutationTests.jl" begin
	fl_design = "https://raw.githubusercontent.com/lindemann09/JuliaDataSets/refs/heads/main/data/cpt1_design.csv"
	fl_epochs = "https://raw.githubusercontent.com/lindemann09/JuliaDataSets/refs/heads/main/data/cpt1_epochs.dat"

	epochs = CSV.Tables.matrix(CSV.File(download(fl_epochs), header = false))
	dat = CPData(epochs, CSV.read(download(fl_design), Table); unit_obs = :subject_id)

	cluster_criterium = ClusterCriterium(threshold = 1.69, min_size = 50) # 10%

	cpt = fit(CPPairedSampleTTest, @formula(y ~ operator_str), dat, cluster_criterium)
	resample!(cpt, 2000; use_threads = false)
	resample!(cpt, 3000; use_threads = true)
	@test length(cluster(cpt)) == 2
	@test cluster_mass_stats(cpt) ≈ [-749.6, -13669.8] atol = 1
	@test cluster_pvalues(cpt)[2] < 0.001
end
