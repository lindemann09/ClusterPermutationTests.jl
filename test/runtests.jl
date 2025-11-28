using Aqua
using ClusterPermutationTests
using Test

using MixedModelsDatasets: dataset
using CSV
using TypedTables
using Downloads: download
using RData
using CodecBzip2

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

@testset "CPData" begin

 	d  = RData.load(download(
	 "https://github.com/dalejbarr/clusterperm/raw/refs/heads/master/data/kb07bins.rda"))
	dat = convert_to_cpdata(d["kb07bins"]; unit_obs = :SubjID, bin=:bin, response=:TAS);
	@test names(dat.design) ==  [:SubjID, :Speaker, :Precedent, :Load]
	@test names_within(dat.design) == [:Speaker, :Precedent, :Load]
	@test names_between(dat.design) == []
	@test unit_observation(dat.design) == :SubjID
	@test epoch_length(dat) == 35
	@test nepochs(dat) == 448
end

@testset "ClusterPermutationTests" begin
	fl_design = "https://raw.githubusercontent.com/lindemann09/JuliaDataSets/refs/heads/main/data/cpt1_design.csv"
	fl_epochs = "https://raw.githubusercontent.com/lindemann09/JuliaDataSets/refs/heads/main/data/cpt1_epochs.dat"

	epochs = CSV.File(download(fl_epochs), header = false)
	dat = CPData(epochs, CSV.read(download(fl_design), Table); unit_obs = :subject_id)

	cluster_criterium = ClusterCriterium(threshold = 1.69, min_size = 50) # 10%

	cpt = fit(CPPairedSampleTTest, @formula(y ~ operator_str), dat, cluster_criterium)
	resample!(cpt, 2000; use_threads = false)
	resample!(cpt, 3000; use_threads = true)
	@test length(cluster(cpt)) == 2
	@test cluster_mass_stats(cpt) ≈ [-749.6, -13669.8] atol = 1
	@test cluster_pvalues(cpt)[2] < 0.001


	cpt_mm = fit(CPMixedModel, @formula(y ~ operator_str + (1|subject_id)), dat,
			cluster_criterium, reml = true)
	resample!(cpt_mm, 10; use_threads = false)
	summary(cpt_mm)
	@test length(cluster(cpt_mm, 1)) == 2
	@test cluster_mass_stats(cpt_mm, 1) ≈ [749.6, 13669.8] atol = 1

	cpt_amm = fit(CPAnovaMixedModel, @formula(y ~ operator_str + (1|subject_id)), dat,
			cluster_criterium)
	resample!(cpt_amm, 10; use_threads = true)
	summary(cpt_amm)
	@test cluster_mass_stats(cpt_amm, 1) ≈ [124.9, 49656.2] atol = 1
end
