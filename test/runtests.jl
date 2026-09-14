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

	epochs = CSV.File(download(fl_epochs), header = false, ntasks=1)
	dat = CPData(epochs, CSV.read(download(fl_design), Table); unit_obs = :subject_id)

	cl_crit = ClusterCriterium(threshold = 1.69, min_size = 50) # 10%
	cp_config = CPConfig(cl_crit; cluster_statistic = :clusterwise)

	cpt = fit(CPPairedSampleTTest, @formula(y ~ operator_str), dat, cp_config)
	@test npermutations(cpt) == 0
	resample!(cpt, 500; use_threads = false)
	resample!(cpt, 2000; use_threads = true)
	@test length(cluster(cpt)) == 2
	@test cluster_mass_stats(cpt) ≈ [-749.6, -13669.8] atol = 2
	@test cluster_pvalues(cpt) ≈ [0.05, 0.001] atol = 0.01

	cp_config2 = CPConfig(cluster_threshold = 1.69,
				cluster_min_size =50,
				cluster_statistic = :maxmass10)

	cpt = fit(CPPairedSampleTTest, @formula(y ~ operator_str), dat, cp_config2)
	@test cpt.config.cluster_statistic == "maxmass10"
	@test cpt.config.mxms == 10
	resample!(cpt, 1000; use_threads = false)
	resample!(cpt, 2000; use_threads = true)

	@test npermutations(cpt) == 3000
	@test cluster_pvalues(cpt) ≈ [0.017, 0.00] atol = 0.02

	cp_config3 = CPConfig(cluster_threshold = 1.69, cluster_min_size = 50)
	@test cp_config3.cluster_statistic == "maxmass"
	@test cp_config3.mxms == 2
	# FIXME missing test predefined cluster
	cpt_mm = fit(CPMixedModel, @formula(y ~ operator_str + (1|subject_id)), dat,
			cp_config3, reml = true)

	resample!(cpt_mm, 10; use_threads = false)
	summary(cpt_mm)
	@test length(cluster(cpt_mm, 1)) == 2
	@test cluster_mass_stats(cpt_mm, 1) ≈ [749.6, 13669.8] atol = 2

	cpt_amm = fit(CPAnovaMixedModel, @formula(y ~ operator_str + (1|subject_id)), dat,
			cp_config)
	resample!(cpt_amm, 10; use_threads = true)
	summary(cpt_amm)
	@test cluster_mass_stats(cpt_amm, 1) ≈ [124.9, 49656.2] atol = 2
end
