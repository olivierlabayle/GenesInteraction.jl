module TestTMLEInputs

using Test
using CSV
using DataFrames
using TargeneCore
using YAML
using BGEN

function cleanup()
    for file in readdir()
        if startswith(file, "final.")
            rm(file)
        end
    end
end

#####################################################################
##################           UNIT TESTS            ##################
#####################################################################

@testset "Test asb_snps / trans_actors" begin
    bQTLs = TargeneCore.asb_snps(joinpath("data", "asb_files", "asb_"))
    @test bQTLs == ["RSID_17", "RSID_99", "RSID_198"]

    eQTLs = TargeneCore.trans_actors(joinpath("data", "trans_actors_fake.csv"))
    @test eQTLs == ["RSID_102", "RSID_2"]
end

@testset "Test genotypes_encoding" begin
    b = Bgen(BGEN.datadir("example.8bits.bgen"))
    v = variant_by_rsid(b, "RSID_10")
    minor_allele_dosage!(b, v)
    # The minor allele is the first one
    @test minor_allele(v) == alleles(v)[1]
    @test TargeneCore.genotypes_encoding(v) == [2, 1, 0]

    # The minor allele is the second one
    v = variant_by_rsid(b, "RSID_102")
    minor_allele_dosage!(b, v)
    @test minor_allele(v) == alleles(v)[2]
    @test TargeneCore.genotypes_encoding(v) == [0, 1, 2]
end

@testset "Test call_genotypes for a single SNP" begin
    probabilities = [NaN 0.3 0.2 0.9;
                     NaN 0.5 0.2 0.05;
                     NaN 0.2 0.6 0.05]
    variant_genotypes = [2, 1, 0]

    threshold = 0.9
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2] === genotypes[3] === missing
    @test genotypes[4] == 2

    threshold = 0.55
    genotypes = TargeneCore.call_genotypes(
        probabilities, 
        variant_genotypes, 
        threshold)
    @test genotypes[1] === genotypes[2]  === missing
    @test genotypes[3] == 0
    @test genotypes[4] == 2
end

@testset "Test call_genotypes for all SNPs" begin
    snp_list = ["RSID_10", "RSID_100"]
    genotypes = TargeneCore.call_genotypes(joinpath("data", "ukbb", "imputed" , "ukbb"), snp_list, 0.95)
    # I only look at the first 10 rows
    # SAMPLE_ID    
    @test genotypes[1:9, "SAMPLE_ID"] == ["sample_00$i" for i in 1:9]
    # RSID_10
    @test genotypes[1:10, "RSID_10"] == ones(10)
    # RSID_100
    @test all(genotypes[1:10, "RSID_100"] .=== [1, 2, 1, missing, 1, 1, missing, 1, 0, 1])
    # Test column order
    @test DataFrames.names(genotypes) == ["SAMPLE_ID", "RSID_10", "RSID_100"]
end

#####################################################################
###############           END-TO-END TESTS            ###############
#####################################################################

@testset "Test tmle_inputs with-param-files: scenario 1" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - genetic and extra confounders
    # - no covariates
    # - no extra treatments
    # - no batch size
    # - no positivity constraint
    parsed_args = Dict(
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => joinpath("data", "extra_confounders.csv"), 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )
    @test_logs(
        (:warn, "Some treatment variables could not be read from the data files and associated parameter files will not be processed: TREAT_1"), 
        tmle_inputs(parsed_args)
    )

    # Data Files
    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2", "21003", "22001"]
    @test size(confounders) == (490, 5)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 3)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_198"]
    @test Set(unique(treatments[!, "RSID_2"])) == Set([missing, 0, 1, 2])
    @test Set(unique(treatments[!, "RSID_198"])) == Set([0, 1, 2])

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    # Parameter files: untouched because no phenotype batches were specified
    # They are deduplicated for both continuous and binary phenotypes
    @test YAML.load_file(joinpath("config", "param_1.yaml")) == YAML.load_file("final.binary.parameter_1.yaml") == YAML.load_file("final.continuous.parameter_1.yaml")
    @test YAML.load_file(joinpath("config", "param_2.yaml")) == YAML.load_file("final.binary.parameter_2.yaml") == YAML.load_file("final.continuous.parameter_2.yaml")
    
    cleanup()
end

@testset "Test tmle_inputs with-param-files: scenario 2" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - no extra confounders
    # - covariates
    # - SNP and extra treatments
    # - batch size: 1
    # - no positivity constraint
    parsed_args = Dict(
        "with-param-files" => Dict{String, Any}("param-prefix" => joinpath("config", "param_1")), 
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => joinpath("data", "extra_treatments.csv"), 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => nothing, 
        "%COMMAND%" => "with-param-files", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => joinpath("data", "covariates.csv"),
        "phenotype-batch-size" => 1,
        "positivity-constraint" => 0.
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2"]
    @test size(confounders) == (490, 3)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 3)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "TREAT_1"]

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    covariates = CSV.read("final.covariates.csv", DataFrame)
    @test names(covariates) == ["SAMPLE_ID", "COV_1"]
    @test size(covariates) == (490, 2)
    # Parameter files: modified because phenotype batches was specified
    # They are deduplicated for both continuous and binary phenotypes
    origin_1 = YAML.load_file(joinpath("config", "param_1.yaml"))
    binary_1_1 = YAML.load_file("final.binary.parameter_1.yaml")
    binary_1_2 = YAML.load_file("final.binary.parameter_2.yaml")
    continuous_1_1 = YAML.load_file("final.continuous.parameter_1.yaml")
    continuous_1_2 = YAML.load_file("final.continuous.parameter_2.yaml")
    # Parameters and Treatments sections unchanged
    @test binary_1_1["Parameters"] == binary_1_2["Parameters"] == continuous_1_1["Parameters"] == origin_1["Parameters"]
    @test binary_1_1["Treatments"] == binary_1_2["Treatments"] == continuous_1_2["Treatments"] == origin_1["Treatments"]
    # Phenotypes sections changed
    @test binary_1_1["Phenotypes"] == ["BINARY_1"]
    @test binary_1_2["Phenotypes"] == ["BINARY_2"]
    @test continuous_1_1["Phenotypes"] == ["CONTINUOUS_1"]
    @test continuous_1_2["Phenotypes"] == ["CONTINUOUS_2"]

    origin_2 = YAML.load_file(joinpath("config", "param_1_with_extra_treatment.yaml"))
    binary_2_1 = YAML.load_file("final.binary.parameter_3.yaml")
    binary_2_2 = YAML.load_file("final.binary.parameter_4.yaml")
    continuous_2_1 = YAML.load_file("final.continuous.parameter_3.yaml")
    continuous_2_2 = YAML.load_file("final.continuous.parameter_4.yaml")
    # Parameters and Treatments sections unchanged
    @test binary_2_1["Parameters"] == binary_2_2["Parameters"] == continuous_2_1["Parameters"] == origin_2["Parameters"]
    @test binary_2_1["Treatments"] == binary_2_2["Treatments"] == continuous_2_2["Treatments"] == origin_2["Treatments"]
    # Phenotypes sections changed
    @test binary_2_1["Phenotypes"] == ["BINARY_1"]
    @test binary_2_2["Phenotypes"] == ["BINARY_2"]
    @test continuous_2_1["Phenotypes"] == ["CONTINUOUS_1"]
    @test continuous_2_2["Phenotypes"] == ["CONTINUOUS_2"]

    cleanup()
end

@testset "Test tmle_inputs with-asb-trans: scenario 1" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - genetic and extra confounders
    # - no covariates
    # - no extra treatments
    # - no batch size
    # - no param
    # - no positivity constraint
    parsed_args = Dict(
        "with-asb-trans" => Dict{String, Any}(
            "asb-prefix" => joinpath("data", "asb_files", "asb"), 
            "trans-actors" => joinpath("data", "trans_actors_fake.csv"),
            "param-prefix" => nothing
            ),
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => nothing, 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => joinpath("data", "extra_confounders.csv"), 
        "%COMMAND%" => "with-asb-trans", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => nothing,
        "phenotype-batch-size" => nothing,
        "positivity-constraint" => 0.
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2", "21003", "22001"]
    @test size(confounders) == (490, 5)

    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 6)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_102", "RSID_17", "RSID_198", "RSID_99"]

    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]

    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]

    snp_combinations = Set(Iterators.product(["RSID_102", "RSID_2"], ["RSID_17", "RSID_198", "RSID_99"]))
    for index in 1:6
        binary_file = YAML.load_file(string("final.binary.parameter_$index.yaml"))
        continuous_file = YAML.load_file(string("final.continuous.parameter_$index.yaml"))
        @test binary_file == continuous_file
        eqtl, bqtl = Tuple(binary_file["Treatments"])
        setdiff!(snp_combinations, [(eqtl, bqtl)])
        if bqtl == "RSID_198"
            @test size(binary_file["Parameters"], 1) == 5
        else
            @test size(binary_file["Parameters"], 1) == 3
        end
    end
    @test snp_combinations == Set()
    cleanup()
end

@testset "Test tmle_inputs with-asb-trans: scenario 2" begin
    # Scenario:
    # - binary and continuous phenotypes
    # - no extra confounders
    # - covariates
    # - extra treatments
    # - batch size
    # - param
    # - no positivity constraint
    parsed_args = Dict(
        "with-asb-trans" => Dict{String, Any}(
            "asb-prefix" => joinpath("data", "asb_files", "asb"), 
            "trans-actors" => joinpath("data", "trans_actors_fake.csv"),
            "param-prefix" => joinpath("config", "template")
            ),
        "binary-phenotypes" => joinpath("data", "binary_phenotypes.csv"), 
        "call-threshold" => 0.8, 
        "extra-treatments" => joinpath("data", "extra_treatments.csv"), 
        "continuous-phenotypes" => joinpath("data", "continuous_phenotypes.csv"), 
        "extra-confounders" => nothing, 
        "%COMMAND%" => "with-asb-trans", 
        "bgen-prefix" => joinpath("data", "ukbb", "imputed" ,"ukbb"), 
        "genetic-confounders" => joinpath("data", "genetic_confounders.csv"), 
        "out-prefix" => "final", 
        "covariates" => joinpath("data", "covariates.csv"),
        "phenotype-batch-size" => 1,
        "positivity-constraint" => 0.0
    )
    tmle_inputs(parsed_args)

    confounders = CSV.read("final.confounders.csv", DataFrame)
    @test names(confounders) == ["SAMPLE_ID", "PC1", "PC2"]
    @test size(confounders) == (490, 3)
    
    treatments = CSV.read("final.treatments.csv", DataFrame)
    @test size(treatments) == (490, 7)
    @test names(treatments) == ["SAMPLE_ID", "RSID_2", "RSID_102", "RSID_17", "RSID_198", "RSID_99", "TREAT_1"]
    
    binary_phenotypes = CSV.read("final.binary-phenotypes.csv", DataFrame)
    @test size(binary_phenotypes) == (490, 3)
    @test names(binary_phenotypes) == ["SAMPLE_ID", "BINARY_1", "BINARY_2"]
    
    continuous_phenotypes = CSV.read("final.continuous-phenotypes.csv", DataFrame)  
    @test size(continuous_phenotypes) == (490, 3)
    @test names(continuous_phenotypes) == ["SAMPLE_ID", "CONTINUOUS_1", "CONTINUOUS_2"]
    
    covariates = CSV.read("final.covariates.csv", DataFrame)
    @test names(covariates) == ["SAMPLE_ID", "COV_1"]
    @test size(covariates) == (490, 2)
    
    # There are two parameter files, one with extra treatments and one without
    # For each phenotype type (e.g. continuous or binary):
    # We thus expect a maximum of 2 * nb_phenotypes * n_eqtls * n_bqtls = 24 parameter files 
    for index in 1:24
        binary_file = YAML.load_file(string("final.binary.parameter_$index.yaml"))
        continuous_file = YAML.load_file(string("final.continuous.parameter_$index.yaml"))
        @test binary_file["Parameters"] == continuous_file["Parameters"]
        @test binary_file["Treatments"] == continuous_file["Treatments"]
    end
    cleanup()

    # Adding positivity constraint, only 20 files are generated
    parsed_args["positivity-constraint"] = 0.01
    tmle_inputs(parsed_args)

    for index in 1:20
        binary_file = YAML.load_file(string("final.binary.parameter_$index.yaml"))
        continuous_file = YAML.load_file(string("final.continuous.parameter_$index.yaml"))
        @test binary_file["Parameters"] == continuous_file["Parameters"]
        @test binary_file["Treatments"] == continuous_file["Treatments"]
    end
    @test !isfile(string("final.binary.parameter_21.yaml"))

    cleanup()
end




end

true