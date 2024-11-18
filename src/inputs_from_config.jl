pcnames(pcs) = filter(!=("SAMPLE_ID"), names(pcs))

confounders_from_pcs(pcs, extraW::Nothing) = Symbol.(pcnames(pcs))
confounders_from_pcs(pcs, extraW) = Symbol.(unique(vcat(pcnames(pcs), extraW)))

retrieve_variants_list(config::Dict) = vcat((retrieve_variants_list(x) for x in values(config))...)
retrieve_variants_list(variants::AbstractVector) = variants

"""
    treatment_tuples_from_groups(treatments_lists, orders)

Generate treatment combinations for all order in orders. The final list is sorted 
to encourage reuse of nuisance function in downstreatm estimation.
"""
function treatment_tuples_from_groups(treatments_lists, orders)
    treatment_combinations = []
    for order in orders
        for treatments_lists_at_order in Combinatorics.combinations(treatments_lists, order)
            for treatment_comb in Iterators.product(treatments_lists_at_order...)
                push!(treatment_combinations, treatment_comb)
            end
        end
    end
    return sort(treatment_combinations)
end

function try_append_new_estimands!(
    estimands, 
    dataset, 
    estimand_constructor, 
    treatments, 
    outcomes, 
    confounders;
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1
    )
    local Ψ
    try
        Ψ = factorialEstimands(
        estimand_constructor, treatments, outcomes; 
        confounders=confounders, 
        dataset=dataset,
        outcome_extra_covariates=outcome_extra_covariates,
        positivity_constraint=positivity_constraint, 
        verbosity=verbosity-1
    )
    catch e
        if !(e == ArgumentError("No component passed the positivity constraint."))
            throw(e)
        end
    else
        append!(estimands, Ψ)
    end
end

function estimands_from_variants(
    variants,
    dataset, 
    estimand_constructor, 
    outcomes, 
    confounders;
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1
    )
    estimands = []
    for variant in variants
       
        if isempty(extra_treatments)
            treatments = treatments_from_variant(variant, dataset)
        elseif length(extra_treatments) == 1
            treatments = Dict(treatments_from_variant(variant, dataset)..., treatments_from_variant(string(extra_treatments[1]), dataset)...)
        else
            error("GWIS mode only supports pairwise interaction with one extra treatment.")
        end
        
        local Ψ
        try
            Ψ = factorialEstimands(
            estimand_constructor, treatments, outcomes; 
            confounders=confounders, 
            dataset=dataset,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint, 
            verbosity=verbosity-1)

        catch e
            if !(e == ArgumentError("No component passed the positivity constraint."))
                throw(e)
            end
        else
            append!(estimands, Ψ)
        end
    end
    return estimands
end

function estimands_from_groups(estimands_configs, dataset, variants_config, outcomes, confounders; 
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1
    )
    estimands = []
    for (groupname, variants_dict) ∈ variants_config
        verbosity > 0 && @info(string("Generating estimands for group: ", groupname))
        for estimands_config in estimands_configs
            estimand_constructor = eval(Symbol(estimands_config["type"]))
            orders = haskey(estimands_config, "orders") ? estimands_config["orders"] : default_order(estimand_constructor)
            treatments_lists = [Symbol.(variant_list) for variant_list in values(variants_dict)]
            isempty(extra_treatments) || push!(treatments_lists, extra_treatments)
            for treatments ∈ treatment_tuples_from_groups(treatments_lists, orders)
                try_append_new_estimands!(
                    estimands, 
                    dataset, 
                    estimand_constructor, 
                    treatments, 
                    outcomes, 
                    confounders;
                    outcome_extra_covariates=outcome_extra_covariates,
                    positivity_constraint=positivity_constraint,
                    verbosity=verbosity-1
                )
            end
        end
    end
    return estimands
end

default_order(::Union{typeof(ATE), typeof(CM)}) = [1]
default_order(::typeof(AIE)) = [2]

treatment_tuples_from_single_list(treatment_variables, orders) =
    reduce(vcat, collect(Combinatorics.combinations(treatment_variables, order)) for order in orders)

function estimands_from_flat_list(estimands_configs, dataset, variants, outcomes, confounders; 
    extra_treatments=[],
    outcome_extra_covariates=[],
    positivity_constraint=0.,
    verbosity=1,
    )
    estimands = []
    treatment_variables = isempty(extra_treatments) ? Symbol.(variants) : Symbol.(vcat(variants, extra_treatments))
    for estimands_config in estimands_configs
        estimand_constructor = eval(Symbol(estimands_config["type"]))
        verbosity > 0 && @info(string("Generating estimands: ", estimand_constructor))
        orders = haskey(estimands_config, "orders") ? estimands_config["orders"] : default_order(estimand_constructor)
        for treatments ∈ treatment_tuples_from_single_list(treatment_variables, orders)
            try_append_new_estimands!(
                estimands, 
                dataset, 
                estimand_constructor, 
                treatments, 
                outcomes, 
                confounders;
                outcome_extra_covariates=outcome_extra_covariates,
                positivity_constraint=positivity_constraint,
                verbosity=verbosity-1
            )
        end
    end
    return estimands
end

function treatments_from_variant(variant::Symbol, dataset)
    variant_levels = sort(levels(dataset[!, variant], skipmissing=true))
    return NamedTuple{(variant,)}([variant_levels])
end

"""
    treatment_from_variant(variant, dataset)

    Generate a key-value pair (dicitionary) for treatment structs.
"""
function treatments_from_variant(variant::String, dataset::DataFrame)
    variant_levels = sort(levels(dataset[!, variant], skipmissing=true))
    return Dict{Symbol, Vector{UInt8}}(Symbol(variant)=>variant_levels)
end

function estimands_from_gwas(dataset, variants, outcomes, confounders; 
    extra_treatments=extra_treatments,
    outcome_extra_covariates = [],
    positivity_constraint=0.,
    verbosity=0
    )
    variants_groups = Iterators.partition(variants, length(variants) ÷ Threads.nthreads())
    estimands_tasks = map(variants_groups) do variants
        Threads.@spawn estimands_from_variants(variants, dataset, ATE, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    end
    estimands_partitions = fetch.(estimands_tasks)
    return vcat(estimands_partitions...)
end

function estimands_from_gwis(dataset, variants, outcomes, confounders;
    extra_treatments=extra_treatments,
    outcome_extra_covariates = [],
    positivity_constraint=0.,
    verbosity=0
    )
    variants_groups = Iterators.partition(variants, length(variants) ÷ Threads.nthreads())
    estimands_tasks = map(variants_groups) do variants
        Threads.@spawn estimands_from_variants(variants, dataset, AIE, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    end
    estimands_partitions = fetch.(estimands_tasks)
    return vcat(estimands_partitions...)
end

get_only_file_with_suffix(files, suffix) = files[only(findall(x -> endswith(x, suffix), files))]

function read_bed_chromosome(bedprefix)
    bed_files = files_matching_prefix(bedprefix)
    fam_file = get_only_file_with_suffix(bed_files, "fam")
    bim_file = get_only_file_with_suffix(bed_files, "bim")
    bed_file = get_only_file_with_suffix(bed_files, "bed")[1:end-4]
    return SnpData(bed_file, famnm=fam_file, bimnm=bim_file)
end

function get_genotypes_from_beds(bedprefix)
    snpdata = read_bed_chromosome(bedprefix)
    genotypes = DataFrame(convert(Matrix{UInt8}, snpdata.snparray), snpdata.snp_info."snpid")
    genotype_map = Union{UInt8, Missing}[0, missing, 1, 2]
    for col in names(genotypes)
        genotypes[!, col] = [genotype_map[x+1] for x in genotypes[!, col]]
    end
    insertcols!(genotypes, 1, :SAMPLE_ID => snpdata.person_info."iid")
    return genotypes
end

function make_genotypes(genotype_prefix, config, call_threshold)
    genotypes = if config["type"] == "gwas"|| config["type"] == "gwis"
        get_genotypes_from_beds(genotype_prefix)
    else
        variants_set = Set(retrieve_variants_list(config["variants"]))
        call_genotypes(genotype_prefix, variants_set, call_threshold)
    end
    return genotypes
end

function inputs_from_config(config_file, genotypes_prefix, traits_file, pcs_file;
    outprefix="final",
    batchsize=nothing,
    call_threshold=0.9,
    positivity_constraint=0.01,
    verbosity=0)

    traits = read_csv_file(traits_file)
    pcs = load_flash_pca_results(pcs_file)
    config = YAML.load_file(config_file)

    # Variables
    verbosity > 0 && @info("Parsing configuration file.")
    extra_treatments = haskey(config, "extra_treatments") ? Symbol.(config["extra_treatments"]) : []
    outcome_extra_covariates = haskey(config, "outcome_extra_covariates") ? Symbol.(config["outcome_extra_covariates"]) : []
    extra_confounders = haskey(config, "extra_confounders") ? Symbol.(config["extra_confounders"]) : []
    confounders = confounders_from_pcs(pcs, extra_confounders)
    nonoutcomes = Set(vcat(:SAMPLE_ID, extra_confounders, outcome_extra_covariates, extra_treatments))
    outcomes = filter(x -> x ∉ nonoutcomes, Symbol.(names(traits)))

    # Genotypes and final dataset
    verbosity > 0 && @info("Building and writing dataset.")
    genotypes = make_genotypes(genotypes_prefix, config, call_threshold)
    dataset = merge(traits, pcs, genotypes)
    Arrow.write(string(outprefix, ".data.arrow"), dataset)

    # Estimands
    config_type = config["type"]
    estimands = if config_type == "flat"
        estimands_from_flat_list(config["estimands"], dataset, config["variants"], outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity)
    elseif config_type == "groups"
        estimands_from_groups(config["estimands"], dataset, config["variants"], outcomes, confounders; 
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    elseif config_type == "gwas"
        variants = filter(!=("SAMPLE_ID"), names(genotypes))
        estimands_from_gwas(dataset, variants, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )
    elseif config_type == "gwis"
        variants = filter(!=("SAMPLE_ID"), names(genotypes))
        estimands_from_gwis(dataset, variants, outcomes, confounders;
            extra_treatments=extra_treatments,
            outcome_extra_covariates=outcome_extra_covariates,
            positivity_constraint=positivity_constraint,
            verbosity=verbosity
        )  
    else
        throw(ArgumentError(string("Unknown extraction type: ", config_type, ", use any of: (flat, groups, gwas)")))
    end

    @assert length(estimands) > 0 "No estimands left, probably due to a too high positivity constraint."

    write_estimation_inputs(outprefix, dataset, groups_ordering(estimands); batchsize=batchsize)

    verbosity > 0 && @info("Done.")

    return 0
end