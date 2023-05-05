actors_error() = throw(ArgumentError("At least two of (bqtls, env, trans-actors) should be specified"))

treatments_from_actors(::Nothing, ::Nothing, _) = actors_error()
treatments_from_actors(_, ::Nothing, ::Nothing) = actors_error()
treatments_from_actors(::Nothing, _, ::Nothing) = actors_error()

function treatments_from_actors(bqtl_file, env_file, trans_actors_prefix)
    bqtls = read_snps_from_csv(bqtl_file)
    extraT  = read_txt_file(env_file)
    transactors = trans_actors_from_prefix(trans_actors_prefix)
    bqtls, transactors, extraT
end


combine_trans_actors(trans_actors::Vector{DataFrame}, extraT::DataFrame, order) = combinations([trans_actors..., extraT], order)
combine_trans_actors(trans_actors::Vector{DataFrame}, extraT::Nothing, order) = combinations(trans_actors, order)
combine_trans_actors(trans_actors::Nothing, extraT::DataFrame, order) = [[extraT]]


function combine_by_bqtl(bqtls::DataFrame, trans_actors::Union{Vector{DataFrame}, Nothing}, extraT::Union{DataFrame, Nothing}, order::Int)
    treatment_combinations = Vector{Symbol}[]
    if order == 1
        for bqtl in bqtls.ID
            push!(treatment_combinations, [Symbol(bqtl)])
        end
    else
        for comb in combine_trans_actors(trans_actors, extraT, order - 1)
            comb_interactions = crossjoin(bqtls, comb..., makeunique=true)
            for row in eachrow(comb_interactions)
                treatment_combination = Symbol.(vcat(row["ID"], [row[string("ID_", i)] for i in 1:order-1]))
                push!(treatment_combinations, treatment_combination)
            end
        end
    end
    return treatment_combinations
end

parse_orders(orders::String) = parse.(Int, split(orders, ","))


all_variants(bqtls::DataFrame, transactors::Nothing) = Set(bqtls.ID)
all_variants(bqtls::Nothing, transactors::Vector{DataFrame}) = Set(vcat((x.ID for x in transactors)...))
all_variants(bqtls::DataFrame, transactors::Vector{DataFrame}) = Set(vcat(bqtls.ID, (x.ID for x in transactors)...))


read_snps_from_csv(path::Nothing) = nothing
read_snps_from_csv(path::String) = unique(CSV.read(path, DataFrame; select=[:ID, :CHR]), :ID)


trans_actors_from_prefix(trans_actors_prefix::Nothing) = nothing
function trans_actors_from_prefix(trans_actors_prefix::String)
    dir, prefix = splitdir(trans_actors_prefix)
    dir_ = dir == "" ? "." : dir
    trans_actors = DataFrame[]
    for filename in readdir(dir_)
        if startswith(filename, prefix)
            push!(trans_actors, read_snps_from_csv(joinpath(dir, filename)))
        end
    end
    return trans_actors
end

sorted_unique_values(v) = sort(unique(skipmissing(v)))

"""
    control_case_settings(::Type{IATE}, treatment_tuple, data::DataFrame)

Generates all possible IATE parameter treatment settings where control ≠ case. 
The procedure works as follows:

    - For each treatment, control ≠ case are generated from pairwise combinations of unique values
    - The cartesian product of all those treatment case/control values is returned
"""
function control_case_settings(::Type{IATE}, treatments, data::DataFrame)
    control_cases_ = []
    # For each treatment variable generate the case/control options
    # where case and control are different
    for t in treatments
        unique_vals = sorted_unique_values(data[!, t])
        push!(control_cases_, [(case=case, control=control) for (case,control) ∈ combinations(unique_vals, 2)])
    end
    # Generate the product of all treatment case/control settings
    return (
        NamedTuple{tuple(Symbol.(treatments)...)}(x) 
        for x in Iterators.product(control_cases_...)
    )
end

"""
    control_case_settings(::Type{ATE}, treatment_tuple, data)

Generates all possible ATE parameter treatment settings where only the bQTL has control ≠ case, 
for the rest of the treatment variables: case=control. The bQTl is assumed to 
be the first variable in the `treatment_tuple`.

The procedure works as follows:

    - For the bQTL, control ≠ case are generated from pairwise combinations of unique alleles
    - For all other treatment variables, control = case are the unique values
    - The cartesian product of all those treatment case/control values is returned
"""
function control_case_settings(::Type{ATE}, treatments, data)
    control_cases_ = []
    # The bQTL is the first treatment variable
    # We get all case/control values
    bQTL_values = TargeneCore.sorted_unique_values(data[!, treatments[1]])
    push!(control_cases_, [(case=case, control=control) for (case,control) in combinations(bQTL_values, 2)])
    # The remainder of the treatment are kept fixed
    for t in treatments[2:end]
        unique_vals = TargeneCore.sorted_unique_values(data[!, t])
        push!(control_cases_, [(case=uv, control=uv) for uv in unique_vals])
    end
    # Generate the product of all treatment case/control settings
    return (
        NamedTuple{tuple(Symbol.(treatments)...)}(x) 
        for x in Iterators.product(control_cases_...)
    )
end

function addParameters!(parameters, treatments, variables, data; positivity_constraint=0.)
    freqs = TargeneCore.frequency_table(data, treatments)
    # This loop adds all ATE parameters where all other treatments than
    # the bQTL are fixed, at the order 1, this is the simple bQTL's ATE
    for setting in control_case_settings(ATE, treatments, data)
        Ψ = ATE(target=Symbol("*"), treatment=setting, confounders=variables.confounders, covariates=variables.covariates)
        if satisfies_positivity(Ψ, freqs; positivity_constraint=positivity_constraint)
            update_parameters_from_targets!(parameters, Ψ, variables.targets)
        end
    end
    # This loop adds all IATE parameters that pass the positivity threshold
    if size(treatments, 1) >= 2
        for setting in control_case_settings(IATE, treatments, data)
            Ψ = IATE(target=Symbol("*"), treatment=setting, confounders=variables.confounders, covariates=variables.covariates)
            if satisfies_positivity(Ψ, freqs; positivity_constraint=positivity_constraint)
                update_parameters_from_targets!(parameters, Ψ, variables.targets)
            end
        end
    end
end

function get_variables(pcs, traits, extraW, extraC, extraT)
    nontargets = Set(vcat("SAMPLE_ID", extraW, extraC, extraT))
    return (
        extra_treatments = extraT isa Nothing ? nothing : Symbol.(extraT),
        covariates = extraC isa Nothing ? [] : Symbol.(extraC),
        confounders = all_confounders(pcs, extraW),
        targets = targets_from_traits(traits, nontargets)
    )
end

function parameters_from_actors(bqtls, transactors, data, variables, orders; positivity_constraint=0.)
    parameters = TMLE.Parameter[]
    extraT_df = variables.extra_treatments isa Nothing ? nothing : DataFrame(ID=variables.extra_treatments)
    # For each interaction order generate parameter files
    for order in orders
        # First generate the `T` section
        treatment_combinations = TargeneCore.combine_by_bqtl(bqtls, transactors, extraT_df, order)
        for treatments in treatment_combinations
            # Generate `Parameters` section
            addParameters!(parameters, treatments, variables, data; positivity_constraint=positivity_constraint)
        end
    end

    length(parameters) > 0 || throw(NoRemainingParamsError(positivity_constraint))
    optimize_ordering!(parameters)

    return parameters
end


function tmle_inputs_from_actors(parsed_args)
    batch_size = parsed_args["batch-size"]
    outprefix = parsed_args["out-prefix"]
    call_threshold = parsed_args["call-threshold"]
    bgen_prefix = parsed_args["bgen-prefix"]
    positivity_constraint = parsed_args["positivity-constraint"]
    traits = TargeneCore.read_data(parsed_args["traits"])
    pcs = TargeneCore.read_data(parsed_args["pcs"])
    orders = TargeneCore.parse_orders(parsed_args["from-actors"]["orders"])
    extraW = TargeneCore.read_txt_file(parsed_args["from-actors"]["extra-confounders"])
    extraC = TargeneCore.read_txt_file(parsed_args["from-actors"]["extra-covariates"])
    genotypes_asint = parsed_args["from-actors"]["genotypes-as-int"]

    # Retrieve SNPs and environmental treatments
    bqtls, transactors, extraT = TargeneCore.treatments_from_actors(
        parsed_args["from-actors"]["bqtls"], 
        parsed_args["from-actors"]["extra-treatments"], 
        parsed_args["from-actors"]["trans-actors-prefix"]
    )
    # Genotypes and final dataset
    variants = TargeneCore.all_variants(bqtls, transactors)
    genotypes = TargeneCore.call_genotypes(bgen_prefix, variants, call_threshold; asint=genotypes_asint)
    data = TargeneCore.merge(traits, pcs, genotypes)

    # Parameter files
    variables = TargeneCore.get_variables(pcs, traits, extraW, extraC, extraT)
    parameters = TargeneCore.parameters_from_actors(
        bqtls, transactors, data, variables, orders; 
        positivity_constraint=positivity_constraint
    )

    # write data and parameter files
    write_tmle_inputs(outprefix, data, parameters; batch_size=batch_size)
end