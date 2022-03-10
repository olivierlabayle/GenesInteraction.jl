###############################################################################
# BUILD TMLE FROM .TOML

function buildmodels(config)
    models = Dict()
    for (modelname, hyperparams) in config
        if !(modelname in ("resampling", "measures"))
            modeltype = eval(Symbol(modelname))
            paramnames = Tuple(Symbol(x[1]) for x in hyperparams)
            counter = 1
            for paramvals in Base.Iterators.product(values(hyperparams)...)
                model = modeltype(;NamedTuple{paramnames}(paramvals)...)
                models[Symbol(modelname*"_$counter")] = model
                counter += 1
            end
        end
    end
    return models
end


function stack_from_config(config::Dict, metalearner; adaptive_cv=true)
    # Define the resampling strategy
    resampling = config["resampling"]
    resampling = eval(Symbol(resampling["type"]))(nfolds=resampling["nfolds"])
    if adaptive_cv
        resampling = AdaptiveCV(resampling)
    end

    # Define the internal cross validation measures to report
    measures = config["measures"]
    measures = (measures === nothing || size(measures, 1) == 0) ? nothing : 
        [getfield(MLJBase, Symbol(fn)) for fn in measures]

    # Define the models library
    models = buildmodels(config)

    # Define the Stack
    Stack(;metalearner=metalearner, resampling=resampling, measures=measures, models...)
end


function estimators_from_toml(config::Dict, queries, outcome_type; adaptive_cv=true)
    # Parse estimator for the propensity score
    metalearner = LogisticClassifier(fit_intercept=false)
    if length(first(queries).case) > 1
        G = FullCategoricalJoint(stack_from_config(config["G"], metalearner, adaptive_cv=adaptive_cv))
    else
        G = stack_from_config(config["G"], metalearner, adaptive_cv=adaptive_cv)
    end
    
    # Parse estimator for the outcome regression
    if outcome_type <: AbstractFloat
        metalearner =  LinearRegressor(fit_intercept=false)
        Q̅ = stack_from_config(config["Qcont"], metalearner, adaptive_cv=adaptive_cv)
    elseif outcome_type <: Bool
        metalearner = LogisticClassifier(fit_intercept=false)
        Q̅ = stack_from_config(config["Qcat"], metalearner, adaptive_cv=adaptive_cv)
    else
        throw(ArgumentError("The type of the outcomes: $outcome_type, should be either a Float or a Bool"))
    end

    return TMLEstimator(Q̅, G, queries...)
end
