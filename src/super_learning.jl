using MLJ, MLJBase
using MLJLinearModels: LogisticClassifier
using DataFrames


###################################################
###   SuperLearner and Constructors             ###
###################################################

"""
    SuperLearner(library, metalearner, nfolds::Bool, shuffle::Int)

Implements the SuperLearner estimator described in :
van der Laan, Mark J.; Polley, Eric C.; and Hubbard, Alan E., "Super Learner" (July 2007).
U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 222.
https://biostats.bepress.com/ucbbiostat/paper222

Current assumptions:
   - Y is binary and a stratified cross validation is used
   - The estimated quantity is E[Y|X]

# Arguments

- `library`: An iterable over models to be trained
- `metalearner`: Model that will learn an optimal combination of the library members
- `nfolds::Int`: Number of folds for internal cross validation
- `shuffle::Bool`: Should the dataset be shuffled.

# Examples

TODO
"""
# mutable struct SuperLearner <: DeterministicComposite
#     library
#     metalearner
#     nfolds::Int
#     shuffle::Bool
# end


macro superlearner(exs...)

    nfolds_ex = :(nfolds::Int=10)
    shuffle_ex = :(shuffle::Bool=false)
    library_exs = []
    for ex in exs
        # Parse shuffle and nfolds values
        if ex.args[1] == :nfolds
            nfolds_ex.args[2] = ex.args[2]
        elseif ex.args[1] == :shuffle
            shuffle_ex.args[2] = ex.args[2]
        # Parse metalearner
        elseif ex.args[1] == :metalearner
            model = eval(ex.args[2])
            metatype = typeof(model)
            metalearner_ex = :(metalearner::$metatype=$model)
            
        # Parse library
        else
            println(ex)
        end
    end
    

    struct_ex = :(mutable struct SuperLearnerBis <: DeterministicComposite
        $(shuffle_ex)
        $(nfolds_ex)
        $(metalearner_ex)
    end)

    __module__.eval(MLJBase.Parameters.with_kw(struct_ex, __module__, false))

    esc(quote
        SuperLearnerBis()
        end)
end

metalearner = LogisticClassifier()
sl = @superlearner(LogisticClassifier(), metalearner=metalearner, nfolds=3, shuffle=false)

"""
    SuperLearner(library, metalearner;nfolds=10, shuffle=false)

"""
SuperLearner(library, metalearner;nfolds=10, shuffle=false) = SuperLearner(library, metalearner, nfolds, shuffle)


###################################################
###   Helper functions                          ###
###################################################


function is_uppercase(char::Char)
    i = Int(char)
    i > 64 && i < 91
end

"""
    snakecase(str, del='_')

Taken from StatisticalTraits.jl
Return the snake case version of the abstract string or symbol, `str`, as in

    snakecase("TheLASERBeam") == "the_laser_beam"

"""
function snakecase(str::AbstractString; delim='_')
    snake = Char[]
    n = length(str)
    for i in eachindex(str)
        char = str[i]
        if is_uppercase(char)
            if i != 1 && i < n &&
                !(is_uppercase(str[i + 1]) && is_uppercase(str[i - 1]))
                push!(snake, delim)
            end
            push!(snake, lowercase(char))
        else
            push!(snake, char)
        end
    end
    return join(snake)
end

snakecase(s::Symbol) = Symbol(snakecase(string(s)))


"""
Taken from MLJBase.jl
"""
function generate_name!(M, existing_names)
    str = split(string(M), '{') |> first
    candidate = split(str, '.') |> last |> snakecase |> Symbol
    candidate in existing_names ||
        (push!(existing_names, candidate); return candidate)
    n = 2
    new_candidate = candidate
    while true
        new_candidate = string(candidate, n) |> Symbol
        new_candidate in existing_names || break
        n += 1
    end
    push!(existing_names, new_candidate)
    return new_candidate
end


expected_value(X::AbstractNode) = node(XX->XX.prob_given_ref[2], X)

train_test_pairs(X::AbstractNode, cv, rows) = node(XX->MLJBase.train_test_pairs(cv, rows, XX), X)

train_restrict(X::AbstractNode, f::AbstractNode, i) = node((XX, ff) -> restrict(XX, ff[i], 1), X, f)

test_restrict(X::AbstractNode, f::AbstractNode, i) = node((XX, ff) -> restrict(XX, ff[i], 2), X, f)


###################################################
###   Fit function                              ###
###################################################


"""
    MLJ.fit(m::SuperLearner, verbosity::Int, X, y)
"""
function MLJ.fit(m::SuperLearner, verbosity::Int, X, y)
    stratified_cv = StratifiedCV(; nfolds=m.nfolds, shuffle=m.shuffle)
    n, = size(y)
    registry = Dict()

    X = source(X)
    y = source(y)

    Zval = []
    yval = []
    folds = train_test_pairs(y, stratified_cv, 1:n)
    for nfold in 1:m.nfolds
        Xtrain = train_restrict(X, folds, nfold)
        ytrain = train_restrict(y, folds, nfold)
        Xtest = test_restrict(X, folds, nfold)
        ytest = test_restrict(y, folds, nfold)

        Zfold = []
        for model in m.library
            mach = machine(model, Xtrain, ytrain)
            push!(Zfold, expected_value(predict(mach, Xtest)))
            registry[(name(model), nfold)] = mach
        end

        Zfold = hcat(Zfold...)
        
        push!(Zval, Zfold)
        push!(yval, ytest)
    end

    Zval = MLJ.table(vcat(Zval...))
    yval = vcat(yval...)

    metamach = machine(m.metalearner, Zval, yval)
    registry["metamachine"] = metamach

    Zpred = []
    for model in m.library
        mach = machine(model, X, y)
        push!(Zpred, expected_value(predict(mach, X)))
        registry[name(model)] = mach
    end

    Zpred = MLJ.table(hcat(Zpred...))
    ŷ = expected_value(predict(metamach, Zpred))

    mach = machine(Deterministic(), X, y; predict=ŷ)

    return!(mach, m, verbosity)

end
