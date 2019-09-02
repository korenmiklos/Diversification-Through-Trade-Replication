# load ImpvolEquilibrium first so that methods are accessible
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium, Logging
@everywhere Logging.configure(level=INFO)

@everywhere using FileIO, JLD2
@everywhere parameters = load("../common_parameters.jld2")["parameters"]

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")

@time results = pmap(t -> (t, ImpvolEquilibrium.period_wrapper(parameters, t)), 1:parameters[:T])
include("../../../save.jl")