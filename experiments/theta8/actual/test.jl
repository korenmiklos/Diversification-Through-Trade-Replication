# load ImpvolEquilibrium first so that methods are accessible
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium, Logging
@everywhere Logging.configure(level=INFO)

@everywhere using FileIO
@everywhere parameters = load("../common_parameters.jld2")["parameters"]

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")
@everywhere parameters[:S] = 2

@everywhere srand(7094)

T = 3
@time results = pmap(t -> (t, ImpvolEquilibrium.period_wrapper(parameters, t)), 1:T)
GDP = zeros(T,)
for t=1:T
	GDP[t] = sum(results[t][2][:real_GDP][1,end,:,1])
end
display(GDP)