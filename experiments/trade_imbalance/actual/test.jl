# load ImpvolEquilibrium first so that methods are accessible
@everywhere include("../../../equilibrium.jl")
@everywhere using ImpvolEquilibrium

@everywhere include("../init_parameters.jl")

# parameters that govern counterfactual
@everywhere include("change_parameters.jl")

display([parameters[:S_nt][1,:,1,23] parameters[:w_njt][1,:,1,23]])

@time results = ImpvolEquilibrium.period_wrapper(parameters, 23)