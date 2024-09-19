# load ImpvolEquilibrium first so that methods are accessible
include("../../equilibrium.jl")
using .ImpvolEquilibrium, FileIO, JLD2

## these are needed for data -> parameters mapping
parameters = Dict{Symbol, Any}()

# CES parameters
parameters[:sigma] = 0.999
parameters[:theta] = 4.0
parameters[:eta] = 4.0

########## parameters common across scenarios
## these are function of data
# inverse of adjustment cost, 0 if cannot readjust
parameters[:one_over_rho] = 0.0

include("../config.jl")
using FileIO
data = load("../../data/impvol_data.jld2")
# change parameters after reading data, but common across scenarios
## Balanced trade
parameters[:beta_j] = ones(size(parameters[:beta_j]))
parameters[:gamma_jk] = zeros(size(parameters[:gamma_jk]))
parameters[:B] = CalibrateParameters.calculate_B(parameters)
parameters[:A] = CalibrateParameters.calculate_A(parameters, data)
CalibrateParameters.decompose_shocks!(parameters, parameters[:importance_weight])
CalibrateParameters.draw_productivity_shocks!(parameters)

parameters[:nu_njt] = CalibrateParameters.compute_alpha(parameters, data)
# change parameters after reading data, but common across scenarios
jldopen("common_parameters.jld2", "w") do file
	file["parameters"] = parameters
end
