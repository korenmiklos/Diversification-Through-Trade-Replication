# load ImpvolEquilibrium first so that methods are accessible
include("../../equilibrium.jl")
using ImpvolEquilibrium, FileIO, JLD2

## these are needed for data -> parameters mapping
parameters = Dict{Symbol, Any}()

# CES parameters
parameters[:sigma] = 0.999
parameters[:theta] = 4.0
parameters[:eta] = 4.0

########## parameters common across scenarios
## these are function of data
# inverse of adjustment cost, 0 if cannot readjust
parameters[:one_over_rho] = 0.001

include("../config.jl")
# change parameters after reading data, but common across scenarios
# change parameters after reading data, but common across scenarios
jldopen("common_parameters.jld2", "w") do file
	file["parameters"] = parameters
end
