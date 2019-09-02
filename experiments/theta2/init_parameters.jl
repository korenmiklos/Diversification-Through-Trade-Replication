# load ImpvolEquilibrium first so that methods are accessible
include("../../equilibrium.jl")
using ImpvolEquilibrium, FileIO, JLD2

## these are needed for data -> parameters mapping
parameters = Dict{Symbol, Any}()

# CES parameters
parameters[:sigma] = 0.999
parameters[:theta] = 2.0
# we also have to change eta so that Gamma function is finite. eta only affects overall level of productivity
parameters[:eta] = 2.0

########## parameters common across scenarios
## these are function of data
# inverse of adjustment cost, 0 if cannot readjust
parameters[:one_over_rho] = 0.0

include("../config.jl")
# change parameters after reading data, but common across scenarios
jldopen("common_parameters.jld2", "w") do file
	file["parameters"] = parameters
end
