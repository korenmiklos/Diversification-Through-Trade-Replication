# load ImpvolEquilibrium first so that methods are accessible
include("../../equilibrium.jl")
using ImpvolEquilibrium, FileIO, JLD2

## these are needed for data -> parameters mapping
parameters = Dict{Symbol, Any}()

# CES parameters
sigmas = [0.3 0.4 0.5 0.7 0.8 0.9 0.999 1.1 1.2 1.3 1.5 1.75 2.0]
parameters[:theta] = 4.0
parameters[:eta] = 4.0

########## parameters common across scenarios
## these are function of data
# inverse of adjustment cost, 0 if cannot readjust
parameters[:one_over_rho] = 0.0

for sigma in sigmas
	parameters[:sigma] = sigma
	include("../config.jl")

	try
		mkdir("$sigma")
	catch
	end
	jldopen("$sigma/common_parameters.jld2", "w") do file
		file["parameters"] = parameters
	end
end
