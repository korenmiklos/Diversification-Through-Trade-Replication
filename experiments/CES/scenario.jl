fname = ARGS[1]

# load ImpvolEquilibrium first so that methods are accessible
include("../../equilibrium.jl")
using ImpvolEquilibrium, Logging
Logging.configure(level=INFO)

using FileIO, JLD2

sigma = parse(split(fname, "/")[1])
println("Sigma = $sigma")
parameters = load(fname)["parameters"]

@time results = map(t -> (t, ImpvolEquilibrium.period_wrapper(parameters, t)), 1:parameters[:T])
for t=1:parameters[:T], (k, v) in results[t][2]
    results[t][2][k] = results[t][2][k][:,:,:,1:1]
end
jldopen("$sigma/results.jld2", "w") do file
		file["results"] = results
end
