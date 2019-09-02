for t=1:parameters[:T], (k, v) in results[t][2]
    results[t][2][k] = results[t][2][k][:,:,:,1:1]
end
jldopen("results.jld2", "w") do file
		file["results"] = results
end
