module ImpvolOutput
	using JLD2
	using FileIO
	using CSV
	using Missings
	using Plots
	using HypothesisTests
	using Distributions
	using DataFrames
	using Logging
	using FilePathsBase

	include("calibration_utils.jl")
	global_logger(ConsoleLogger(stderr, Logging.Debug))

	function read_results(path = "experiments/baseline/actual/results.jld2")
		file = jldopen(path, "r")
		return file["results"]
	end

	function sort_results(results)
		return sort(collect(results), by = x -> x[1])
	end

	function list_keys(results)
		for (key, value) in results[1][2]
			println(key)
		end
	end

	function make_series(results, key = :real_GDP)
		series = zeros(size(results[1][2][key])[1], size(results[1][2][key])[2], size(results[1][2][key])[3], length(results))
		for t in 1:length(results)
			# Only the first element in the shock dimension is interesting, the rest are there only for optimizaton purposes used in the algorithm
			series[:,:,:,t] = results[t][2][key][:,:,:,1]
		end
		return series
	end

	function calculate_volatilities(x, parameters, bool_detrend::Bool, range=:)
		if bool_detrend
			x_c, x_t = DetrendUtilities.detrend(log.(x),parameters[:bp_weights])
		else
			x_c = log.(x)
		end

		# time is the last dimension
		return var(x_c[:,:,:,range],dims=ndims(x_c))
	end

	function plot_model_vs_data(plot_data::Tuple, title::String)
		# data and model are expected to be of dimension N x T
		# length of label should be equal N
		data = plot_data[1]
		model = plot_data[2]
		label = plot_data[3]

		data = permutedims(data,(2,1))
		model = permutedims(model,(2,1))

		size(data,2) == size(model,2) || error("You can only compare matching time series between the model outcome and data")

		colors = distinguishable_colors(size(data,2))

		ENV["GKSwstype"] = "gksqt"
		fig = Plots.plot()
		for i in 1:size(data,2)
			Plots.plot!([data[:,i] model[:,i]], color = colors[i,1], ls = [:solid :dash], label = [label[i,1] ""], title = title)
		end
		return fig
	end


	function plot_data(path = "experiments/baseline/actual/results.jld2", key = :real_GDP, country_range = ":")
		# output: "data, model, label", as an input for the function 'plot_model_vs_data'
		data = load("data/impvol_data.jld2")
		gdp_d = squeeze(sum(data["va"],3), (1,3))
		gdp_d = Float64.(collect(Missings.replace(gdp_d, NaN)))
		gdp_d = gdp_d[eval(parse(country_range)),:]

		pwt = squeeze(data["pwt"], (1,3))
		pwt = Float64.(collect(Missings.replace(pwt, NaN)))
		pwt = pwt[eval(parse(country_range)),:]

		cpi = CSV.read("data/cpi.csv", header = true)
		cpi = permutedims(convert(Array, cpi)[:,2:end], (2,1))
		cpi = cpi ./ cat(2, cpi[:,24]) # 24 = 1995-base
		cpi = Float64.(collect(Missings.replace(cpi, NaN)))
		cpi = cpi[eval(parse(country_range)),:]

		xr = CSV.read("data/exchange_rates.csv", header = true)
		xr = permutedims(convert(Array, xr)[:,2:end], (2,1))
		xr = Float64.(collect(Missings.replace(xr, NaN)))
		xr = xr[eval(parse(country_range)),:]

		country_names = CSV.read("data/country_name.txt", header = false, types = [String])
		country_names = country_names[:1][eval(parse(country_range))]

		results = sort_results(read_results(path))
		gdp_m = make_series(results, key)
		gdp_m = sum(gdp_m[1,eval(parse(country_range)),:,:],2)
		gdp_m = squeeze(gdp_m, 2)

		gdp_d = gdp_d ./ cpi .* xr
		#gdp_d = gdp_d ./ pwt

		# normalization: model(1972) = data(1972)
		gdp_m = gdp_m .* (gdp_d[:,1] ./ gdp_m[:,1])

		return log.(gdp_d), log.(gdp_m), country_names
	end

	################ Running the 'plot_data' function ################
	# include("output.jl")
	# plt_dta = ImpvolOutput.plot_data() # with the desired arguments if needed
	# Base.print_matrix(IOContext(STDOUT, :limit => false), plt_dta[3]) # to check the countries' position
	# idx = zeros(length(plt_dta[3]))
	# idx[6], idx[10], idx[11], idx[13], idx[15], idx[16], idx[24], idx[25] = 1, 1, 1, 1, 1, 1, 1, 1
	# idx = idx .== 1
	# ImpvolOutput.plot_model_vs_data((plt_dta[1][idx,:], plt_dta[2][idx,:], plt_dta[3][idx]), "GDP")
	# # solid line = data, dashed line = model
	##################################################################


	function write_results(parameters, rootpath = "experiments/baseline/", key = :real_GDP, bool_detrend = true, dirsdown = 1, pattern = r"jld2$")
		# Volatility of the desired variable found in the last folders of depth 'dirsdown' from the 'rootpath' is calculated by this function
		# Read country names
		country_names = CSV.read("data/country_name.txt", DataFrame; header = false, types = [String])
		println("Dimensions of country_names: ", size(country_names))
	
		# Create 'stats' DataFrame to store volatilities with the correct number of rows
		
		stats = DataFrame(
			country_names = country_names[!, 1],
			actual = Vector{Float64}(undef, size(country_names, 1)),
			kappa1972 = Vector{Float64}(undef, size(country_names, 1)),
			nosectoral = Vector{Float64}(undef, size(country_names, 1)),
			nosectoral_kappa1972 = Vector{Float64}(undef, size(country_names, 1)),
			onlyglobal = Vector{Float64}(undef, size(country_names, 1)),
			onlyglobal_kappa1972 = Vector{Float64}(undef, size(country_names, 1)),
			trade_barriers = Vector{Float64}(undef, size(country_names, 1)),
			diversification = Vector{Float64}(undef, size(country_names, 1)),
			specialization = Vector{Float64}(undef, size(country_names, 1))
		)
		println("Dimensions of stats: ", size(stats))
	
	
		# Fill 'stats' DataFrame
		for (root, dirs, files) in walkdir(dirname(rootpath))
			dir_depth = length(collect(eachmatch(r"/", root))) - 1 #count(x -> x == '/', root) - 1 
	
			if dir_depth == dirsdown
				for file in files
					if occursin(pattern, file)
						sectoral_series = make_series(sort_results(read_results(joinpath(root, file))), key)

						# Summation over the sectors
						series = cat(dims=ndims(sectoral_series), sum(sectoral_series, dims=3))
						symbol_name = Symbol(SubString(root, findall(x -> x == '/', root)[end] + 1, length(root)))
						stats[!,symbol_name] = dropdims(calculate_volatilities(series, parameters, bool_detrend), dims=(1,3,4))
					end
				end
			end
		end
	
		@debug stats
	
		# Trade barriers
		stats.trade_barriers = 100 * (stats.actual .- stats.kappa1972) ./ stats.kappa1972
	
		# Diversification
		stats.diversification = 100 * (stats.nosectoral .- stats.nosectoral_kappa1972) ./ stats.kappa1972
	
		# Specialization
		stats.specialization = 100 * (stats.actual .- stats.kappa1972 .- stats.nosectoral .+ stats.nosectoral_kappa1972) ./ stats.kappa1972
	
		CSV.write(rootpath * "/output_table.csv", stats)
	end
	

#= 
	function write_results(parameters, rootpath = "experiments/baseline/", key = :real_GDP, bool_detrend = true, dirsdown = 1, pattern = r"jld2$")
		# Volatility of the desired variable found in the last folders of depth 'dirsdown' from the 'rootpath' is calculated by this function

		# Create 'stats' array to store volatilities
		stats = DataFrame([String, Real, Real, Real, Real, Real, Real, Real, Real, Real], [:country_names, :actual, :kappa1972, :nosectoral, :nosectoral_kappa1972, :onlyglobal, :onlyglobal_kappa1972, :trade_barriers, :diversification, :specialization], 25)
		println("Dimensions of stats: ", size(stats))
		stats[:country_names] = CSV.read("data/country_name.txt", header = false, types = [String], nullable=false)[1]

		# Fill 'stats' array
		for (root, dirs, files) in walkdir(dirname(rootpath))
			dir_depth = length(matchall(r"/", root)) - 1

			if dir_depth == dirsdown

				for file in files
					if ismatch(pattern, file)
						sectoral_series = make_series(sort_results(read_results(joinpath(root, file))), key)

						# Summation over the sectors
						series = cat(ndims(sectoral_series), sum(sectoral_series, 3))
						stats[eval(parse(":" * SubString(root,findin(root,'/')[end] + 1, length(root))))] = squeeze(calculate_volatilities(series, parameters, bool_detrend), (1,3,4))
					end
				end
			end
		end

		Logging.debug(stats)

		# Trade barriers
		stats[:trade_barriers] = 100 * (stats[:actual] - stats[:kappa1972]) ./ stats[:kappa1972]

		# Diversification
		stats[:diversification] = 100 * (stats[:nosectoral] - stats[:nosectoral_kappa1972]) ./ stats[:kappa1972]

		# Specialization
		stats[:specialization] = 100 * (stats[:actual] - stats[:kappa1972] - stats[:nosectoral] + stats[:nosectoral_kappa1972]) ./ stats[:kappa1972]

		CSV.write(rootpath * "/output_table.csv", stats)
	end
 =#
	################ Running the 'write_results' function ################
	# include("output.jl")
	# parameters = Dict{Symbol, Any}()
	# parameters[:N] = 25
	# parameters[:bp_weights] = [0.774074394803123, -0.201004684236153, -0.135080548288772,-0.050951964876636]
	# stats = ImpvolOutput.write_results(parameters)
	######################################################################


end

