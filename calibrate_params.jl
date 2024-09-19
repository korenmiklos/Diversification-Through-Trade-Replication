module CalibrateParameters
	include("calibration_utils.jl")
	using JLD2, FileIO, ..ImpvolEquilibrium, Statistics, LinearAlgebra, Logging
	import SpecialFunctions
	const gamma_function = SpecialFunctions.gamma

	function eye(n)
		return Matrix{Float64}(I, n, n)
	end

	function calibrate_parameters!(parameters, fname="../../data/impvol_data.jld2")
		data = load(fname)

		_, N, J, T = size(data["beta"])
		parameters[:N], parameters[:J], parameters[:T] = N, J, T
		parameters[:beta_j] = mean(data["beta"]; dims=(1,2,4))

		parameters[:gamma_jk] = compute_gamma(parameters, data)

		parameters[:S_nt] = zeros(1,N,1,T)

		parameters[:d] = expenditure_shares(parameters, data)

		parameters[:kappa_mnjt] = trade_costs(parameters)

		final_expenditure_shares = calculate_expenditure_shares(parameters, data)
		parameters[:final_expenditure_shares] = final_expenditure_shares

		# broad country weights for final expenditure
		country_weights = sum(data["va"], dims=(1,3,4))
		country_weights = country_weights ./ sum(country_weights, dims=2)

		calculate_p_and_nu!(parameters, data, final_expenditure_shares, country_weights)

		parameters[:w_njt] = calculate_nominal_wages(parameters, data)

		parameters[:B_j] = calculate_B(parameters)

		parameters[:xi] = calculate_xi(parameters)

		parameters[:A] = calculate_A(parameters, data)

		# total world expenditure in the data .- needed to get reasonable starting values
		parameters[:nominal_world_expenditure] = sum(data["va"] ./ parameters[:beta_j], dims=(1,2,3))
		# deflate trade imbalance to 1972 dollars
		deflator = CES_price_index(parameters[:nu_njt][:,end:end,:,:], parameters[:p_sectoral][:,end:end,:,:], parameters[:sigma])
		@info deflator[:]

		parameters[:S_nt_data] = (data["trade_balance"] .- mean(data["trade_balance"]; dims=2)) ./ deflator

		# global, all-time average of sector final expenditure shares
		importance_weight = mean(parameters[:nu_njt], dims=(1, 2, 4))
		# special-case CES, when nu does not have direct meaning
		if abs(parameters[:sigma]-1)>0.1
			share = data["va"] ./ sum(data["va"], dims=3)
			importance_weight = mean(share, dims=(1, 2, 4))
		end
		parameters[:importance_weight] = importance_weight
		decompose_shocks!(parameters, importance_weight)
		draw_productivity_shocks!(parameters)
	end

	function compute_gamma(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		beta = parameters[:beta_j]
		io_values = data["io_values"]
		total_output = data["total_output"]
		output_shares = data["output_shares"]
		intermediate_input_shares = data["intermediate_input_shares"]

		# Summing sectors
		# Agriculture & mining
		io_values[1,:,:,:]  = sum(io_values[1:2,:,:,:], dims=1)
		io_values         = io_values[setdiff(1:end,2),:,:,:]
		io_values[:,1,:,:]  = sum(io_values[:,1:2,:,:], dims=2)
		io_values         = io_values[:,setdiff(1:end,2),:,:]
		total_output[:,1,:,:] = sum(total_output[:,1:2,:,:], dims=2)
		total_output      = total_output[:,setdiff(1:end,2),:,:]

		# Services
		io_values[18,:,:,:]  = sum(io_values[18:end,:,:,:], dims=1)
		io_values          = io_values[setdiff(1:end,19:end),:,:,:]
		io_values[:,18,:,:]  = sum(io_values[:,18:end,:,:], dims=2)
		io_values          = io_values[:,setdiff(1:end,19:end),:,:]
		total_output[:,18,:,:] = sum(total_output[:,18:end,:,:], dims=2)
		total_output       = total_output[:,setdiff(1:end,19:end),:,:]

		# Split rows
		dupl_idx = [1, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 13, 14, 15, 16, 17, 18]
		io_values_dupl = io_values[dupl_idx,:,:,:]

		output_shares_full = ones(size(io_values_dupl))
		split_idx = [2, 3, 4, 5, 6, 8, 9, 17, 18, 19]
		output_shares_full[split_idx,:,:,:] = permutedims(repeat(output_shares, outer = [1,1,size(io_values_dupl,2),1]),(2,3,1,4))

		io_values_new = io_values_dupl .* output_shares_full

		# Split columns
		io_values_dupl = io_values_new[:,dupl_idx,:,:]

		intermediate_input_shares_full = ones(size(io_values_dupl))
		intermediate_input_shares_full[:,split_idx,:,:] = repeat(intermediate_input_shares, outer = [size(io_values_dupl,1),1,1,1])

		io_values_new = io_values_dupl .* intermediate_input_shares_full

		total_output = total_output[:,dupl_idx,:,:]
		output_shares_full = ones(size(total_output))
		output_shares_full[:,split_idx,:,:] = output_shares
		total_output = total_output .* output_shares_full

		# Correct order of sectors
		io_values_new[18,:,:,:], io_values_new[19,:,:,:], io_values_new[20,:,:,:] = io_values_new[20,:,:,:], io_values_new[18,:,:,:], io_values_new[19,:,:,:]
		io_values_new[:,18,:,:], io_values_new[:,19,:,:], io_values_new[:,20,:,:] = io_values_new[:,20,:,:], io_values_new[:,18,:,:], io_values_new[:,19,:,:]
		total_output[:,18,:,:], total_output[:,19,:,:], total_output[:,20,:,:] = total_output[:,20,:,:], total_output[:,18,:,:], total_output[:,19,:,:]

		# Compute gamma
		gamma = io_values_new ./ repeat(total_output, outer = [size(io_values_new,1),1,1,1])
		gamma = mean(gamma, dims=4)
		gamma = gamma .* permutedims(1 .- beta,(1,3,2,4)) ./ sum(gamma, dims=1)
		return gamma = dropdims(gamma, dims=(3,4))
	end

	function compute_alpha(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		va = data["va"]
		beta = parameters[:beta_j]
		gamma = parameters[:gamma_jk]
		weights = parameters[:bp_weights]

		alpha = zeros(J,T)

		for t in 1:T
			va_t = transpose(sum(va[1,:,:,t], dims=1))
			alpha[:,t] = (eye(J) .- gamma) * diagm(1 ./ beta[:]) * va_t / sum(va_t)
		end

		# Replace negative elements with 0
		alpha = (alpha .+ abs.(alpha)) / 2

		# Smooth the series
		alpha_c, alpha_t = DetrendUtilities.detrend(alpha, weights)

		# Normalization
		alpha = alpha_t ./ sum(alpha_t, dims=1)

		return alpha = permutedims(cat(alpha, dims=ndims(alpha) + 2), (3,4,1,2))
	end

	function trade_costs(parameters)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		d = parameters[:d]
		theta = parameters[:theta]
		n_zero = parameters[:numerical_zero]

		# smooth kappa over time so that it does not introduce shocks
		weights = parameters[:bp_weights]
		_, share_trend = DetrendUtilities.detrend(d, weights)
		# DEBUG: undo all this
		share_trend = d

		kappa = zeros(size(d))
		for j in 1:(J-1)
			for t = 1:T
				kappa[:,:,j,t] = ((share_trend[:,:,j,t] .* transpose(share_trend[:,:,j,t])) ./ (diag(share_trend[:,:,j,t]) * transpose(diag(share_trend[:,:,j,t])))).^(1 / (2 * theta))
			end
		end

		kappa[kappa .< n_zero] .= n_zero
		kappa[:,:,end:end,:] .= repeat(eye(N), outer = [1,1,1,T]) # Services

		kappa = min.(kappa,1)
		# smooth kappa over time so that it does not introduce shocks
		_, kappa_trend = DetrendUtilities.detrend(kappa, weights)
		# DEBUG: return kappa, not kappa_trend
		return kappa
	end

	function expenditure_shares(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		import_shares = data["import_shares"]
		n_zero = parameters[:numerical_zero]

		d = import_shares

		within_import = d ./ sum(d, dims=2)
		domestic_per_import = 1 ./ sum(d, dims=2) .- 1
		domestic_per_import[domestic_per_import .< n_zero] .= n_zero
		d = within_import ./ (1 .+ domestic_per_import)

		for n in 1:N
			d[n,n,:,:] = ones(J,T) .- dropdims(sum(d[n,:,:,:], dims=1), dims=1)
		end

		d[d .< n_zero] .= n_zero

		return d ./ sum(d, dims=2)
	end

	function calculate_xi(parameters)
		theta = parameters[:theta]
		eta = parameters[:eta]

		return gamma_function((theta + 1 - eta)/theta)
	end

	function calculate_B(parameters)
		beta = parameters[:beta_j]
		gamma = parameters[:gamma_jk]

		gamma = permutedims(cat(gamma, dims=ndims(gamma) + 2), [1,3,2,4])
		return B = (beta .^ -beta) .* prod(gamma .^ -gamma, dims=1)
	end

	function calculate_nominal_wages(parameters, data)
		nulla = parameters[:numerical_zero]
		weights = parameters[:bp_weights]
		value_added_shares = data["va"] ./ sum(data["va"], dims=3)
		V_c, V_t = DetrendUtilities.detrend(value_added_shares, weights)

		if parameters[:one_over_rho]>0.0
			trend = 0.5*(V_t .- parameters[:one_over_rho])
			labor_share = trend .+ (trend .^2 .+ parameters[:one_over_rho]*value_added_shares) .^ 0.5
			wage_ratio = value_added_shares ./ labor_share
			@info("Unweighted wage ratio should be 1: ", mean(wage_ratio))
		else
			# if no labor adjustment, the ratio of value added = the ratio of wages
			wage_ratio = value_added_shares ./ V_t
		end
		nominal_GDP = sum(data["va"], dims=3)
		return nominal_GDP .* wage_ratio
	end


	function calculate_p_and_nu!(parameters, data, final_expenditure_shares, country_weights)
		N = parameters[:N]
		J = parameters[:J]
		T = parameters[:T]

		p_sectoral_data = data["p_sectoral_data"]
		d = parameters[:d]
		kappa = parameters[:kappa_mnjt]
		theta = parameters[:theta]
		sigma = parameters[:sigma]
		nulla = parameters[:numerical_zero]

		p_sectoral_base = p_sectoral_data ./ p_sectoral_data[:,:,:,1]

		# step 1: calculate nu and price index for base country
		# NB: US as the last country in the matrix
		p_sectoral_US = p_sectoral_base[:,end:end,:,:]
		nu_US = final_expenditure_shares[:,end:end,:,:] .* p_sectoral_US .^ (sigma-1)
		nu_US = nu_US ./ sum(nu_US, dims=3)
		P_US = CES_price_index(nu_US, p_sectoral_US, sigma)
		@assert p_sectoral_US[1,1,:,1] ≈ ones(J) atol=1e-9
		@assert P_US[1,1,1,1] ≈ 1.0 atol=1e-9

		# step 2: calculate sectoral prices from market shares relative to US
		# US is assumed to be chosen as a base country (US = end), else pwt should be used to do the conversion
		# normalization: p_sectoral[1,end,:,1] = 1.0
		# Permute the dimensions of the last slice of d and kappa
		d_permuted = permutedims(cat(d[end, :, :, :], dims=ndims(d)), [4, 1, 2, 3])
		kappa_permuted = permutedims(cat(kappa[end, :, :, :], dims=ndims(kappa)), [4, 1, 2, 3])

		# Compute the mean after performing element-wise operations
		log_term = log.(d ./ d_permuted) .- log.(kappa ./ kappa_permuted)
		mean_log_term = mean(1 / theta * log_term, dims=2)

		# Permute the dimensions of the log-transformed p_sectoral_base and repeat it
		p_sectoral_base_permuted = permutedims(cat(log.(p_sectoral_base[:, end, :, :]), dims=ndims(p_sectoral_base)), [1, 4, 2, 3])

		# Compute the exponential of the mean log term and add the repeated p_sectoral_base
		# no need to repeat p_sectoral_base_permuted, this is not MATLAB
		result = exp.(mean_log_term .+ p_sectoral_base_permuted)

		# Transpose the resulting array
		p_sectoral = array_transpose(result)

		@assert any(isnan, p_sectoral[:,:,1:end-1,:]) == false
		# step 3: calculate tradable nu and infer nontradable nu
		nu = final_expenditure_shares .* (p_sectoral ./ (data["pwt"] .* P_US)) .^ (sigma-1)
		@assert any(isnan, nu[:,:,1:end-1,:]) == false

		nontradable_nu = 1 .- sum(final_expenditure_shares[:,:,1:end-1,:], dims=3)
		# Replace negative elements with second smallest positive
		nu[:,:,end:end,:] .= DetrendUtilities.winsorize(nontradable_nu, 1)
		nu .= nu ./ sum(nu, dims=3)

		# demand shifter only varies across sectors and over time, not across countries
		parameters[:nu_njt] = sum(country_weights .* nu, dims=(1,2))

		# enforce comformity of model with data
		final_expenditure_shares = parameters[:nu_njt] .* (p_sectoral ./ (data["pwt"] .* P_US)) .^ (1-sigma)
		nontradable_nu = 1 .- sum(final_expenditure_shares[:,:,1:end-1,:], dims=3)
		# Replace negative elements with second smallest positive
		final_expenditure_shares[:,:,end:end,:] .= DetrendUtilities.winsorize(nontradable_nu, 1)
		final_expenditure_shares .= final_expenditure_shares ./ sum(final_expenditure_shares, dims=3)

		# step 4: calculate nontradable prices
		# NB: DO NOT recalibrate tradable prices, expenditure_shares are very noisy for small sectors
		if abs(sigma-1)>0.01
			p_sectoral[:,:,end:end,:] = data["pwt"] .* P_US .* (parameters[:nu_njt][:,:,end:end,:] ./ final_expenditure_shares[:,:,end:end,:]) .^ (1/(sigma-1))
		else
			p_sectoral[:,:,end:end,:] = (data["pwt"] .* P_US ./ (prod(p_sectoral[:,:,1:end-1,:] .^ parameters[:nu_njt][:,:,1:end-1,:], dims=3))) .^ (1 ./ parameters[:nu_njt][:,:,end:end,:])
		end
		parameters[:p_sectoral] = p_sectoral
		parameters[:final_expenditure_shares] = final_expenditure_shares
	end

	function calculate_A(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		p_njt = parameters[:p_sectoral]
		p_mjt = array_transpose(p_njt)
		beta_j = parameters[:beta_j]
		gamma = parameters[:gamma_jk]
		kappa_mnjt = parameters[:kappa_mnjt]
		w_njt = parameters[:w_njt]
		B = parameters[:B_j]
		d_mnjt = parameters[:d]
		xi = parameters[:xi]
		theta = parameters[:theta]

		z = zeros(1, N, J, T)

		# use eq 15 in algorithm.pdf
		rho_mnjt = kappa_mnjt .* p_mjt .* d_mnjt .^ (-1/theta)
		rho_njt = exp.(mean(log.(rho_mnjt), dims=1))
		# nontradable input price equals output price
		rho_njt[1,:,end,:] = p_njt[1,:,end,:]
		input_price_index = exp.(rotate_sectors(gamma', log.(p_njt)))

		A_njt = xi .* B ./ rho_njt .* w_njt .^ beta_j .* input_price_index
		return A_njt
	end

	function calculate_expenditure_shares(parameters, data)
		N, J, T = parameters[:N], parameters[:J], parameters[:T]

		nulla = parameters[:numerical_zero]
		weights = parameters[:bp_weights]

		beta = parameters[:beta_j]
		gamma = parameters[:gamma_jk]
		d = parameters[:d]
		# service import shares are NaN
		for t=1:T
			d[:,:,J,t] = eye(N)
		end
		va = data["va"]

		#beta = dropdims(beta,(1,2,4))
		revenue = va ./ beta
		expenditure = zeros(size(revenue))
		for j=1:J
			for t=1:T
				expenditure[1,:,j,t]  = revenue[1,:,j,t]' * inv(d[:,:,j,t])
			end
		end
		intermediate = rotate_sectors(gamma, revenue)
		final_expenditure = expenditure .- intermediate

		# Smooth the series
		_, nu_guess = DetrendUtilities.detrend(final_expenditure, weights)

		# Replace negative elements with smallest positive
		nu_guess .= DetrendUtilities.winsorize(nu_guess ./ sum(nu_guess, dims=3), 0)
		# Smooth the series
		nu_c, nu_t = DetrendUtilities.detrend(nu_guess, weights)

		# Normalization
		return nu_t ./ sum(nu_t, dims=3)	
end

	function estimate_AR1(data)
		# data is M,N,J,T
		_, N, J, T = size(data)
		current = data[:,:,:,2:T]
		lag = data[:,:,:,1:T-1]

		constant = zeros(1, N, J, 1)
		rho = zeros(1, N, J, 1)
		sigma = zeros(1, N, J, 1)

		# estimate a separate AR(1) for each series
		for n=1:N
			for j=1:J
				y = current[1,n,j,:]
				X = cat(ones(T-1), lag[1,n,j,:], dims=2)

				constant[1,n,j,1], rho[1,n,j,1] = X \ y
				sigma[1,n,j,1] = std(y .- X * [constant[1,n,j,1], rho[1,n,j,1]])
			end
		end

		return (constant, rho, sigma)
	end

	function draw_productivity_shocks!(parameters)
		S, T = parameters[:S], parameters[:T]

		parameters[:global_sectoral_shock_njs] = draw_random_realizations(parameters[:global_sectoral_shock], S)
		parameters[:country_shock_njs] = draw_random_realizations(parameters[:country_shock], S)
		parameters[:idiosyncratic_shock_njs] = draw_random_realizations(parameters[:idiosyncratic_shock], S)

		parameters[:A_njs] = map(t ->
			exp.(ImpvolEquilibrium.non_random_variable(parameters[:productivity_trend], t)
				.+ parameters[:global_sectoral_shock_njs][t]
				.+ parameters[:country_shock_njs][t]
				.+ parameters[:idiosyncratic_shock_njs][t]),
				 1:T)
	end

	function draw_random_realizations(data, S)
		# data is M,N,J,T
		_, N, J, T = size(data)
		constant, rho, sigma = estimate_AR1(data)

		draws = Array{Array{Float64, 4}}(undef, T)
		draws[1] = ImpvolEquilibrium.non_random_variable(data, 1)
		for t=2:T
			innovation = sigma .* randn(1,N,J,S - 1)
			random_realization = ImpvolEquilibrium.non_random_variable(data, t)
			past_productivity = ImpvolEquilibrium.non_random_variable(data, t-1)
			# reversion towards mean
			draws[t] = cat(random_realization, constant .* (1 .- rho) .+ past_productivity .* rho .+ innovation, dims=4)
		end
		return draws
	end

	function decompose_shocks!(parameters, sectoral_weights)
		# M,N,J,T
		# Smooth the series
		weights = parameters[:bp_weights]
		detrended_log_productivity, parameters[:productivity_trend] = DetrendUtilities.detrend(log.(parameters[:A]), weights)

		global_sectoral_shock = mean(detrended_log_productivity, dims=2)
		# weighted by sector importance, see https://github.com/ceumicrodata/impvol/commit/91d92905678df96d7068b8dd729e6f6d7cf470d8
		country_shock = sum(sectoral_weights .* (detrended_log_productivity .- global_sectoral_shock), dims=3) ./ sum(sectoral_weights, dims=3)
		idiosyncratic_shock = detrended_log_productivity .- global_sectoral_shock .- country_shock

		parameters[:global_sectoral_shock] = global_sectoral_shock
		parameters[:country_shock] = country_shock
		parameters[:idiosyncratic_shock] = idiosyncratic_shock
	end
	function jld_saver(data, file_name="results.jld2")
		jldopen(file_name, "w") do file
	   		file["results"] = data
		end
	end
	function jld_loader(file_name="results.jld2")
		file = jldopen(file_name, "r")
		return file["results"]
	end
end
