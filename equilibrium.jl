module ImpvolEquilibrium

export period_wrapper, coerce_parameters!, rotate_sectors, CES_price_index, array_transpose, remove_shock!

using Logging, LinearAlgebra, Random, Statistics

include("utils.jl")
parameters = Dict{Symbol, Any}()

# for matrix conformity, store all variables in a 4-dimensional array:
# mnjt: destination, source, sector, time

# per-period random variables are stored as
# mnjs: destination, source, sector, state

function report(A)
	display(A[1,1:2,end-1:end,1])
	sleep(1)
end

function where_in_simplex(labor_share, s)
	l = labor_share[:,:,:,s]
	return (minimum(l), ind2sub(l,indmin(l)))
end

function biggest_gap(x1, x2, s)
	y = (x1 .- x2)[:,:,:,:] .^ 2
	return (maximum(y), ind2sub(y,indmax(y)))
end

function expected_value(y)
	_, _, _, S = size(y)
	if S==1
		return y[:,:,:,1]
	else
		# ignore first realization of S
		return mean(y[:,:,:,2:end], dims=4)
	end
end

function rotate_sectors(A, y)
	# y may be a time series variable or a random variable
	M, N, J, T = size(y)
	K = size(A,1)
	B = zeros(M,N,K,T)
	for m=1:M
		for n=1:N
			for t=1:T
				B[m,n,:,t] .= A * y[m,n,:,t]
			end
		end
	end
	return B
end

function p_distance(p1, p2, theta)
	s1 = p1 .^ (-theta) ./ sum(p1 .^ (-theta), dims=2)
	s2 = p2 .^ (-theta) ./ sum(p2 .^ (-theta), dims=2)
	return share_distance(s1, s2)
end

function share_distance(p1, p2)
	return mean((p1 .- p2) .^ 2) .^ 0.5
end

function deflate_all_nominal_variables!(random_variables, parameters, t)
	compute_price_index!(random_variables, parameters, t)
	US_price_index = random_variables[:P_ns][1:1,end:end,1:1,:]
	for variable in [:R_njs, :rho_njs, :w_njs, :P_ns, :input_price_njs, :P_njs, :E_mjs]
		random_variables[variable] = random_variables[variable] ./ US_price_index
	end
end

function free_trade_sector_shares!(parameters)
	N, J, T = parameters[:N], parameters[:J], parameters[:T]
	gamma_jk = parameters[:gamma_jk]
	beta_j = parameters[:beta_j]
	alpha_jt = parameters[:nu_njt] ./ sum(parameters[:nu_njt], dims=3)

	revenue_shares = zeros(1,1,J,T)
	for t=1:T
		revenue_shares[1,1,:,t] .= eigen_share(alpha_jt[1,1,:,t]*beta_j[1,1,:,1]' .+ gamma_jk)
	end
	parameters[:sector_shares] = revenue_shares
end

function input_price_index!(random_variables, parameters)
	P = random_variables[:P_njs]
	random_variables[:input_price_njs] = exp.(rotate_sectors(parameters[:gamma_jk]', log.(P)))
end

function compute_price!(random_variables, parameters, t)
	theta = parameters[:theta]
	kappa = non_random_variable(parameters[:kappa_mnjt], t)
	rho_njs = random_variables[:rho_njs]

	random_variables[:P_njs] = array_transpose(sum((rho_njs ./ kappa) .^ (-theta), dims=2) .^ (-1/theta))
end

function CES_price_index(alpha, P_njs, sigma)
	return sum(alpha .* P_njs .^ (1-sigma), dims=3) .^ (1/(1-sigma))
end

function compute_price_index!(random_variables, parameters, t)
	nu = non_random_variable(parameters[:nu_njt], t)
	alpha = nu ./ sum(nu, dims=3)
	P_njs = random_variables[:P_njs]
	sigma = parameters[:sigma]

	# Cobb-Douglas is a special case when sigma ~ 1
	random_variables[:P_ns] = CES_price_index(alpha, P_njs, sigma)
end

function free_trade_country_shares!(random_variables, parameters)
	free_trade_sector_shares!(parameters)
	A_njs = random_variables[:A_njs]
	L_njs = random_variables[:L_njs]
	theta = parameters[:theta]
	beta_j = parameters[:beta_j]

	d_njs = A_njs .^ (theta ./ (1 .+ theta*beta_j)) .* L_njs .^ (theta .* beta_j ./ (1 .+ theta*beta_j))
	random_variables[:d_njs_free] = d_njs ./ sum(d_njs, dims=2)
end

function free_trade_labor_shares!(random_variables, parameters, t)
	free_trade_country_shares!(random_variables, parameters)

	d_njs = random_variables[:d_njs_free]
	beta_j = parameters[:beta_j]
	E_jt = non_random_variable(parameters[:sector_shares], t)

	y_njs = beta_j .* d_njs .* E_jt

	random_variables[:L_njs_free] = y_njs ./ sum(y_njs, dims=3)
end

function free_trade_wages!(random_variables, parameters, t)
	free_trade_country_shares!(random_variables, parameters)

	d_njs = random_variables[:d_njs_free]
	L_njs = random_variables[:L_njs]
	beta_j = parameters[:beta_j]
	E_wt = non_random_variable(parameters[:sector_shares], t) .* non_random_variable(parameters[:nominal_world_expenditure], t)

	random_variables[:w_njs] = beta_j .* d_njs .* E_wt ./ L_njs
end

function free_trade_prices!(random_variables, parameters, t)
	free_trade_country_shares!(random_variables, parameters)
	beta_j = parameters[:beta_j]
	theta = parameters[:theta]
	d_njs = random_variables[:d_njs_free]
	w_njs = random_variables[:w_njs]
	xi = parameters[:xi]
	B_j = parameters[:B_j]
	A_njs = random_variables[:A_njs]
	gamma_jk = parameters[:gamma_jk]'

	random_variables[:P_njs] = exp.(rotate_sectors(inv(eye(size(gamma_jk, 1)) .- gamma_jk), log.(xi * d_njs .^(1/theta) .* B_j .* (w_njs .^beta_j) ./ A_njs)))
	random_variables[:rho_njs] = random_variables[:P_njs] ./ d_njs .^(1/theta)
	input_price_index!(random_variables, parameters)
end

function compute_trade_shares!(random_variables, parameters, t)
	theta = parameters[:theta]
	kappa = non_random_variable(parameters[:kappa_mnjt], t)
	rho_njs = random_variables[:rho_njs]
	P_mjs = array_transpose(random_variables[:P_njs])

	random_variables[:d_mnjs] = (rho_njs ./ kappa ./ P_mjs) .^ (-theta)
end

function compute_wage!(random_variables, parameters)
	input_price_index!(random_variables, parameters)

	rho_njs = random_variables[:rho_njs]
	input_price_njs = random_variables[:input_price_njs]

	beta_j = parameters[:beta_j]
	B_j = parameters[:B_j]
	xi = parameters[:xi]
	A_njs = random_variables[:A_njs]

	random_variables[:w_njs] = (rho_njs .* A_njs ./ input_price_njs ./ (xi * B_j)) .^ (1 ./ beta_j)
end

function compute_revenue!(random_variables, parameters)
	w_njs = random_variables[:w_njs]
	L_njs = random_variables[:L_njs]
	beta_j = parameters[:beta_j]

	random_variables[:R_njs] = w_njs .* L_njs ./ beta_j
end

function CES_share(nu, price, sigma)
	temp = nu .* price .^ (1-sigma)
	return temp ./ sum(temp, dims=3)
end

function compute_expenditure_shares!(random_variables, parameters, t)
	R_nks = random_variables[:R_njs]
	beta_j = parameters[:beta_j]
	gamma_jk = parameters[:gamma_jk]
	nu = non_random_variable(parameters[:nu_njt], t)
	# encompass CES and Cobb-Douglas
	alpha_njt = CES_share(nu, random_variables[:P_njs], parameters[:sigma])
	S_nt = non_random_variable(parameters[:S_nt], t)
	expenditure = sum(R_nks, dims=3) .- S_nt

	wagebill_ns = rotate_sectors(beta_j[:]', R_nks)
	intermediate_njs = rotate_sectors(gamma_jk, R_nks)
	e_njs = alpha_njt .* wagebill_ns .+ intermediate_njs .- alpha_njt .* S_nt
	# trade imbalance may make it negative
	e_njs = max.(parameters[:numerical_zero], e_njs)
	random_variables[:e_mjs] = array_transpose(e_njs ./ sum(e_njs, dims=3))
end

function compute_real_gdp!(random_variables, parameters, t)
	compute_price_index!(random_variables, parameters, t)
	w_njs = random_variables[:w_njs]
	L_njs = random_variables[:L_njs]
	P_ns = random_variables[:P_ns]

	random_variables[:real_GDP] = w_njs .* L_njs ./ P_ns
end

function fixed_expenditure_shares!(random_variables, parameters, t)
	R_mjs = array_transpose(random_variables[:R_njs])
	expenditure = sum(R_mjs, dims=3) .- array_transpose(non_random_variable(parameters[:S_nt], t))
	random_variables[:E_mjs] = random_variables[:e_mjs] .* max.(parameters[:numerical_zero], expenditure)
end

function shadow_price_step(random_variables, parameters, t)
	theta = parameters[:theta]
	E_mjs = random_variables[:E_mjs]
	R_njs = random_variables[:R_njs]
	kappa_mnjt = non_random_variable(parameters[:kappa_mnjt], t)
	P_mjs = array_transpose(random_variables[:P_njs])

	return sum( (kappa_mnjt .* P_mjs) .^ theta .* E_mjs ./ R_njs, dims=1) .^ (1/theta)
end

function starting_values!(random_variables, parameters, t)
	free_trade_wages!(random_variables, parameters, t)
	free_trade_prices!(random_variables, parameters, t)
	free_trade_labor_shares!(random_variables, parameters, t)
	compute_revenue!(random_variables, parameters)
	compute_expenditure_shares!(random_variables, parameters, t)
	fixed_expenditure_shares!(random_variables, parameters, t)
	deflate_all_nominal_variables!(random_variables, parameters, t)
end

function inner_loop!(random_variables, parameters, t)
	@debug("------ BEGIN Inner loop")
	lambda = parameters[:inner_step_size]
	dist = 999
	k = 1

	while (dist > parameters[:inner_tolerance]) && (k <= parameters[:max_iter_inner])
		new_rho = shadow_price_step(random_variables, parameters, t)
		dist = p_distance(new_rho, random_variables[:rho_njs], parameters[:theta])
		@debug("-------- Inner ", k, ": ", dist)
		random_variables[:rho_njs] = new_rho
		compute_price!(random_variables, parameters, t)
		compute_wage!(random_variables, parameters)
		old_R = random_variables[:R_njs]
		compute_revenue!(random_variables, parameters)
		random_variables[:R_njs] = lambda*random_variables[:R_njs] .+ (1-lambda)*old_R
		fixed_expenditure_shares!(random_variables, parameters, t)
		deflate_all_nominal_variables!(random_variables, parameters, t)

		k = k+1
	end
	@warn("inner: ", k-1)
	@debug("------ END Inner loop")
end

function middle_loop!(random_variables, parameters, t)
	@debug("---- BEGIN Middle loop")
	dist = 999
	k = 1
	old_expenditure_shares = random_variables[:e_mjs]

	while (dist > parameters[:middle_tolerance]) && (k <= parameters[:max_iter_middle])
		inner_loop!(random_variables, parameters, t)
		compute_expenditure_shares!(random_variables, parameters, t)

		dist = share_distance(random_variables[:e_mjs][:,:,1:end-1,:], old_expenditure_shares[:,:,1:end-1,:])

		random_variables[:e_mjs] = parameters[:middle_step_size]*random_variables[:e_mjs]+(1-parameters[:middle_step_size])*old_expenditure_shares
		@info("------ Middle ", k, ": ", dist)

		old_expenditure_shares = random_variables[:e_mjs]
		k = k+1
	end
	@debug("---- END Middle loop")
end

function adjustment_loop!(random_variables, L_nj_star, parameters, t)

	function max_step_size(x0, direction)
		nulla = parameters[:numerical_zero]
		negative_entries = min.(direction, -nulla)
		step_sizes = nulla .- x0 ./ negative_entries
        return minimum(step_sizes, dims=3)
	end

	function evaluate_utility(random_variables, L_nj_star, parameters, t)
		random_variables[:L_njs] = max.(parameters[:numerical_zero], random_variables[:L_njs])
		random_variables[:L_njs] = random_variables[:L_njs] ./ sum(random_variables[:L_njs], dims=3)

		middle_loop!(random_variables, parameters, t)

		w_ns = sum(random_variables[:w_njs] .* random_variables[:L_njs], dims=3)
		wage_gap = random_variables[:w_njs] ./ w_ns
		return sum(parameters[:one_over_rho]*log.(w_ns) .- 0.5*sum((random_variables[:L_njs] .- L_nj_star).^2, dims=3))
	end

	function calculate_derivative(random_variables, L_nj_star, parameters)
		w_njs = random_variables[:w_njs]
		L_njs = random_variables[:L_njs]
		w_ns = sum(w_njs .* random_variables[:L_njs], dims=3)
		wage_gap = w_njs ./ w_ns
		gradient = parameters[:one_over_rho]*wage_gap .- (L_njs .- L_nj_star)

		# ensure that steps sum to zero: we can only reallocate labor across sectors
		return gradient .- mean(gradient, dims=3)
	end

	@debug("-- BEGIN Adjustment loop")
	random_variables[:L_njs] = L_nj_star
	free_trade_labor_shares!(random_variables, parameters, t)
	L_njs_free = random_variables[:L_njs_free]
	k = 1

	nulla = ones(1,1,1,1)*parameters[:numerical_zero]

	utility = evaluate_utility(random_variables, L_nj_star, parameters, t)
	gradient = calculate_derivative(random_variables, L_nj_star, parameters)
	dist = mean(gradient .^ 2) .^ 0.5
	@info("Starting from $dist")

	while (dist > parameters[:adjustment_tolerance]) && (k <= parameters[:max_iter_adjustment])
		snapshot = copy(random_variables)

		previous_utility = utility
		previous_L = random_variables[:L_njs]
		step_size = min.(0.33*max_step_size(previous_L, gradient), 0.1)

		snapshot[:L_njs] = previous_L .+ gradient .* step_size

		utility = evaluate_utility(snapshot, L_nj_star, parameters, t)
		gradient = calculate_derivative(snapshot, L_nj_star, parameters)
		dist = mean(gradient .^ 2) .^ 0.5

		difference = utility - previous_utility
		proportional_increase = difference / sum(gradient .^ 2 .* step_size)
		@debug("Difference: ", difference)
		@debug("Proportional increase: ", proportional_increase)

		@info("---- Adjustment $k: $dist")

		k = k+1
		random_variables = snapshot
	end
	@debug("-- END Adjustment loop")
end

function expected_wage_share(random_variables)
	w_njs = random_variables[:w_njs]
	L_njs = random_variables[:L_njs]

	wage_bill = w_njs .* L_njs
	wage_share = wage_bill ./ sum(wage_bill, dims=3)

	return expected_value(wage_share)
end

function outer_loop!(random_variables, parameters, t, L_nj_star)
	N, J = parameters[:N], parameters[:J]
	lambda = parameters[:outer_step_size]
	random_variables[:L_njs] = copy(L_nj_star)

	@info("BEGIN Outer loop")
	starting_values!(random_variables, parameters, t)

	dist = 999
	k = 1

	while (dist > parameters[:outer_tolerance]) && (k <= parameters[:max_iter_outer])
		old_wage_share = L_nj_star ./ sum(L_nj_star, dims=3)
		adjustment_loop!(random_variables, L_nj_star, parameters, t)
		wage_share = expected_wage_share(random_variables)
		@assert sum(wage_share, dims=3) ≈ ones(1,N,1,1) atol=1e-9

		dist = share_distance(wage_share, old_wage_share)
		@info("-- Outer ", k, ": ", dist)

		L_nj_star = (1-lambda)*old_wage_share .+ lambda*wage_share
		k = k+1
	end
	@info("END Outer loop")
	return L_nj_star
end

function period_wrapper(parameters, t)
	N, J = parameters[:N], parameters[:J]

	@info("--- Period ", t, " ---")
	A_njs = parameters[:A_njs][t]
	random_variables = Dict{Symbol, Any}()
	random_variables[:A_njs] = A_njs

	random_variables[:L_njs] = zeros(1,N,J,1)

	for n=1:N
		# start from expenditure labor weights
		random_variables[:L_njs][1,n,:,1] .= parameters[:importance_weight][:]
	end
	free_trade_labor_shares!(random_variables, parameters, t)
	stv = random_variables[:L_njs_free]
	L_nj_star = outer_loop!(random_variables, parameters, t, stv)

	compute_trade_shares!(random_variables, parameters, t)
	compute_real_gdp!(random_variables, parameters, t)
	return random_variables
end

end

