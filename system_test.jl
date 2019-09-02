include("equilibrium.jl")
using ImpvolEquilibrium
include("calibrate_params.jl")
using .CalibrateParameters
using Base.Test
using Logging
Logging.configure(level=INFO)

function test_data!(parameters, data)
    N, J, T, S = parameters[:N], parameters[:J], parameters[:T], parameters[:S]
    parameters[:beta_j] = zeros(1,1,J,1)
    parameters[:beta_j][1,1,:,1] = 0.5

    parameters[:gamma_jk] = [0.25 0.125 0.125; 0.125 0.25 0.125; 0.125 0.125 0.25]
    @test size(parameters[:gamma_jk]) == (J,J)
    C = diagm(parameters[:beta_j][1,1,:,1])+parameters[:gamma_jk]'-eye(J)
    @test sum(C, 2) ≈ zeros(J,1) atol=1e-9

    parameters[:S_nt] = zeros(1,N,1,T)
    parameters[:S_nt_data] = zeros(parameters[:S_nt])
    @test size(parameters[:S_nt_data]) == (1,N,1,T)
    @test sum(parameters[:S_nt_data], 2) ≈ zeros(1,1,1,T) atol=1e-9

    parameters[:d] = zeros(N,N,J,T)
    parameters[:d][:,:,1,1] = [0.9 0.1; 0.3 0.7]
    parameters[:d][:,:,2,1] = [0.7 0.3; 0.1 0.9]
    parameters[:d][:,:,3,1] = [1.0 0.0; 0.0 1.0]
    @test size(parameters[:d]) == (N,N,J,T)
    @test sum(parameters[:d], 2) ≈ ones(N,1,J,T) atol=1e-9

    parameters[:kappa_mnjt] = CalibrateParameters.trade_costs(parameters)
    @test size(parameters[:kappa_mnjt]) == (N,N,J,T)

    final_expenditure_shares = zeros(1,1,J,1)
    final_expenditure_shares[1,1,:,1] = [0.3 0.1 0.6]
    @test sum(final_expenditure_shares, 3) ≈ ones(1,1,1,1) atol=1e-9

    # broad country weights for final expenditure
    country_weights = ones(1,N,1,1)
    country_weights = country_weights ./ sum(country_weights, 2)

    CalibrateParameters.calculate_p_and_nu!(parameters, data, final_expenditure_shares, country_weights)
    @test size(parameters[:nu_njt]) == (1,1,J,T)
    @test sum(parameters[:nu_njt], 3) ≈ ones(1,1,1,T) atol=1e-9
    @test size(parameters[:p_sectoral]) == (1,N,J,T)
    @test any(isnan.(parameters[:p_sectoral])) == false
    @test parameters[:p_sectoral][1,end,:,1] ≈ ones(J) atol=1e-9

    display(parameters[:p_sectoral][1,:,:,1])
    parameters[:w_njt] = CalibrateParameters.calculate_nominal_wages(parameters, data)

    parameters[:B_j] = CalibrateParameters.calculate_B(parameters)
    @test size(parameters[:B_j]) == (1,1,J,1)

    parameters[:xi] = CalibrateParameters.calculate_xi(parameters)
    @test typeof(parameters[:xi]) == Float64

    parameters[:A] = CalibrateParameters.calculate_A(parameters, data)
    @test size(parameters[:A]) == (1, N, J, T)
    display(parameters[:A][1,:,:,1])

    info("US wage rate: ", sum(data["va"], 3)[1,end,1,1])

    # total world expenditure in the data - needed to get reasonable starting values
    parameters[:nominal_world_expenditure] = sum(data["va"] ./ parameters[:beta_j], (1,2,3))
    @test size(parameters[:nominal_world_expenditure]) == (1,1,1,T)

    # global, all-time average of sector final expenditure shares
    importance_weight = mean(parameters[:nu_njt], (1, 2, 4))
    parameters[:importance_weight] = importance_weight
    parameters[:A_njs] = Array{Array{Float64, 4}}(T)
    parameters[:A_njs][1] = parameters[:A] .* exp.(0.2*randn(1,N,J,S))
    @test size(parameters[:A_njs]) == (T,)
    @test size(parameters[:A_njs][1]) == (1,N,J,S)

end

function init_parameters()
    parameters = Dict{Symbol, Any}()
    N, J, T = 2, 3, 1
    parameters[:N], parameters[:J], parameters[:T] = N, J, T
    # CES parameters
	parameters[:sigma] = 0.999
	parameters[:theta] = 4.0
	parameters[:eta] = 4.0

	########## parameters common across scenarios
	## these are function of data
	# inverse of adjustment cost, 0 if cannot readjust
	parameters[:one_over_rho] = 0.1
    parameters[:S] = 100

    parameters[:numerical_zero] = 1e-12

    parameters[:bp_weights] = [0.774074394803123; -0.201004684236153; -0.135080548288772; -0.0509519648766360]

    # adaptive step size. large lambda means large steps
    parameters[:inner_step_size] = exp(-0.10*(parameters[:J]-1)^0.75)
    # large substitution needs more dampening
    parameters[:middle_step_size] = exp(-0.275*max(1,parameters[:sigma]))
    parameters[:adjustment_step_size] = 0.25
    # any deviation from sigma=1 needs more dampening
    parameters[:outer_step_size] = exp(-0.5*abs(log(parameters[:sigma])))
    # this is log points of average input price differences
    parameters[:inner_tolerance]  = 0.001
    parameters[:middle_tolerance] = 0.001
    parameters[:adjustment_tolerance] = 0.0005
    parameters[:outer_tolerance] = 0.001

    # maximum number of iterations in each loop
    parameters[:max_iter_inner] = 1000
    parameters[:max_iter_middle] = 50
    parameters[:max_iter_adjustment] = 50
    parameters[:max_iter_outer] = 50
    return parameters
end

function init_data(parameters)
    N, J, T = parameters[:N], parameters[:J], parameters[:T]
    data = Dict{String,Any}()

    data["p_sectoral_data"] = ones(1,N,J,T)
    data["pwt"] = ones(1,N,1,T)
    data["va"] = ones(1,N,J,T)

    return data
end

srand(7094)

parameters = init_parameters()
data = init_data(parameters)
test_data!(parameters, data)

results = ImpvolEquilibrium.period_wrapper(parameters, 1)
