include("equilibrium.jl")
include("calibrate_params.jl")
using CalibrateParameters
using Logging

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

function test_parameters(parameters)
    function positive(args...)
        for item in args
            @assert minimum(item) > 0.0
        end
    end
    @assertset "CalibrateParameters tests" begin
    N, J, T = parameters[:N], parameters[:J], parameters[:T] 
    @assert size(parameters[:beta_j]) == (1,1,J,1)
    @assert size(parameters[:gamma_jk]) == (J,J)

    C = diagm(parameters[:beta_j][1,1,:,1])+parameters[:gamma_jk]'-eye(J)
    @assert sum(C, 2) ≈ zeros(J,1) atol=1e-9

    @assert size(parameters[:S_nt_data]) == (1,N,1,T)
    @assert sum(parameters[:S_nt_data], 2) ≈ zeros(1,1,1,T) atol=1e-9

    @assert size(parameters[:d]) == (N,N,J,T)
    @assert sum(parameters[:d], 2)[:,:,1:end-1,:] ≈ ones(N,1,J-1,T) atol=1e-9

    @assert size(parameters[:kappa_mnjt]) == (N,N,J,T)

    @assert size(parameters[:final_expenditure_shares]) == (1,N,J,T)
    @assert sum(parameters[:final_expenditure_shares], 3) ≈ ones(1,N,1,T) atol=1e-9

    @assert size(parameters[:nu_njt]) == (1,1,J,T)
    @assert sum(parameters[:nu_njt], 3) ≈ ones(1,1,1,T) atol=1e-9
    @assert size(parameters[:p_sectoral]) == (1,N,J,T)
    @assert any(isnan.(parameters[:p_sectoral])) == false
    @assert parameters[:p_sectoral][1,end,:,1] ≈ ones(J) atol=1e-9

    @assert size(parameters[:w_njt]) == (1,N,J,T)

    @assert size(parameters[:B_j]) == (1,1,J,1)

    @assert typeof(parameters[:xi]) == Float64

    @assert size(parameters[:A]) == (1, N, J, T)
    # productivity does not change more than 50% per year
    logA = log.(parameters[:A])
    @assert maximum(abs.(logA[:,:,:,T] .- logA[:,:,:,1])) < log(1.50)*T

    @assert size(parameters[:nominal_world_expenditure]) == (1,1,1,T)

    @assert size(parameters[:A_njs]) == (T,)
    @assert size(parameters[:A_njs][1]) == (1,N,J,1)
    @assert size(parameters[:A_njs][2]) == (1,N,J,parameters[:S])

    for t=1:T
        @assert parameters[:A_njs][t][:,:,:,1] ≈ parameters[:A][:,:,:,t]
    end

    positive(parameters[:final_expenditure_shares], parameters[:nu_njt], parameters[:A_njs][2], parameters[:p_sectoral], parameters[:A])
    end
end

Random.seed!(7094)

parameters = init_parameters()
CalibrateParameters.calibrate_parameters!(parameters, "data/impvol_data.jld2")
test_parameters(parameters)