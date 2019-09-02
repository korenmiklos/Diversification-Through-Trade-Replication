module DetrendUtilities

	function winsorize(data, trim=1)
		M, N, J, T = size(data)

		nu_guess = copy(data)
		for m = 1:M, j = 1:J, t = 1:T
			vector = nu_guess[m,:,j,t]
			smallest = sort(vector[vector .> 0])[1+trim]
			largest = sort(vector[vector .< 1])[end-trim]
			nu_guess[m,vector .< smallest,j,t] = smallest
			nu_guess[m,vector .> largest,j,t] = largest
		end
	    return nu_guess
	end

	function detrend(X, weights::Array{Float64,1})
		(ndims(X) == 2 || ndims(X) == 4) || error("The function 'detrend' can only be applied to 2 or 4-dimensional arrays")
		
		dims = ndims(X)

		if dims == 2
			I, S = size(X)
			J, K = 1, 1
		elseif dims == 4
			I, J, K, S = size(X)
			X = reshape(X, (I*J*K, S))
		end

		X_c = zeros(I*J*K,S)
		X_t = zeros(I*J*K,S)

		# Maximum lag
		P = size(weights,1) - 1
		w = [flipdim(weights,1); weights[2:end]]

		# Extended series with P 0's at both ends
		X_ext = zeros(I*J*K, S + 2P)
		X_ext[:,(P + 1):(P + S)] = X

		# Reflect head/tail
		X_ext[:,1:P] = repeat(X[:,1], outer = [1,P]) + flipdim(repeat(X[:,1], outer = [1,P]) - X[:,2:(P + 1)],2)
		X_ext[:,(S + P + 1):(S + 2P)] = repeat(X[:,S], outer = [1,P]) + flipdim(repeat(X[:,S],outer = [1,P]) - X[:,(S - P):(S - 1)],2)

		for s in 1:S
			X_c[:,s] = X_ext[:,s:(s + 2P)] * w
			X_t[:,s] = X_ext[:,s + P] - X_c[:,s]
		end

		if dims == 4
			X_c = reshape(X_c, (I, J, K, S))
			X_t = reshape(X_t, (I, J, K, S))
		end

		return X_c, X_t
	end
end