macro assert0(expression)
	:(@assert all($expression .<= 1e-12))
end

macro unpack_dictionary(expression)
	lines = []
	dict = eval(expression)
	for (key, value) in dict
		:(local $key = $value)
	end 
end

function eye(n)
	return Matrix{Float64}(I, n, n)
end

function fill_dict(;kwargs...)
	return Dict(kwargs)
end

function fill_dict!(dict; kwargs...)
	merge!(dict, Dict(kwargs))
end

function array_transpose(array)
	index = [i for i in 1:ndims(array)]
	index[1], index[2] = 2, 1
	return permutedims(array, index)	
end

function rotate_along_dimension(rotator, array, n)
	index = [i for i in 1:ndims(array)]
	index[1], index[n] = n, 1
	temp = permutedims(array, index)
	S1 = size(temp)
	A = zeros()
	for i in CartesianRange(S1[2:end])
		@show i
		A[:,i] = rotator*temp[:,i]
	end
	return A
end

function share_within_dimension(X, i)
	return X ./ sum(X,i)
end

function eigen_share(A)
	# the solution to x = Ax with x'1 = 1
	D, V = eigen(A)
	i = (1:length(D))[abs.(D .- 1.0) .< 1e-6][1]
	return V[:,i] ./ sum(V[:,i])
end

function non_random_variable(y, t)
	# array coersion is going to take care of rest
	B = y[:,:,:,t]
	return cat(B, dims=ndims(B)+1)
end

function remove_shock!(parameters, symbol)
	for t=1:length(parameters[:A_njs])
		shock_to_remove = exp.(parameters[symbol][t])
		parameters[:A_njs][t] = parameters[:A_njs][t] ./ shock_to_remove
	end
end
