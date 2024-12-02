# 2024-09-28
## Look for major changes in code between 0.6 and 1.10
Version of record is 3a882eab0979a7dedee5804c7aa9625ae25241de

- Commit 6888ba19e31dcd0ca9464ff084270096b0f29be7 has a lot of reshapes that can potentially change functionality. 
- `read_data.jl` is also changed in f2781b14a0d9687372f333fe422d95f117ee376a
- `p_sectoral` changed in 8bd5838963b15bd713a2c270ff0d1d7a7ee77769
- Other changes seem like API changes only.

## Check 8bd5838963b15bd713a2c270ff0d1d7a7ee77769
### Old code
```julia
p_sectoral = array_transpose(exp.( mean(1 / theta * log.(d ./ permutedims(cat(ndims(d),d[end,:,:,:]),[4,1,2,3])) .- log.(kappa ./ permutedims(cat(ndims(kappa),kappa[end,:,:,:]),[4,1,2,3])), dims=2) .+ repeat(permutedims(cat(ndims(p_sectoral_base), log.(p_sectoral_base[:,end,:,:])), [1,4,2,3]), outer = [size(d,1),1,1,1]) ))
```

Prettify:
```julia
p_sectoral = array_transpose(
    exp.( 
        mean(
            1 / theta * log.(d ./ permutedims(
                    cat(ndims(d),
                        d[end,:,:,:]
                    ),
                    [4,1,2,3]
                )) .- log.(kappa ./ permutedims(
                    cat(ndims(kappa),
                        kappa[end,:,:,:]
                    ),
                    [4,1,2,3]
                )), 
        dims=2) 
        .+ 
        repeat(
            permutedims(
                cat(ndims(p_sectoral_base), 
                    log.(p_sectoral_base[:,end,:,:])
                    ), 
                    [1,4,2,3]), 
        outer = [size(d,1),1,1,1]) 
    )
)
```

### New code
```julia
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
```

> In the new code, `1/theta` applies to both d and kappa. In the old code, only d.

## Check this