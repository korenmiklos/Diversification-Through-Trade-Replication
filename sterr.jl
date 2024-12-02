using CSV
using DataFrames
using Statistics

seeds = [
    7018
    1201
    6503
    5695
    1499
]

function readfile(seed)
    CSV.read("output/$seed.csv", DataFrame)
end

function extract(dfs, column)
    K = length(dfs)
    N = nrow(dfs[1])
    out = zeros(N, K)
    for k = 1:K
        out[:, k] = dfs[k][:, column]
    end
    out
end

stats(A) = [mean(A, dims=2) std(A, dims=2)]

function main()
    published = readfile("table1-published")
    dfs = readfile.(seeds)
    trade_barriers = extract(dfs, :trade_barriers) |> stats
    diversification = extract(dfs, :diversification) |> stats

    df = DataFrame(country_names=dfs[1].country_names, 
        trade_barriers_mean=trade_barriers[:, 1], 
        trade_barriers_se=trade_barriers[:, 2],
        diversification_mean=diversification[:, 1], 
        diversification_se=diversification[:, 2],
    )
    a = innerjoin(df, 
        published[!, [:country_names, :trade_barriers, :diversification]], 
        on=:country_names)
    a[!, "trade_barriers_z"] = (a.trade_barriers .- a.trade_barriers_mean) ./ a.trade_barriers_se
    a[!, "diversification_z"] = (a.diversification .- a.diversification_mean) ./ a.diversification_se 
    CSV.write("output/sterr.csv", a)
end

main()