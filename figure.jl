include("output.jl")
plt_dta = ImpvolOutput.plot_data("experiments/labor_adjustment/actual/results.jld2") # with the desired arguments if needed

Base.print_matrix(IOContext(STDOUT, :limit => false), plt_dta[3]) # to check the countries' position
idx = zeros(length(plt_dta[3]))
idx[6], idx[10], idx[11], idx[13], idx[15], idx[16], idx[24], idx[25] = 1, 1, 1, 1, 1, 1, 1, 1
idx = idx .== 1
ImpvolOutput.plot_model_vs_data((plt_dta[1][idx,:], plt_dta[2][idx,:], plt_dta[3][idx]), "GDP")
# solid line = data, dashed line = model