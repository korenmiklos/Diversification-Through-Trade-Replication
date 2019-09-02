using FileIO
using CSV
using DataFrames

function manipulate_import_shares(data::DataFrames.DataFrame, dim_size::NTuple)
	# This function extends the imported import shares, so that it will be full dimensional
	# Adds service sector (with 0 import share) and adds own country import share (0)
	insert!(data, length(names(data))  + 1, 0, :s)
	N = dim_size[2]
	J = dim_size[4] + 1
	T = dim_size[3]
	for t in 1:T
		y = 1971 + t
		for n in 1:N
			vec = zeros(1,J + 3)
			vec[1,1:3] = [y,n,n]
			push!(data, vec)
		end
	end
	return sort!(data,(1,2,3)), map(+,dim_size,(1,0,0,1))
end

function read_data(relativepath::String, dim_size::NTuple, in_pos::Array{Bool,1}, out_pos::Array{Int64,1}, delim::Char, header::Bool, drop::Int64, indic::Bool, convfloat::Bool)
	# relativepath - gives the relative path of the data to be imported
	# dim_size     - gives the size of each dimension the data should be stored in julia,
	#                their product should correspond to the number of datapoints in the dataset
	#                to be imported, order should be appropriate
	# in_pos       - gives existing dimensions in the input data matrix 
	#                corresponding to mnjt ordering notation as expected (but not restricted) in the output,
	#                its length nevertheless defines the dimensionality of the output matrix
	# out_pos      - gives the position of each dimension in the output matrix
	#                as in the order of argument 'dim_size',
	#                the order prefarabely corresponds to mnjt
	#                m - importer country, n - exporter country, j - sector, t - period
	# delim        - delimiter of the CSV.read command
	# header       - header of the CSV.read command
	# drop         - drop the first 'drop' number of columns in the input data matrix
	# indic        - indicator for reading import shares
	#
	# example: read_data("data/raw_imputed/beta_panel.txt",(36,25,24),[false,true,true,true],[2,3,1],'\t',true,2,false)

	length(dim_size) == count(in_pos) || error("Each dimension should be present in the input matrix")
	any(x -> x <= length(in_pos), out_pos) || error("The output matrix cannot be more than $(length(in_pos))-dimensional")
	unique(out_pos) == out_pos || error("Each output dimension should be unique")

	length(dim_size) == length(out_pos) || error("Each dimension should have its own output position and vice versa")

	data = CSV.read(relativepath, header = header, delim = delim)

	# Part purely for manipulating import shares
	if indic
		data, dim_size = manipulate_import_shares(data,dim_size)
	end
	# End of manipulation

	data = convert(Array,data)
	if drop > 0
		data = data[:,(drop + 1):end]
	end
	prod(size(data)) == prod(dim_size) || error("The number of datapoints ($(prod(size(data)))) should match the product of elements in argument 'dim_size' ($(prod(dim_size)))")

	data = reshape(data, dim_size)

	pos = zeros(Int64,size(in_pos))
	pos[in_pos] = out_pos
	if length(in_pos) > length(dim_size)
		pos[.!in_pos] = convert(Array{Int64,1},length(in_pos):-1:(length(dim_size) + 1))
	end
	pos = convert(Array{Int64,1},pos)
	data = permutedims(cat(ndims(data) + (length(in_pos) - length(dim_size)), data), pos)

	if convfloat
		data = collect(Missings.replace(data,NaN))
		data = convert(Array{Float64,length(in_pos)},data)
	end

	return data
end

country_names   = read_data("data/country_name.txt",(25,),[false,true,false,false],[1],'\t',false,0,false,false)
beta            = read_data("data/beta_panel.txt",(36,25,24),[false,true,true,true],[2,3,1],'\t',true,2,false,true)                
pwt             = read_data("data/aggregate_price_relative_to_US.csv",(36,25),[false,true,false,true],[2,1],',',false,2,false,true)
va              = read_data("data/sectoral_value_added.csv",(36,25,24),[false,true,true,true],[2,3,1],',',true,2,false,true)
import_shares   = read_data("data/import_share.txt",(24,25,36,23),[true,true,true,true],[2,1,4,3],'\t',false,3,true,true)
io_values       = read_data("data/oecd_io_values.csv",(34,34,13),[true,true,false,true],[2,1,3],',',true,3,false,true)
total_output    = read_data("data/oecd_total_output.csv",(34,13),[false,true,false,true],[1,2],',',true,2,false,true)
output_shares   = read_data("data/output_shares.csv",(13,10),[false,true,false,true],[2,1],',',true,1,false,true)
intermediate_input_shares = read_data("data/intermediate_input_shares.csv",(13,10),[false,true,false,true],[2,1],',',true,1,false,true)
trade_balance   = read_data("data/trade_balance_new.csv",(25,36),[false,true,false,true],[1,2],',',true,1,false,true)
p_sectoral_data = read_data("data/sectoral_price_index.csv",(36,18,24),[false,true,true,true],[2,3,1],',',false,0,false,true)

save("data/impvol_data.jld2", "beta", beta, "pwt", pwt, "va", va, "import_shares", import_shares, "io_values", io_values, "total_output", total_output, "output_shares", output_shares, "intermediate_input_shares", intermediate_input_shares, "trade_balance", trade_balance, "p_sectoral_data", p_sectoral_data)
