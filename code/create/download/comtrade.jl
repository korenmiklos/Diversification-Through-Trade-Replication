#!/usr/bin/env julia
include("utils.jl")

function download_comtrade_year(year::Int; api_key::String = "")
    output_dir = "input/comtrade/raw/$year"
    mkpath(output_dir)
    
    if isempty(api_key)
        try
            api_key = read_api_key("COMTRADE_API_KEY")
        catch e
            @error "No Comtrade API key found. Please add COMTRADE_API_KEY to .env file"
            return
        end
    end
    
    base_url = "https://comtradeapi.un.org/data/v1/get/C/A"
    
    headers = Dict(
        "Ocp-Apim-Subscription-Key" => api_key
    )
    
    query = Dict(
        "period" => year,
        "reporterCode" => "all",
        "partnerCode" => "all", 
        "flowCode" => "M,X",
        "classificationCode" => "S1",
        "includeDesc" => true
    )
    
    checkpoint_file = joinpath(output_dir, "checkpoint.txt")
    
    if isfile(checkpoint_file)
        @info "Found checkpoint for year $year, resuming..."
    end
    
    max_records = 250000
    offset = 0
    batch = 1
    
    while true
        output_path = joinpath(output_dir, "comtrade_$(year)_batch_$(batch).json")
        
        if isfile(output_path)
            @info "Batch $batch already downloaded, skipping..."
            batch += 1
            offset += max_records
            continue
        end
        
        query["maxRecords"] = max_records
        query["offset"] = offset
        
        try
            @info "Downloading Comtrade data for $year (batch $batch)..."
            
            sleep(1.5)
            
            response = api_request(base_url, headers=headers, query=query)
            
            open(output_path, "w") do f
                JSON3.write(f, response)
            end
            
            record_count = length(get(response, :data, []))
            @info "Downloaded $record_count records for year $year (batch $batch)"
            
            if record_count < max_records
                @info "Completed downloading year $year"
                break
            end
            
            batch += 1
            offset += max_records
            
        catch e
            @error "Failed to download Comtrade data" year=year batch=batch error=e
            
            open(checkpoint_file, "w") do f
                println(f, "year=$year,batch=$batch,offset=$offset")
            end
            
            throw(e)
        end
    end
    
    if isfile(checkpoint_file)
        rm(checkpoint_file)
    end
end

function download_comtrade_all()
    years = 1972:2007
    
    api_key = try
        read_api_key("COMTRADE_API_KEY")
    catch e
        @error "No Comtrade API key found. Please add COMTRADE_API_KEY to .env file"
        return
    end
    
    for year in years
        try
            download_comtrade_year(year, api_key=api_key)
        catch e
            @error "Failed to download Comtrade data for year $year" error=e
            @info "You can resume from year $year"
            break
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 0
        year = parse(Int, ARGS[1])
        @info "Downloading Comtrade data for year $year..."
        download_comtrade_year(year)
    else
        @info "Downloading all Comtrade data (1972-2007)..."
        download_comtrade_all()
    end
    @info "Comtrade download complete"
end