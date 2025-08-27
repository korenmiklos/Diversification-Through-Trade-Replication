#!/usr/bin/env julia
include("utils.jl")

function download_wdi_archive()
    output_dir = "input/wdi"
    mkpath(output_dir)
    
    urls = [
        "http://databank.worldbank.org/data/download/archive/WDI_excel_2015_10.zip",
        "http://databank.worldbank.org/data/download/archive/WDI_csv_2015_10.zip"
    ]
    
    for (idx, url) in enumerate(urls)
        format = idx == 1 ? "excel" : "csv"
        output_path = joinpath(output_dir, "WDI_$(format)_2015_10.zip")
        
        try
            @info "Downloading WDI October 2015 archive ($format format)..."
            download_with_retry(url, output_path)
            log_download(url, output_path, "success")
            @info "Downloaded WDI archive: $format format"
            break
        catch e
            log_download(url, output_path, "failed")
            if idx == length(urls)
                @error "Failed to download WDI archive from all sources" error=e
            else
                @warn "Failed to download $format format, trying alternative..." error=e
            end
        end
    end
end

function download_wdi_api_data()
    output_dir = "input/wdi"
    mkpath(output_dir)
    
    indicators = [
        "NY.GDP.MKTP.KN",
        "NY.GDP.MKTP.CN", 
        "FP.CPI.TOTL.ZG",
        "FP.CPI.TOTL",
        "PA.NUS.FCRF"
    ]
    
    base_url = "http://api.worldbank.org/v2/country/all/indicator"
    
    for indicator in indicators
        output_path = joinpath(output_dir, "wdi_$indicator.json")
        url = "$base_url/$indicator?date=1970:2007&format=json&per_page=20000"
        
        try
            @info "Downloading WDI indicator: $indicator"
            download_with_retry(url, output_path)
            log_download(url, output_path, "success")
        catch e
            log_download(url, output_path, "failed")
            @warn "Failed to download WDI indicator via API" indicator=indicator error=e
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Downloading World Development Indicators..."
    download_wdi_archive()
    download_wdi_api_data()
    @info "WDI download complete"
end