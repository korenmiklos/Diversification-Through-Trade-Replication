#!/usr/bin/env julia
include("utils.jl")

function download_pwt71()
    base_url = "https://www.rug.nl/ggdc/productivity/pwt/pwt-releases/pwt-7.1"
    output_dir = "input/pwt"
    mkpath(output_dir)
    
    urls = Dict(
        "pwt71.zip" => "$base_url/pwt71_11302012version.zip",
        "pwt71_data.csv" => "$base_url/pwt71.csv"
    )
    
    for (filename, url) in urls
        output_path = joinpath(output_dir, filename)
        try
            download_with_retry(url, output_path)
            log_download(url, output_path, "success")
            @info "Downloaded PWT 7.1: $filename"
        catch e
            @warn "Could not download from primary source, trying alternative" filename=filename
            alt_url = "https://www.rug.nl/ggdc/docs/pwt71.xlsx"
            try
                download_with_retry(alt_url, output_path)
                log_download(alt_url, output_path, "success")
                @info "Downloaded PWT 7.1 from alternative source: $filename"
            catch e2
                log_download(url, output_path, "failed")
                @error "Failed to download PWT 7.1" filename=filename error=e2
            end
        end
    end
end

function download_pwt56()
    base_url = "https://www.rug.nl/ggdc/productivity/pwt/pwt-releases/pwt-5.6"
    output_dir = "input/pwt"
    mkpath(output_dir)
    
    url = "$base_url/pwt56.zip"
    output_path = joinpath(output_dir, "pwt56.zip")
    
    try
        download_with_retry(url, output_path)
        log_download(url, output_path, "success")
        @info "Downloaded PWT 5.6"
    catch e
        log_download(url, output_path, "failed")
        @error "Failed to download PWT 5.6" error=e
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Downloading Penn World Table data..."
    download_pwt71()
    download_pwt56()
    @info "PWT download complete"
end