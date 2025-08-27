#!/usr/bin/env julia
include("utils.jl")

function download_un_national_accounts()
    base_url = "http://data.un.org/ws/rest/data"
    output_dir = "input/unna"
    mkpath(output_dir)
    
    dataset_code = "SNA"
    series_code = "201"
    
    query_params = Dict(
        "dataset" => dataset_code,
        "series" => series_code,
        "format" => "csv"
    )
    
    countries = [
        "USA", "CAN", "MEX", "BRA", "ARG", "CHL", "COL", "PER", "VEN",
        "GBR", "DEU", "FRA", "ITA", "ESP", "NLD", "BEL", "CHE", "SWE", "NOR", "DNK", "FIN",
        "JPN", "CHN", "IND", "KOR", "IDN", "THA", "MYS", "PHL", "SGP", "VNM",
        "AUS", "NZL", "ZAF", "EGY", "TUR", "ISR", "SAU", "ARE", "RUS", "POL"
    ]
    
    for country in countries
        output_path = joinpath(output_dir, "sna_201_$country.csv")
        
        url = "$base_url/$dataset_code/$series_code/all/$country/all"
        
        try
            @info "Downloading UN SNA data for $country"
            download_with_retry(url, output_path)
            log_download(url, output_path, "success")
        catch e
            @warn "Failed to download UN SNA data for $country, trying alternative approach" error=e
            
            alt_url = "http://data.un.org/Export.aspx?d=$dataset_code&f=group_code:$series_code;country_code:$country&c=2,3,4,5,6&s=country_name:asc,year:desc&v=1"
            try
                download_with_retry(alt_url, output_path)
                log_download(alt_url, output_path, "success")
            catch e2
                log_download(url, output_path, "failed")
                @error "Failed to download UN SNA data" country=country error=e2
            end
        end
    end
    
    @info "Downloaded UN National Accounts data"
end

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Downloading UN National Accounts data..."
    download_un_national_accounts()
    @info "UN National Accounts download complete"
end