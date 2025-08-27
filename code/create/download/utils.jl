using Downloads
using HTTP
using JSON3
using DataFrames
using CSV
using Dates
using Printf

struct DownloadConfig
    cache_dir::String
    max_retries::Int
    initial_delay::Float64
    max_delay::Float64
    timeout::Int
end

DownloadConfig() = DownloadConfig(
    "temp/download_cache",
    3,
    1.0,
    60.0,
    300
)

function ensure_cache_dir(config::DownloadConfig)
    mkpath(config.cache_dir)
end

function get_cache_path(config::DownloadConfig, url::String)
    hash_str = string(hash(url))
    joinpath(config.cache_dir, "cache_$hash_str.tmp")
end

function download_with_retry(url::String, output_path::String; config::DownloadConfig = DownloadConfig())
    ensure_cache_dir(config)
    cache_path = get_cache_path(config, url)
    
    if isfile(cache_path)
        @info "Using cached file for $url"
        cp(cache_path, output_path, force=true)
        return output_path
    end
    
    delay = config.initial_delay
    
    for attempt in 1:config.max_retries
        try
            @info "Download attempt $attempt for $url"
            Downloads.download(url, output_path, timeout=config.timeout)
            cp(output_path, cache_path, force=true)
            @info "Successfully downloaded $url"
            return output_path
        catch e
            if attempt == config.max_retries
                @error "Failed to download after $config.max_retries attempts" url=url error=e
                throw(e)
            else
                @warn "Download failed, retrying..." attempt=attempt delay=delay error=e
                sleep(delay)
                delay = min(delay * 2, config.max_delay)
            end
        end
    end
end

function api_request(url::String; 
                    method::String = "GET",
                    headers::Dict = Dict(),
                    query::Dict = Dict(),
                    config::DownloadConfig = DownloadConfig())
    
    delay = config.initial_delay
    
    for attempt in 1:config.max_retries
        try
            response = HTTP.request(
                method,
                url,
                headers,
                query=query,
                timeout=config.timeout
            )
            
            if response.status == 200
                return JSON3.read(String(response.body))
            else
                error("API request failed with status $(response.status)")
            end
        catch e
            if attempt == config.max_retries
                @error "API request failed after $config.max_retries attempts" url=url error=e
                throw(e)
            else
                @warn "API request failed, retrying..." attempt=attempt delay=delay error=e
                sleep(delay)
                delay = min(delay * 2, config.max_delay)
            end
        end
    end
end

function read_api_key(key_name::String)
    env_file = ".env"
    if !isfile(env_file)
        error("No .env file found. Please create one with API keys.")
    end
    
    for line in eachline(env_file)
        if startswith(line, "$key_name=")
            return strip(split(line, "=", limit=2)[2])
        end
    end
    
    error("API key $key_name not found in .env file")
end

function validate_download(file_path::String, expected_size::Union{Int, Nothing} = nothing)
    if !isfile(file_path)
        error("Downloaded file does not exist: $file_path")
    end
    
    actual_size = filesize(file_path)
    if actual_size == 0
        error("Downloaded file is empty: $file_path")
    end
    
    if expected_size !== nothing && actual_size != expected_size
        @warn "File size mismatch" expected=expected_size actual=actual_size
    end
    
    return true
end

function log_download(url::String, output_path::String, status::String)
    log_file = "temp/download_log.csv"
    timestamp = Dates.now()
    
    if !isfile(log_file)
        df = DataFrame(
            timestamp = DateTime[],
            url = String[],
            output_path = String[],
            status = String[]
        )
    else
        df = CSV.read(log_file, DataFrame)
    end
    
    push!(df, (timestamp, url, output_path, status))
    CSV.write(log_file, df)
end