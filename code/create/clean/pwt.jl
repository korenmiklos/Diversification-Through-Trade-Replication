#!/usr/bin/env julia
using DataFrames
using CSV
using ZipFile

function extract_pwt_zip(version::String)
    zip_path = "input/pwt/pwt$version.zip"
    extract_dir = "temp/pwt$version"
    
    if !isfile(zip_path)
        @warn "PWT $version zip file not found, trying CSV directly"
        return nothing
    end
    
    mkpath(extract_dir)
    
    reader = ZipFile.Reader(zip_path)
    for file in reader.files
        output_path = joinpath(extract_dir, file.name)
        mkpath(dirname(output_path))
        open(output_path, "w") do out
            write(out, read(file))
        end
    end
    close(reader)
    
    return extract_dir
end

function clean_pwt71()
    extract_dir = extract_pwt_zip("71")
    
    if extract_dir === nothing
        csv_path = "input/pwt/pwt71_data.csv"
        if !isfile(csv_path)
            error("No PWT 7.1 data found")
        end
        df = CSV.read(csv_path, DataFrame)
    else
        csv_files = filter(f -> endswith(f, ".csv"), readdir(extract_dir, join=true))
        if isempty(csv_files)
            xls_files = filter(f -> endswith(f, ".xls") || endswith(f, ".xlsx"), readdir(extract_dir, join=true))
            if !isempty(xls_files)
                @warn "Found Excel files, need to convert to CSV first"
                error("Excel conversion not implemented yet")
            end
        end
        df = CSV.read(csv_files[1], DataFrame)
    end
    
    required_cols = [:country, :isocode, :year, :pop, :rgdpch, :rgdpl, :p]
    
    for col in required_cols
        if !(col in names(df, Symbol))
            @warn "Missing column in PWT 7.1: $col"
        end
    end
    
    df_clean = select(df, 
        :country => :country_name,
        :isocode => :iso3,
        :year => :year,
        :pop => :population,
        :rgdpch => :gdp_per_capita,
        :rgdpl => :gdp_laspeyres,
        :p => :price_level
    )
    
    df_clean = dropmissing(df_clean, [:year, :iso3])
    
    df_clean = filter(row -> 1970 <= row.year <= 2007, df_clean)
    
    output_path = "temp/pwt71_clean.csv"
    CSV.write(output_path, df_clean)
    @info "Cleaned PWT 7.1 data saved to $output_path"
    
    return df_clean
end

function clean_pwt56()
    extract_dir = extract_pwt_zip("56")
    
    if extract_dir === nothing
        @warn "PWT 5.6 not available, skipping"
        return nothing
    end
    
    csv_files = filter(f -> endswith(f, ".csv"), readdir(extract_dir, join=true))
    if isempty(csv_files)
        @warn "No CSV files found in PWT 5.6 archive"
        return nothing
    end
    
    df = CSV.read(csv_files[1], DataFrame)
    
    former_countries = ["USSR", "CSK", "YUG"]
    df_former = filter(row -> row.country in former_countries, df)
    
    if nrow(df_former) > 0
        df_clean = select(df_former,
            :country => :country_name,
            :year => :year,
            :pop => :population,
            :rgdpch => :gdp_per_capita,
            :p => :price_level
        )
        
        df_clean = filter(row -> 1970 <= row.year <= 1991, df_clean)
        
        output_path = "temp/pwt56_former_clean.csv"
        CSV.write(output_path, df_clean)
        @info "Cleaned PWT 5.6 data for former countries saved to $output_path"
        
        return df_clean
    else
        @warn "No data found for former USSR, Czechoslovakia, or Yugoslavia in PWT 5.6"
        return nothing
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    @info "Cleaning Penn World Table data..."
    clean_pwt71()
    clean_pwt56()
    @info "PWT cleaning complete"
end