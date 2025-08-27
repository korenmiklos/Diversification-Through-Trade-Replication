# Data Pipeline

This directory contains the data pipeline for reproducing the Diversification Through Trade analysis from original data sources.

## Structure

- `download/` - Scripts to download raw data from original sources
- `clean/` - Scripts to clean and standardize raw data
- `transform/` - Scripts to transform data into final analysis format

## Data Sources

### Required API Access

1. **UN Comtrade** (Required)
   - Register at https://comtradeplus.un.org/
   - Add API key to `.env` file

2. **World Bank** (Optional, recommended)
   - Register at https://datahelpdesk.worldbank.org/
   - Provides faster access to WDI data

3. **IMF Data** (Optional)
   - Register at https://data.imf.org/
   - Used for exchange rates

### Data Sources Status

| Source | Original URL | Current Status | Action Required |
|--------|--------------|----------------|-----------------|
| EU KLEMS 2008 | euklems.net | ❌ Moved | Use Luiss Lab archive |
| UN National Accounts | data.un.org | ✅ Available | None |
| Penn World Table 7.1 | rug.nl/ggdc | ✅ Available | None |
| WDI October 2015 | worldbank.org | ⚠️ Large file | May need manual download |
| UN Comtrade | comtrade.un.org | ⚠️ New API | Requires API key |

## Setup

1. Copy `.env.template` to `.env`:
   ```bash
   cp .env.template .env
   ```

2. Add your API keys to `.env`

3. Install Julia packages:
   ```bash
   make install
   ```

## Running the Pipeline

### Full Pipeline
```bash
make pipeline
```

### Individual Steps
```bash
make download   # Download raw data
make clean      # Clean raw data
make transform  # Transform to analysis format
```

### Download Specific Sources
```bash
make download-pwt       # Penn World Tables
make download-unna      # UN National Accounts
make download-wdi       # World Development Indicators
make download-comtrade  # UN Comtrade (requires API key)
```

## Output Files

Final cleaned data files are saved to:
- `output/gross_output.csv` - Sectoral gross output by country/year
- `output/value_added.csv` - Sectoral value added by country/year
- `output/trade_flows.csv` - Bilateral trade flows by sector
- `output/price_indices.csv` - Price indices by country/year

## Notes

- Downloads are cached in `temp/download_cache/` to avoid repeated downloads
- Progress is logged to `temp/download_log.csv`
- UN Comtrade downloads implement rate limiting (1.5 second delay between requests)
- Large datasets are downloaded in batches with checkpoint support

## Troubleshooting

### API Rate Limits
If you encounter rate limit errors:
1. Check `temp/download_log.csv` for failed downloads
2. Wait before retrying (usually 1 hour for most APIs)
3. Resume downloads will skip already completed batches

### Missing Data Sources
If original sources are unavailable:
1. Check AGENTS.md for alternative sources
2. Use cached data from the replication package as fallback
3. Document any deviations in your analysis

### Memory Issues
For large datasets:
1. Process data in chunks
2. Use `temp/` directory for intermediate files
3. Clear cache periodically with `rm -rf temp/download_cache/`