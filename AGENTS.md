# Data Source Investigation Report
Date: 2025-08-27

## Executive Summary

Investigation of data sources for the Diversification Through Trade replication pipeline revealed significant changes to original 2008-2012 data locations. Most sources have migrated to new platforms with updated access methods.

## Data Source Status

### 1. EU KLEMS (March 2008 Release)
- **Original URL:** `http://www.euklems.net/euk08i.shtml` 
- **Status:** ❌ Not accessible (404)
- **Current Location:** Luiss Lab (https://euklems-intanprod-llee.luiss.it/)
- **Alternative Archive:** University of Groningen with DOI: 10.34894/MGSB4H
- **Data Format:** Excel, CSV, R, STATA
- **Access Method:** Direct download after registration
- **Action Required:** Update download URLs and potentially registration process

### 2. UN National Accounts (SNA Table 2.1)
- **Original URL:** `http://data.un.org/Data.aspx?d=SNA&f=group_code:201`
- **Status:** ✅ Still accessible
- **Data Format:** XML, CSV export available
- **Coverage:** 224 countries, 1946-2023
- **Limitation:** 100,000 records max via web interface
- **Access Method:** Direct web download or bulk API
- **Action Required:** Implement pagination for large queries

### 3. Penn World Table 7.1
- **Original URL:** `https://www.rug.nl/ggdc/productivity/pwt/pwt-releases/pwt-7.1`
- **Status:** ✅ Accessible as archive
- **Data Format:** Zipped Excel files
- **Coverage:** 189 countries, 1950-2010
- **Note:** Current version is PWT 10.01, but 7.1 remains available for replication
- **Access Method:** Direct download
- **Action Required:** None, original source intact

### 4. World Development Indicators (October 2015)
- **Original URL:** `http://databank.worldbank.org/data/download/archive/WDI_excel_2015_10.zip`
- **Status:** ⚠️ File too large to verify directly
- **Current Platform:** World Bank Data Catalog
- **License:** CC BY 4.0
- **Access Method:** API available with registration
- **Action Required:** Obtain API key for programmatic access

### 5. UN Comtrade
- **Original URL:** `https://comtrade.un.org/`
- **Status:** ⚠️ Redirected to new system
- **Current Platform:** https://comtradeplus.un.org/
- **API Changes:** Now subscription-based with free tier limitations
- **Authentication:** API key required
- **Rate Limits:** Free tier has strict limits
- **Action Required:** Register for API key, implement rate limiting

### 6. UNIDO INDSTAT2
- **Status:** 🔒 Paid subscription required
- **Access:** Manual credentials/API with institutional subscription
- **Alternative:** Use replication package's cached data if no access
- **Action Required:** Check for institutional access or use fallback

### 7. IMF Exchange Rates
- **Platform:** IMF Data API
- **Authentication:** API key required
- **Alternative:** OECD or World Bank official FX rates
- **Action Required:** Register for IMF API access

## Implementation Recommendations

### Priority 1: Data Access Setup
1. Register for required API keys:
   - UN Comtrade Plus
   - World Bank Data API
   - IMF Data API
2. Document credentials in `.env` file (not committed)
3. Implement fallback to cached data where original sources unavailable

### Priority 2: Code Architecture
1. Create modular download scripts with retry logic
2. Implement caching to avoid repeated downloads
3. Use bead for version management of downloaded data
4. Add validation checksums for downloaded files

### Priority 3: Rate Limiting Strategy
- UN Comtrade: Implement exponential backoff with 1-second minimum delay
- World Bank: Batch requests to stay within limits
- Add progress indicators for long-running downloads

## Technical Decisions

### Download Strategy
- Use Julia's `Downloads.jl` and `HTTP.jl` for API calls
- Implement parallel downloads where possible
- Cache responses in `temp/` to enable resumption
- Use bead to track data lineage

### Data Format Handling
- Standardize on CSV for intermediate storage
- Use DataFrames.jl for processing
- Convert Excel files immediately to CSV
- Preserve original downloaded files in `input/raw/`

### Error Handling
- Implement retry logic with exponential backoff
- Log all download attempts and failures
- Provide clear error messages with resolution steps
- Fall back to cached data with warnings

## Next Steps

1. Create `code/create/download/` directory structure
2. Implement base download utilities
3. Write source-specific download scripts
4. Update Makefile with download targets
5. Test with small data samples
6. Document API registration process in README

## Notes

- EU KLEMS migration most significant change requiring URL updates
- UN Comtrade API changes may impact download speed due to rate limits
- Consider implementing parallel processing for independent data sources
- May need institutional access for UNIDO data