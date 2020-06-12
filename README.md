# triadic-bodily-synchrony

Raw body-movement data sample and analysis scripts for submission:

-- Dale, R., Bryant, G. A., Manson, J. H., & Gervais, M. M. (resubmitted). Body synchrony in triadic interaction. 

`MAIN_SCRIPT.R` begins the overall analysis, and calls upon a number of scripts, including (1) `functions.R`, containing some key functions, (2) `bodyExtract.R` which stores flow data, (3) `computeCorrelations.R` which calculates cross-correlations, and (4) `combineMansonData.R` which integrates the individual differences measures. (5) `plotCCF.R` generates some detailed plots from the manuscript of the cross-correlation function.

NB: ZIP files of raw data must be unpacked in `raw_flow_data` folder.
