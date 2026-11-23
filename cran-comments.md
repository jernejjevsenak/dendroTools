Dear CRAN Team,

This is a major update of the dendroTools R package, version 2.0.0. The current CRAN version is 1.2.16, so the changes listed in NEWS.md under 1.2.17 are also new on CRAN.

The four response functions, daily_response(), monthly_response(), daily_response_seascorr() and monthly_response_seascorr(), can now detrend the aggregated climate data with the methods of dplR::detrend() and calculate the Standardized Precipitation-Evapotranspiration Index (SPEI) of every window. The new function water_balance() calculates the climatic water balance for the SPEI. The seascorr functions can also return the partial correlations of the control variable and partition the explained variance between both climate variables. New significance criteria account for the number of windows screened, and the bootstrap of correlations is much faster.

The new major version reflects changes that affect existing code and results: the argument remove_insignificant_boot was replaced by remove_insignificant_method, the returned element $boot_ci_filter was renamed $significance_filter, and with boot = TRUE the functions now report the coefficients of the observed data instead of the mean of the bootstrap replicates.

The release also adds example data from a black pine site in Slovenia and three new vignettes.

Best,
Jernej

## Test environments
* local Windows install, R 4.4.1

* win-check oldrelease (https://win-builder.r-project.org/jG8gmJ8v7Lmq/00check.log)
* win-check release (https://win-builder.r-project.org/Ysg9Bz9z9pCR/00check.log)
* win-check devel (https://win-builder.r-project.org/hNKBVJ23fc23/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 1 NOTE.

The NOTE lists possibly misspelled words in DESCRIPTION: Evapotranspiration, Jevsenak, Levanic, Meko, SPEI, al and et. They are spelled correctly: Jevsenak, Levanic and Meko are author names, "et al." is part of the references, and Evapotranspiration and SPEI come from the name of the Standardized Precipitation-Evapotranspiration Index (SPEI).

## Downstream dependencies
There are currently no downstream dependencies for this package.
