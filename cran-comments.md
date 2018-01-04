##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 3.4.0
* Ubuntu precise (Ubuntu 14.04.5 LTS) (on travis-ci), R version 3.4.2 (2017-01-27)
* win-builder R Under development (unstable) (2017-09-12 r73242)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 

## Vignete 

I created vignette for the dendroTools. However, due to the package functionality, it takes around 30 minutes to build vignette. The core function daily_response() runs through daily environmental data and calculates averages of different ranges and afterwards compares this average with some dependent variable and returns some linear / nonlinear metric. By default, it is computationally intensive, and it is meant to be so.
