## Resubmission
This is a resubmission. CRAN devtools results returned the following error: Package required and available but unsuitable version: 'stats'. To remove the error:

* I changed the Depends field in the Description file: Depends: R (>= 3.4). I decided to do that since the error was produced in an earlier version of R (i.e. R 3.3.3 (2017-03-06)). (Solution was suggested by [R-pkg-devel] mailing list.) 
* I removed the specific version of 'stats' in the Imports field in the Description file. 
* I changed the import in the NAMESCPACE for 'stats'. Now: importFrom("stats", "cor", "lm", "qt", "quantile")

I have also changed the version of dendroExtra to 0.0.2

## Test environments
* local OS X install, R 3.4.0
* Ubuntu precise (12.04.5 LTS) (on travis-ci), R version 3.4.1 (2017-06-30)
* win-builder (R Under development (unstable) (2017-08-03 r73028))

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs


## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
