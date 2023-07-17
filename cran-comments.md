Dear CRAN

This is a regular update of the dendroTools R package.

There was an error reported by a dendroTools user, which happened in an 
exceptional case when no significant correlation was found between the proxy and
the climate variables. The error is now fixed.

In previous versions, we agreed to disable vignette checking on some platforms. I believe it makes sense to keep it that way.  

Best,   
Jernej


##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 4.1.1

* rhub Windows Server 2022 (https://builder.r-hub.io/status/original/dendroTools_1.2.10.tar.gz-1e0ff2d94d724acfb85ad3f86bb565f0)
* rhub Fedora Linux (https://builder.r-hub.io/status/original/dendroTools_1.2.10.tar.gz-a6b8c05ef2f546eca99db30eb5111273)
* rhub Ubuntu Linux (https://builder.r-hub.io/status/original/dendroTools_1.2.10.tar.gz-5d8dccb871b940ac8ce38917358ad41b)


* win-check oldrelease (https://win-builder.r-project.org/1PBMTsB7JR2X/00check.log)
* win-check release (https://win-builder.r-project.org/e593Nl4sZeQv/00check.log)
* win-check devel (https://win-builder.r-project.org/1AQYGedjc15l/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
