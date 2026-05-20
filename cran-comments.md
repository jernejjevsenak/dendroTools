Dear CRAN Team,

This is a routine update of the dendroTools R package, version 1.2.16.

This release extends the daily and monthly response-functions by allowing analyses of multiple previous years. 
It also includes corresponding updates to summary and plotting methods, so that results from multi-year climate windows are reported and visualized correctly.

In addition, this version improves the robustness of the brnn method by handling occasional unstable model fits more safely, returning missing values for failed windows rather than interrupting the full analysis.

Best,
Jernej

##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 4.1.1

* win-check oldrelease (https://win-builder.r-project.org/tUSHg4H99Sfz/00check.log)
* win-check release (https://win-builder.r-project.org/91t82wKqDtxh/00check.log)
* win-check devel (https://win-builder.r-project.org/fUk7rBr317B2/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
