Dear CRAN

This is a regular update of the dendroTools R package. Approximately three days ago, the CRAN team requested the removal of certain NOTES that appeared during CRAN checks. Specifically, this NOTE pertained to an undocumented argument in one of the internal functions. I have addressed this issue, ensuring that the argument is now properly documented. Alongside this major fix, the update includes several minor enhancements to improve the overall functionality and user experience of the package.

In previous versions, we agreed to disable vignette checking on some platforms. I believe it makes sense to keep it that way.  

Best,   
Jernej


##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 4.1.1

* rhub Windows Server 2022 (https://builder.r-hub.io/status/original/dendroTools_1.2.11.tar.gz-3dbf92ce144b470883613e0e3705024a0)
* win-check oldrelease (https://win-builder.r-project.org/YluCbl0G5Cwd/00check.log)
* win-check release (https://win-builder.r-project.org/C9x7tvD86IST/00check.log)
* win-check devel (https://win-builder.r-project.org/GybK3Wyed753/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
