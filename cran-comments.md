Dear CRAN, 

I am resubmitting my package dendroTools. Previous version (1.0.8) was rejected by CRAN due to wrong DOI of my publication. I have now corrected wrongly specified DOI, and I also do not get any NOTE from rhub or any other test environment.

I have tested new version on all platforms and all examples were finished successfully. 

Please note extensive vignette examples. In previous versions we agreed to turn off vignette checks on some platforms. I believe it makes sense to keep it so. Thank you very much.

Best,
Jernej 


##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 3.6.2
* rhub platforms:  	
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  - Ubuntu Linux 16.04 LTS, R-release, GCC
  - Fedora Linux, R-devel, clang, gfortran
* win-check oldrelease (https://win-builder.r-project.org/oV3fKgW1JovZ/00check.log)
* win-check release (https://win-builder.r-project.org/BEH01siE6OLd/00check.log)
* win-check devel (https://win-builder.r-project.org/5kfyj33RkWkG)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
