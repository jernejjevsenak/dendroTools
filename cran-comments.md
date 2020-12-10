Dear CRAN, 

I am resubmitting dendroTools R package with corrections related to testthat. There are no other changes.

I have tested new version on all platforms and all examples were finished successfully. 

Please note extensive vignette examples. In previous versions we agreed to turn off vignette checks on some platforms. I believe it makes sense to keep it so. Thank you very much.

Best,
Jernej 


##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 3.6.2
* rhub platforms:  	
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit (https://builder.r-hub.io/status/original/dendroTools_1.1.1.tar.gz-5876aaab93524986a528c2d245b938bd)
  - Ubuntu Linux 16.04 LTS, R-release, GCC (https://builder.r-hub.io/status/original/dendroTools_1.1.1.tar.gz-e8ae635759454daa858dec5d07209fe7)
  - Fedora Linux, R-devel, clang, gfortran (https://builder.r-hub.io/status/original/dendroTools_1.1.1.tar.gz-e5c50c64648241f5814a0a63e54697c5)
* win-check oldrelease (https://win-builder.r-project.org/Fhz8ZbccmYwk/00check.log)
* win-check release (https://win-builder.r-project.org/9mgN44lls6ZT/00check.log)
* win-check devel (https://win-builder.r-project.org/o0I8S8PLF5eC/00check.log)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
