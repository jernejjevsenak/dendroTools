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
  - Windows Server 2008 R2 SP1, R-devel, 32/64 bit (https://builder.r-hub.io/status/original/dendroTools_1.1.0.tar.gz-c3753a4ad64248e3a2f2f77f968e7104)
  - Ubuntu Linux 16.04 LTS, R-release, GCC (https://builder.r-hub.io/status/original/dendroTools_1.1.0.tar.gz-e1e430758d1e4d1f8192e7baf02802e0)
  - Fedora Linux, R-devel, clang, gfortran (https://builder.r-hub.io/status/original/dendroTools_1.1.0.tar.gz-b12498b92aed45838cb685270fdd5834)
* win-check oldrelease (https://win-builder.r-project.org/78aJ4Sy402lW)
* win-check release (https://win-builder.r-project.org/RQnYife9AjO0)
* win-check devel (https://win-builder.r-project.org/Sv2zV7bQdY8I)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
