Dear CRAN
After updating my dendroTools R package yesterday, I have received an email from CRAN with the following content: 

“Dear maintainer,

Dear CRAN
After updating my dendroTools R package on March 8, 2021, I received an email from CRAN with the following content: Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_dendroTools.html>.

Please correct before 2021-03-22 to safely retain your package on CRAN.

It seems we need to remind you of the CRAN policy:

'Packages which use Internet resources should fail gracefully with an informative message
if the resource is not available or has changed (and not give a check warning nor error).'

This needs correction whether or not the resource recovers.

The CRAN Team”

I can see that the problem was related to one of my examples in a vignette. I decided to remove that example entirely, so I assume the problem is now fixed. I've also seen similar issues reported previously for other packages where the examples used data from the internet. However, all of my examples only rely on data within the dendroTools r package. I tested the updated package on all platforms and got 0 ERRORS, 0 WARNING and 1 NOTE: Days since last update: 6.

Since I was asked by CRAN to update my package as soon as possible, I believe NOTE is acceptable this time.

Please also note the extensive vignette examples. In previous versions, we agreed to disable vignette checking on some platforms. I believe it makes sense to keep it that way.

Many thanks at this point. 
Jernej


##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 4.0.4
* rhub platforms (https://builder.r-hub.io/status/original/dendroTools_1.1.3.tar.gz-d3de58f10edd4d5f9e10c2b102f6c95b)
* win-check oldrelease (https://win-builder.r-project.org/fnx831RO88br/00check.log)
* win-check release (https://win-builder.r-project.org/918gc45WjbAp/00check.log)
* win-check devel (https://win-builder.r-project.org/Elf83DXaqNU3/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 1 NOTE (please see the explanation above)

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
