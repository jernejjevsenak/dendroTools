Dear CRAN
After updating my dendroTools R package yesterday, I have received an email from CRAN with the following content: 

“Dear maintainer,

Please see the problems shown on
<https://cran.r-project.org/web/checks/check_results_dendroTools.html>.

Please correct before 2021-03-22 to safely retain your package on CRAN.

It seems we need to remind you of the CRAN policy:

'Packages which use Internet resources should fail gracefully with an informative message
if the resource is not available or has changed (and not give a check warning nor error).'

This needs correction whether or not the resource recovers.

The CRAN Team”

As I can see, similar problems usually occur when the examples/vignettes in the package are really based on external data. However, all my examples and vignettes only rely on data available within my package, so I'm not sure what caused the problem yesterday. However, I have removed the problematic example from my vignette, and I believe I have solved the problem. I have run all the examples and checked the necessary tests and received no ERROR, WARNING or NOTE.

Please note extensive vignette examples. In previous versions we agreed to turn off vignette checks on some platforms. I believe it makes sense to keep it so. Thank you very much.

Best,
Jernej 


##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 3.6.2
* rhub platforms (https://builder.r-hub.io/status/original/dendroTools_1.1.2.tar.gz-dcb5bb950c1047448b1ebf8cb8d73723)
* win-check oldrelease (https://win-builder.r-project.org/OI5FbJHWJYg8/00check.log)
* win-check release (https://win-builder.r-project.org/08U1w71HpKmI/00check.log)
* win-check devel (https://win-builder.r-project.org/o0I8S8PLF5eC/00check.log)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
