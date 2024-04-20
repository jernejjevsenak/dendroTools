Dear CRAN

I am submitting an updated dendroTools R package after receiving an email from Prof. Brian Ripley specifically pointing out the problems associated with progress bars in non-interactive use. I wrapped all progress bars in interacitve(), which should now fix the problem.

Since »Additional issues « appeared only after my last update, I assume they are related to the additional features that were implemented after the last update. Therefore, I have removed all features and dependencies that were introduced in my last update.

Some Additional issues were related to testthat(), which has now been reduced.
Some Additional issues were also related to the extensive examples, which have now been reduced to a minimum.

In addition:
- some less important features were removed from the package to make it more efficient
- based on the user suggestions, new detrending method was introduced, i.e. Simple Linear Detrending (SLD)
- also two new features were implemented that allow for skipping specific windows and thus speed up the calculation of correlation coefficients

After the modifications, the dendroTools R package has been successfully tested on the local machine and on all standard platforms, rhub, win_devel, win_oldrelease, wind_release.

I hope that these updates fulfill the criteria to keep the dendroTools R package on CRAN. 

If there are any further requests, I will be happy to take action immediately.

In previous versions, we agreed to disable vignette checking on some platforms. I believe it makes sense to keep it that way.  

Best,   
Jernej

##  Resubmission
* This is a resubmission of the package dendroTools.

## Test environments
* local OS X install, R 4.1.1

* win-check oldrelease (https://win-builder.r-project.org/w3602oi8zKg8/00check.log)
* win-check release (https://win-builder.r-project.org/7q8SV3vgTPXR/00check.log)
* win-check devel (https://win-builder.r-project.org/PfaxzKj6BnV1/00check.log)

## R CMD check results
There were 0 ERRORs, 0 WARNINGs and 0 NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
