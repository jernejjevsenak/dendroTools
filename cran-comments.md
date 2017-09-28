##  Resubmission

Improvements: 
- CRAN maintainer asked about the reference to be added to the DESCRIPTION field. Currently, there is no reference about our method. But the article for scientific journal is already being prepared. As soon as the article will be published, the reference will be added to the next version of our package. 

- CRAN maintainer asked to add more executable examples in our Rd-files. New examples were added to the calculate_measures(), compare_methods(), daily_response() and years_to_rownames(). However, most of those examples were also wrapped into dontrun{}. The reason for that is: even the most simple ones exceed the RCMD maximum elapsed time for function examples. 
Internal functions, i.e. count_ones() and iter() do not have examples, because we get an error from RCMD check if examples are added. This is logical since those functions are not exported and examples can not be checked. 

Please note: Previously, package was called dendroExtra (already on CRAN), but authors decided to publish the new version under new name dendroTools. There are several reasons for that:  1) there is already a package called dendroextras, which is very close to dendroExtra and 2) the future development of the package will be focused on the implementation of new tools, so the new name is definately more appropriate. We kindly ask CRAN for understanding and to approve the change of the name of this package. Here, we ask CRAN, when dendroTools will be on CRAN, to remove the dendroExtra package from CRAN. 

## Test environments
* local OS X install, R 3.4.0
* Ubuntu precise (Ubuntu 14.04.5 LTS) (on travis-ci), R version 3.4.1 (2017-01-27)
* win-builder R Under development (unstable) (2017-09-12 r73242)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs

## Downstream dependencies
We have also run R CMD check on downstream dependencies of dendroTools
https://github.com/jernejjevsenak/dendroTools/blob/master/revdep/checks.rds

All packages that we could install passed. 
