## R CMD check results
0 errors v | 0 warnings v | 0 notes v

R CMD check succeeded

## Rhub

All platforms worked except for Debian, where I received the following errors:

ERROR: dependency ‘Rcpp’ is not available for package ‘RcppArmadillo’
* removing ‘/home/docker/R/RcppArmadillo’

ERROR: dependencies ‘Rcpp’, ‘RcppArmadillo’ are not available for package ‘bmabasket’
* removing ‘/home/docker/R/bmabasket’

This seems to be an Rhub-specific issue, and not related to the code of the
package.


## Github Actions
I have also run R CMD CHECK via Github Actions:
https://github.com/ethan-alt/bmabasket/actions
