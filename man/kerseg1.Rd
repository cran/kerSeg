\name{kerseg1}
\alias{kerseg1}
\title{Kernel-based change-point detection for single change-point alternatives}
\description{This function finds a break point in the sequence where the underlying distribution changes.}
\usage{
kerseg1(n, K, r1=1.2, r2=0.8, n0=0.05*n, n1=0.95*n,
   pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100)
}
\arguments{
  \item{n}{The number of observations in the sequence.}
  \item{K}{The kernel matrix of observations in the sequence.}
  \item{r1}{The constant in the test statistics \eqn{\textrm{Z}_{W,r1}(t)}.}
  \item{r2}{The constant in the test statistics \eqn{\textrm{Z}_{W,r2}(t)}.}
  \item{n0}{The starting index to be considered as a candidate for the change-point.}
  \item{n1}{The ending index to be considered as a candidate for the change-point.}
  \item{pval.appr}{If it is TRUE, the function outputs the p-value approximation based on asymptotic properties.}
  \item{skew.corr}{This argument is useful only when pval.appr=TRUE. If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
  \item{pval.perm}{If it is TRUE, the function outputs the p-value from doing B permutations, where B is another argument that you can specify. Doing permutation could be time consuming, so use this argument with caution as it may take a long time to finish the permutation.}
  \item{B}{This argument is useful only when pval.perm=TRUE. The default value for B is 100.}
}
\value{
  Returns a list \code{stat} containing the each scan statistic, \code{tauhat} containing the estimated location of change-point, \code{appr} containing the approximated p-values of the fast tests when argument `pval.appr' is TRUE, and \code{perm} containing the permutation p-values of the fast tests and GKCP when argument `pval.perm' is TRUE. See below for more details.
  \item{seq}{A vector of each scan statistic (standardized counts).}
  \item{Zmax}{The test statistics (maximum of the scan statistics).}
  \item{tauhat}{An estimate of the location of the change-point.}
  \item{fGKCP1_bon}{The p-value of \eqn{\textrm{fGKCP}_{1}} obtained by the Bonferroni procedure.}
  \item{fGKCP1_sim}{The p-value of \eqn{\textrm{fGKCP}_{1}} obtained by the Simes procedure.}
  \item{fGKCP2_bon}{The p-value of \eqn{\textrm{fGKCP}_{2}} obtained by the Bonferroni procedure.}
  \item{fGKCP2_sim}{The p-value of \eqn{\textrm{fGKCP}_{2}} obtained by the Simes procedure.}
  \item{GKCP}{The p-value of GKCP obtained by the random permutation.}
}
\seealso{
  \code{\link{kerSeg-package}}, \code{\link{kerseg1}}, \code{\link{gaussiankernel}}, \code{\link{kerseg2}}
}
\examples{
## Sequence 1: change in the mean in the middle of the sequence.
d = 50
mu = 2
tau = 25
n = 50
set.seed(1)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d)), n-tau))
K = gaussiankernel(y) # Gaussian kernel matrix
a = kerseg1(n, K, pval.perm=TRUE, B=1000)
# output results based on the permutation and the asymptotic results.
# the scan statistics can be found in a$scanZ.
# the approximated p-values can be found in a$appr.
# the permutation p-values can be found in a$perm.

## Sequence 2: change in both the mean and variance away from the middle of the sequence.
d = 50
mu = 2
sigma = 0.7
tau = 35
n = 50
set.seed(1)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d),sigma), n-tau))
K = gaussiankernel(y)
a = kerseg1(n, K, pval.perm=TRUE, B=1000)
}
