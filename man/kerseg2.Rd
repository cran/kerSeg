\name{kerseg2}
\alias{kerseg2}
\title{Kernel-based change-point detection for changed-interval alternatives}
\description{This function finds an interval in the sequence where their underlying distribution differs from the rest of the sequence.}
\usage{
kerseg2(n, K, r1=1.2, r2=0.8, l0=0.05*n, l1=0.95*n,
   pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100)
}
\arguments{
  \item{n}{The number of observations in the sequence.}
  \item{K}{The kernel matrix of observations in the sequence.}
  \item{r1}{The constant in the test statistics \eqn{\textrm{Z}_{W,r1}(t_{1},t_{2})}.}
  \item{r2}{The constant in the test statistics \eqn{\textrm{Z}_{W,r2}(t_{1},t_{2})}.}
  \item{l0}{The minimum length of the interval to be considered as a changed interval.}
  \item{l1}{The maximum length of the interval to be considered as a changed interval.}
  \item{pval.appr}{If it is TRUE, the function outputs the p-value approximation based on asymptotic properties.}
  \item{skew.corr}{This argument is useful only when pval.appr=TRUE. If skew.corr is TRUE, the p-value approximation would incorporate skewness correction.}
  \item{pval.perm}{If it is TRUE, the function outputs the p-value from doing B permutations, where B is another argument that you can specify. Doing permutation could be time consuming, so use this argument with caution as it may take a long time to finish the permutation.}
  \item{B}{This argument is useful only when pval.perm=TRUE. The default value for B is 100.}
}
\value{
 Returns a list \code{stat} containing the each scan statistic, \code{tauhat} containing the estimated changed-interval, \code{appr} containing the approximated p-values of the fast tests when argument `pval.appr' is TRUE, and \code{perm} containing the permutation p-values of the fast tests and GKCP when argument `pval.perm' is TRUE. See below for more details.
  \item{seq}{A matrix of each scan statistic (standardized counts).}
  \item{Zmax}{The test statistics (maximum of the scan statistics).}
  \item{tauhat}{An estimate of the two ends of the changed-interval.}
  \item{fGKCP1_bon}{The p-value of \eqn{\textrm{fGKCP}_{1}} obtained by the Bonferroni procedure.}
  \item{fGKCP1_sim}{The p-value of \eqn{\textrm{fGKCP}_{1}} obtained by the Simes procedure.}
  \item{fGKCP2_bon}{The p-value of \eqn{\textrm{fGKCP}_{2}} obtained by the Bonferroni procedure.}
  \item{fGKCP2_sim}{The p-value of \eqn{\textrm{fGKCP}_{2}} obtained by the Simes procedure.}
  \item{GKCP}{The p-value of GKCP obtained by the random permutation.}
}
\seealso{
  \code{\link{kerSeg}}, \code{\link{kerseg2}}, \code{\link{gaussiankernel}}, \code{\link{kerseg1}}
}
\examples{
## Sequence 3: change in both the mean and variance happens on an interval.
d = 50
mu = 2
sigma = 0.5
tau1 = 25
tau2 = 35
n = 50
set.seed(1)
y1 = matrix(rnorm(d*tau1),tau1)
y2 = matrix(rnorm(d*(tau2-tau1),mu/sqrt(d),sigma), tau2-tau1)
y3 = matrix(rnorm(d*(n-tau2)), n-tau2)
y = rbind(y1, y2, y3)
K = gaussiankernel(y)
a = kerseg2(n, K, pval.perm=TRUE, B=1000)
}
