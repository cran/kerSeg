\name{kerSeg}
\alias{kerSeg-package}
\title{
  New kernel-based change-point detection
}
\description{
  This package can be used to detect change-points where the distributions abruptly change. The Gaussian kernel with the median heuristic, which is the median of all pairwise distances among observations, is used. 
}
\details{
  To compute the Gaussian kernel matrix with the median heuristic bandwidth, the function \code{\link{gaussiankernel}} should be used. The main functions are \code{\link{kerseg1}} for the single change-point alternative and \code{\link{kerseg2}} for the changed-interval alternative.
}
\author{
Hoseung Song and Hao Chen

Maintainer: Hoseung Song (hosong@ucdavis.edu)
}
\references{
  Song, H. and Chen, H. (2022). New kernel-based change-point detection. arXiv:2206.01853
}
\seealso{
  \code{\link{kerseg1}}, \code{\link{kerseg2}}, \code{\link{gaussiankernel}}
}
\examples{
## Sequence 1: change in the mean in the middle of the sequence.
d = 50
mu = 2
tau = 15
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
