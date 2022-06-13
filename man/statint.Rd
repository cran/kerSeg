\name{statint}
\alias{statint}
\title{
Compute the test statistics, D and W, for the changed-interval alternatives.
}
\description{
This function provides the test statistics, \eqn{\textrm{D}(t_{1},t_{2})}, \eqn{\textrm{W}(t_{1},t_{2})}, and the weighted \eqn{\textrm{W}(t_{1},t_{2})} for the changed-interval alternatives.
}
\usage{
statint(K, Rtemp, R0, r1, r2)
}
\arguments{
  \item{K}{A kernel matrix of observations in the sequence.}
  \item{Rtemp}{A numeric vector of \eqn{k_{i.}}, the sum of kernel values for each row i.}
  \item{R0}{The term \eqn{R_{0}}, defined in the paper.}
  \item{r1}{The constant in the test statistics \eqn{\textrm{Z}_{W,r1}(t_{1},t_{2})}.}
  \item{r2}{The constant in the test statistics \eqn{\textrm{Z}_{W,r2}(t_{1},t_{2})}.}
}
\value{
  Returns a list of test statistics, \eqn{\textrm{D}(t_{1},t_{2})}, \eqn{\textrm{W}(t_{1},t_{2})}, \eqn{\textrm{W}_{r1}(t_{1},t_{2})}, and \eqn{\textrm{W}_{r2}(t_{1},t_{2})}.
}
\examples{
## Sequence : change in the mean in the middle of the sequence.
d = 50
mu = 2
tau = 50
n = 100
set.seed(1)
y = rbind(matrix(rnorm(d*tau),tau), matrix(rnorm(d*(n-tau),mu/sqrt(d)), n-tau))
K = gaussiankernel(y) # Gaussian kernel matrix
R_temp = rowSums(K)
R0 = sum(K)
a = statint(K, R_temp, R0, r1=1.2, r2=0.8)
}
