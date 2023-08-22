\name{gaussiankernel}
\alias{gaussiankernel}
\title{Compute the Gaussian kernel matrix}
\description{This function provides the Gaussian kernel matrix computed with the median heuristic bandwidth.}
\usage{
gaussiankernel(X)
}
\arguments{
  \item{X}{The samples in the sequence.}
}
\value{
  Returns a numeric matrix, the Gaussian kernel matrix computed with the specified bandwidth.
}
\seealso{
  \code{\link{kerSeg-package}},\code{\link{kerseg1}},\code{\link{kerseg2}}
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
}
