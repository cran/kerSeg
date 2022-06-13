\name{skew}
\alias{skew}
\title{
Compute some components utilized in the third moment fomulas.
}
\description{
This function provides some components used in the third moment fomulas.
}
\usage{
skew(K, Rtemp, Rtemp2, R0, R2)
}
\arguments{
  \item{K}{A kernel matrix of observations in the sequence.}
  \item{Rtemp}{A numeric vector of \eqn{k_{i.}}, the sum of kernel values for each row i.}
  \item{Rtemp2}{A numeric vector, the sum of squared kernel values for each row i.}
  \item{R0}{The term \eqn{R_{0}}, defined in the paper.}
  \item{R2}{The term \eqn{R_{2}}, defined in the paper.}
}
\value{
  Returns a list of components used in the third moment fomulas.
}