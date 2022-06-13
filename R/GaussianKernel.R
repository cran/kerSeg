# function to get the Gaussian kernel
gaussiankernel = function(X) {
  n = nrow(X)
  Dx = dist(X)^2
  sigma = median(Dx)
  if (sigma == 0) {
    sigma = sigma + 0.1
  }
  K = exp(-as.matrix(Dx)/sigma/2) - diag(n) # kernel matrix
  return(K)
}