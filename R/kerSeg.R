# n: Total length of the sequence
# K: Kernel matrix
kerseg1 = function(n, K, r1=1.2, r2=0.8, n0=0.05*n, n1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100) {
  R = list()

  n0 = ceiling(n0)
  n1 = floor(n1)

  n0_us = n0
  n1_us = n1

  if(n0_us<=1){
    message("  Note: Starting index has been set to n0 = 2 as the test statistic is not well-defined for t<2. \n")
  }
  if(n1_us>=n-1){
    message("  Note: Ending index has been set to n1 =", n-2, " as the test statistic is not well-defined for t>",n-2,". \n")
  }

  if(n0<2){
    n0=2
  }
  if(n1>(n-2)){
    n1=n-2
  }

  scanZ = computeSTAT1(n, K, r1, r2, n0, n1)

  R$stat$GPK$seq = scanZ$GKCP$seq
  R$stat$GPK$Zmax = scanZ$GKCP$Zmax
  R$stat$ZD$seq = scanZ$Z_D$seq
  R$stat$ZD$Zmax = scanZ$Z_D$Zmax
  R$stat$ZW$seq = scanZ$Z_W$seq
  R$stat$ZW$Zmax = scanZ$Z_W$Zmax
  R$stat$ZW1$seq = scanZ$Z_W1$seq
  R$stat$ZW1$Zmax = scanZ$Z_W1$Zmax
  R$stat$ZW2$seq = scanZ$Z_W2$seq
  R$stat$ZW2$Zmax = scanZ$Z_W2$Zmax
  R$tauhat = scanZ$GKCP$tauhat

  if (pval.appr==TRUE) {
    mypval1 = pval1(n, K, scanZ, r1, r2, skew.corr, n0, n1)
  }
  if (pval.perm==TRUE){
    mypval2 = permpval1(n, K, scanZ, r1, r2, B, n0, n1)
  }

  if (pval.appr==TRUE) {
    temp_loc = c(mypval1$Z_D, mypval1$Z_W1, mypval1$Z_W2)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 3*c(temp_loc1$x[1], temp_loc1$x[2], temp_loc1$x[3])
    R$appr$fGKCP1_bon = min(temp_loc2)

    temp_loc2 = c(3*temp_loc1$x[1], 1.5*temp_loc1$x[2], temp_loc1$x[3])
    R$appr$fGKCP1_sim= min(temp_loc2)

    temp_loc = c(mypval1$Z_W1, mypval1$Z_W2)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 2*c(temp_loc1$x[1], temp_loc1$x[2])
    R$appr$fGKCP2_bon = min(temp_loc2)

    temp_loc2 = c(2*temp_loc1$x[1], temp_loc1$x[2])
    R$appr$fGKCP2_sim = min(temp_loc2)
  }

  if (pval.perm==TRUE) {
    temp_loc = c(mypval2$Z_D$pval, mypval2$Z_W1$pval, mypval2$Z_W2$pval)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 3*c(temp_loc1$x[1], temp_loc1$x[2], temp_loc1$x[3])
    R$perm$fGKCP1_bon = min(temp_loc2)

    temp_loc2 = c(3*temp_loc1$x[1], 1.5*temp_loc1$x[2], temp_loc1$x[3])
    R$perm$fGKCP1_sim = min(temp_loc2)

    temp_loc = c(mypval2$Z_W1$pval, mypval2$Z_W2$pval)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 2*c(temp_loc1$x[1], temp_loc1$x[2])
    R$perm$fGKCP2_bon = min(temp_loc2)

    temp_loc2 = c(2*temp_loc1$x[1], temp_loc1$x[2])
    R$perm$fGKCP2_sim = min(temp_loc2)

    R$perm$GKCP = mypval2$GKCP$pval
  }

  return(R)
}


# n: Total length of the sequence
# K: Kernel matrix
kerseg2 = function(n, K, r1=1.2, r2=0.8, l0=0.05*n, l1=0.95*n, pval.appr=TRUE, skew.corr=TRUE, pval.perm=FALSE, B=100) {
  R = list()

  l0 = ceiling(l0)
  l1 = floor(l1)

  scanZ = computeSTAT2(n, K, r1, r2, l0, l1)

  R$stat$GPK$seq = scanZ$GKCP$seq
  R$stat$GPK$Zmax = scanZ$GKCP$Zmax
  R$stat$ZD$seq = scanZ$Z_D$seq
  R$stat$ZD$Zmax = scanZ$Z_D$Zmax
  R$stat$ZW$seq = scanZ$Z_W$seq
  R$stat$ZW$Zmax = scanZ$Z_W$Zmax
  R$stat$ZW1$seq = scanZ$Z_W1$seq
  R$stat$ZW1$Zmax = scanZ$Z_W1$Zmax
  R$stat$ZW2$seq = scanZ$Z_W2$seq
  R$stat$ZW2$Zmax = scanZ$Z_W2$Zmax
  R$tauhat = scanZ$GKCP$tauhat

  if (pval.appr==TRUE) {
    mypval1 = pval2(n, K, scanZ, r1, r2, skew.corr, l0, l1)
  }
  if (pval.perm==TRUE){
    mypval2 = permpval2(n, K, scanZ, r1, r2, B, l0, l1)
  }

  if (pval.appr==TRUE) {
    temp_loc = c(mypval1$Z_D, mypval1$Z_W1, mypval1$Z_W2)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 3*c(temp_loc1$x[1], temp_loc1$x[2], temp_loc1$x[3])
    R$appr$fGKCP1_bon = min(temp_loc2)

    temp_loc2 = c(3*temp_loc1$x[1], 1.5*temp_loc1$x[2], temp_loc1$x[3])
    R$appr$fGKCP1_sim = min(temp_loc2)

    temp_loc = c(mypval1$Z_W1, mypval1$Z_W2)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 2*c(temp_loc1$x[1], temp_loc1$x[2])
    R$appr$fGKCP2_bon = min(temp_loc2)

    temp_loc2 = c(2*temp_loc1$x[1], temp_loc1$x[2])
    R$appr$fGKCP2_sim = min(temp_loc2)
  }

  if (pval.perm==TRUE) {
    temp_loc = c(mypval2$Z_D$pval, mypval2$Z_W1$pval, mypval2$Z_W2$pval)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 3*c(temp_loc1$x[1], temp_loc1$x[2], temp_loc1$x[3])
    R$perm$fGKCP1_bon = min(temp_loc2)

    temp_loc2 = c(3*temp_loc1$x[1], 1.5*temp_loc1$x[2], temp_loc1$x[3])
    R$perm$fGKCP1_sim = min(temp_loc2)

    temp_loc = c(mypval2$Z_W1$pval, mypval2$Z_W2$pval)
    temp_loc1 = sort(temp_loc, index.return=TRUE)

    temp_loc2 = 2*c(temp_loc1$x[1], temp_loc1$x[2])
    R$perm$fGKCP2_bon = min(temp_loc2)

    temp_loc2 = c(2*temp_loc1$x[1], temp_loc1$x[2])
    R$perm$fGKCP2_sim = min(temp_loc2)

    R$perm$GKCP$pval = mypval2$GKCP$pval
  }

  return(R)
}


computeSTAT1 = function(n, K, r1, r2, n0=ceiling(0.05*n), n1=floor(0.95*n)) {
  if(n0<2){
    n0=2
  }
  if(n1>(n-2)){
    n1=n-2
  }

  t = 1:n
  temp = n0:n1

  R_temp = rowSums(K)

  Kx = Ky = rep(0, n)
  Kx[2] = sum(K[1:2,1:2])
  Ky[2] = sum(K[3:n,3:n])
  for (i in 3:(n-2)) {
    temp_col = K[,i]
    add = sum(temp_col[1:i])
    subtract = R_temp[i] - add
    Kx[i] = Kx[i-1] + 2*add
    Ky[i] = Ky[i-1] - 2*subtract
  }
  Kx = Kx/t/(t-1)
  Ky = Ky/(n-t)/(n-t-1)

  R0 = sum(K)
  mu_Kx = R0/n/(n-1)
  mu_Ky = R0/n/(n-1)

  R1 = sum(K^2)
  R2 = sum(R_temp^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  p1 = t*(t-1)/n/(n-1)
  p2 = p1*(t-2)/(n-2)
  p3 = p2*(t-3)/(n-3)

  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = q1*(n-t-2)/(n-2)
  q3 = q2*(n-t-3)/(n-3)

  var_Kx = (2*R1*p1 + 4*R2*p2 + R3*p3)/t/t/(t-1)/(t-1) - mu_Kx^2
  var_Ky = (2*R1*q1 + 4*R2*q2 + R3*q3)/(n-t)/(n-t)/(n-t-1)/(n-t-1) - mu_Ky^2
  cov_Kx_Ky = R3/n/(n-1)/(n-2)/(n-3) - mu_Kx*mu_Ky

  scanZ = list()

  # test statistic Z_D
  u.D = t*(t-1)
  v.D = -(n-t)*(n-t-1)
  mean.D = mu_Kx*u.D + mu_Ky*v.D
  var.D = (u.D^2)*var_Kx + (v.D^2)*var_Ky + 2*u.D*v.D*cov_Kx_Ky
  Z.D = (Kx*u.D + Ky*v.D - mean.D)/sqrt(var.D)
  Z.D = abs(Z.D)
  tauhat = temp[which.max(Z.D[n0:n1])]
  ZD = list(tauhat=tauhat, Zmax=Z.D[tauhat], seq=Z.D)
  scanZ$Z_D = ZD

  # test statistic Z_W
  u.W = t*(t-1)*(n-t)/n
  v.W = (n-t)*(n-t-1)*t/n
  mean.W = mu_Kx*u.W + mu_Ky*v.W
  var.W = var_Kx*u.W^2 + var_Ky*v.W^2 + 2*u.W*v.W*cov_Kx_Ky
  Z.W = (Kx*u.W + Ky*v.W - mean.W)/sqrt(var.W)
  tauhat = temp[which.max(Z.W[n0:n1])]
  ZW = list(tauhat=tauhat, Zmax=Z.W[tauhat], seq=Z.W)
  scanZ$Z_W = ZW

  # test statistic Z_W_r1
  u.W1 = r1*t*(t-1)*(n-t)/n
  v.W1 = (n-t)*(n-t-1)*t/n
  mean.W1 = mu_Kx*u.W1 + mu_Ky*v.W1
  var.W1 = var_Kx*u.W1^2 + var_Ky*v.W1^2 + 2*u.W1*v.W1*cov_Kx_Ky
  Z.W1 = (Kx*u.W1 + Ky*v.W1 - mean.W1)/sqrt(var.W1)
  tauhat = temp[which.max(Z.W1[n0:n1])]
  ZW1 = list(tauhat=tauhat, Zmax=Z.W1[tauhat], seq=Z.W1)
  scanZ$Z_W1 = ZW1

  # test statistic Z_W_r2
  u.W2 = r2*t*(t-1)*(n-t)/n
  v.W2 = (n-t)*(n-t-1)*t/n
  mean.W2 = mu_Kx*u.W2 + mu_Ky*v.W2
  var.W2 = var_Kx*u.W2^2 + var_Ky*v.W2^2 + 2*u.W2*v.W2*cov_Kx_Ky
  Z.W2 = (Kx*u.W2 + Ky*v.W2 - mean.W2)/sqrt(var.W2)
  tauhat = temp[which.max(Z.W2[n0:n1])]
  ZW2 = list(tauhat=tauhat, Zmax=Z.W2[tauhat], seq=Z.W2)
  scanZ$Z_W2 = ZW2

  # GKCP
  u.W = t*(t-1)*(n-t)/n
  v.W = (n-t)*(n-t-1)*t/n
  mean.W = mu_Kx*u.W + mu_Ky*v.W
  var.W = var_Kx*u.W^2 + var_Ky*v.W^2 + 2*u.W*v.W*cov_Kx_Ky
  Z.W = (Kx*u.W + Ky*v.W - mean.W)/sqrt(var.W)
  gkcp = Z.D^2 + Z.W^2
  tauhat = temp[which.max(gkcp[n0:n1])]
  GKCP = list(tauhat=tauhat, Zmax=gkcp[tauhat], seq=gkcp)
  scanZ$GKCP = GKCP

  return(scanZ)
}


computeSTAT2 = function(n, K, r1, r2, l0=ceiling(0.05*n), l1=floor(0.95*n)) {
  if(l0<2){
    l0=2
  }
  if(l1>(n-2)){
    l1=n-2
  }

  t = 1:n
  temp = l0:l1

  R_temp = rowSums(K)
  R0 = sum(K)

  a = statint(K, R_temp, R0, r1, r2)
  D = a$D
  W = a$W
  W1 = a$W1
  W2 = a$W2

  dif = matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){
      dif[i,j] = j-i
    }
  }
  difv = as.vector(t(dif))
  ids = which(difv>0)
  ids2 = which((difv>=l0) & (difv<=l1))

  mu_Kx = R0/n/(n-1)
  mu_Ky = R0/n/(n-1)

  R1 = sum(K^2)
  R2 = sum(R_temp^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  p1 = t*(t-1)/n/(n-1)
  p2 = p1*(t-2)/(n-2)
  p3 = p2*(t-3)/(n-3)

  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = q1*(n-t-2)/(n-2)
  q3 = q2*(n-t-3)/(n-3)

  var_Kx = (2*R1*p1 + 4*R2*p2 + R3*p3)/t/t/(t-1)/(t-1) - mu_Kx^2
  var_Ky = (2*R1*q1 + 4*R2*q2 + R3*q3)/(n-t)/(n-t)/(n-t-1)/(n-t-1) - mu_Ky^2
  cov_Kx_Ky = R3/n/(n-1)/(n-2)/(n-3) - mu_Kx*mu_Ky

  scanZ = list()

  # test statistic Z_D
  u.D = t*(t-1)
  v.D = -(n-t)*(n-t-1)
  mean.D = mu_Kx*u.D + mu_Ky*v.D
  var.D = (u.D^2)*var_Kx + (v.D^2)*var_Ky + 2*u.D*v.D*cov_Kx_Ky
  D_v = as.vector(t(D))
  Z.D = rep(0,n*n)
  Z.D[ids] = (D_v[ids]-mean.D[difv[ids]])/sqrt(var.D[difv[ids]])
  Z.D = matrix(abs(Z.D),n,byrow=T)
  Z.D[,(n-1):n] = 0
  Z.Dv = as.vector(t(Z.D))
  Z.D_max = max(Z.Dv[ids2])
  tauhat_0 = which(Z.Dv == Z.D_max)
  tauhat = c(floor(tauhat_0/n)+1, (tauhat_0-1)%%n+1)
  ZD = list(tauhat=tauhat, Zmax=Z.D_max, seq=Z.D)
  scanZ$Z_D = ZD

  # test statistic Z_W
  u.W = t*(t-1)*(n-t)/n
  v.W = (n-t)*(n-t-1)*t/n
  mean.W = mu_Kx*u.W + mu_Ky*v.W
  var.W = var_Kx*u.W^2 + var_Ky*v.W^2 + 2*u.W*v.W*cov_Kx_Ky
  W_v = as.vector(t(W))
  Z.W = rep(0,n*n)
  Z.W[ids] = (W_v[ids]-mean.W[difv[ids]])/sqrt(var.W[difv[ids]])
  Z.W = matrix(Z.W,n,byrow=T)
  Z.W[,(n-1):n] = 0
  Z.Wv = as.vector(t(Z.W))
  Z.W_max = max(Z.Wv[ids2])
  tauhat_0 = which(Z.Wv == Z.W_max)
  tauhat = c(floor(tauhat_0/n)+1, (tauhat_0-1)%%n+1)
  ZW = list(tauhat=tauhat, Zmax=Z.W_max, seq=Z.W)
  scanZ$Z_W = ZW

  # test statistic Z_W_r1
  u.W1 = r1*t*(t-1)*(n-t)/n
  v.W1 = (n-t)*(n-t-1)*t/n
  mean.W1 = mu_Kx*u.W1 + mu_Ky*v.W1
  var.W1 = var_Kx*u.W1^2 + var_Ky*v.W1^2 + 2*u.W1*v.W1*cov_Kx_Ky
  W1_v = as.vector(t(W1))
  Z.W1 = rep(0,n*n)
  Z.W1[ids] = (W1_v[ids]-mean.W1[difv[ids]])/sqrt(var.W1[difv[ids]])
  Z.W1 = matrix(Z.W1,n,byrow=T)
  Z.W1[,(n-1):n] = 0
  Z.W1v = as.vector(t(Z.W1))
  Z.W1_max = max(Z.W1v[ids2])
  tauhat_0 = which(Z.W1v == Z.W1_max)
  tauhat = c(floor(tauhat_0/n)+1, (tauhat_0-1)%%n+1)
  ZW1 = list(tauhat=tauhat, Zmax=Z.W1_max, seq=Z.W1)
  scanZ$Z_W1 = ZW1

  # test statistic Z_W_r2
  u.W2 = r2*t*(t-1)*(n-t)/n
  v.W2 = (n-t)*(n-t-1)*t/n
  mean.W2 = mu_Kx*u.W2 + mu_Ky*v.W2
  var.W2 = var_Kx*u.W2^2 + var_Ky*v.W2^2 + 2*u.W2*v.W2*cov_Kx_Ky
  W2_v = as.vector(t(W2))
  Z.W2 = rep(0,n*n)
  Z.W2[ids] = (W2_v[ids]-mean.W2[difv[ids]])/sqrt(var.W2[difv[ids]])
  Z.W2 = matrix(Z.W2,n,byrow=T)
  Z.W2[,(n-1):n] = 0
  Z.W2v = as.vector(t(Z.W2))
  Z.W2_max = max(Z.W2v[ids2])
  tauhat_0 = which(Z.W2v == Z.W2_max)
  tauhat = c(floor(tauhat_0/n)+1, (tauhat_0-1)%%n+1)
  ZW2 = list(tauhat=tauhat, Zmax=Z.W2_max, seq=Z.W2)
  scanZ$Z_W2 = ZW2

  # GKCP
  u.W = t*(t-1)*(n-t)/n
  v.W = (n-t)*(n-t-1)*t/n
  mean.W = mu_Kx*u.W + mu_Ky*v.W
  var.W = var_Kx*u.W^2 + var_Ky*v.W^2 + 2*u.W*v.W*cov_Kx_Ky
  W_v = as.vector(t(W))
  Z.W = rep(0,n*n)
  Z.W[ids] = (W_v[ids]-mean.W[difv[ids]])/sqrt(var.W[difv[ids]])
  Z.W = matrix(Z.W,n,byrow=T)
  Z.W[,(n-1):n] = 0
  Z.Wv = as.vector(t(Z.W))
  GKCPv = Z.Dv^2 + Z.Wv^2
  GKCP_max = max(GKCPv[ids2])
  tauhat_0 = which(GKCPv == GKCP_max)
  tauhat = c(floor(tauhat_0/n)+1, (tauhat_0-1)%%n+1)
  GKCP = list(tauhat=tauhat, Zmax=GKCP_max, seq=GKCPv)
  scanZ$GKCP = GKCP

  return(scanZ)
}

# the Nu function
Nu = function(x){
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

# p value approximation for single change-point
pval1 = function(n, K, scanZ, r1, r2, skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)) {
  output = list()

  if(lower<2){
    lower = 2
    #print(lower)
  }
  if(upper>(n-2)){
    upper = n-2
    #print(upper)
  }

  t = 1:n
  t = as.numeric(t)

  R0 = sum(K)
  R1 = sum(K^2)
  R2 = sum(rowSums(K)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  if (skew.corr==FALSE) {
    bD = scanZ$Z_D$Zmax
    if (bD>0) {
      integrandD = function(t){
        c1 = rhoD(n, t)
        c1*Nu(sqrt(2*bD^2*c1))
      }
      pval.ZD = 2*dnorm(bD)*bD*integrate(integrandD, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZD = 1
    }
    output$Z_D = min(pval.ZD,1)

    bW1 = scanZ$Z_W1$Zmax
    if (bW1>0) {
      integrandW1 = function(t){
        c1 = rhoW(n, t, r1, R0, R1, R2)
        c1*Nu(sqrt(2*bW1^2*c1))
      }
      pval.ZW1 = dnorm(bW1)*bW1*integrate(integrandW1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW1 = 1
    }
    output$Z_W1 = min(pval.ZW1,1)

    bW2 = scanZ$Z_W2$Zmax
    if (bW2>0) {
      integrandW2 = function(t){
        c1 = rhoW(n, t, r2, R0, R1, R2)
        c1*Nu(sqrt(2*bW2^2*c1))
      }
      pval.ZW2 = dnorm(bW2)*bW2*integrate(integrandW2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW2 = 1
    }
    output$Z_W2 = min(pval.ZW2,1)

    return(output)
  }


  temp = skewcorr(n, K, r1, r2)
  ED.3 = temp$ED.3
  EW1.3 = temp$EW1.3
  EW2.3 = temp$EW2.3

  mu_Kx = R0/n/(n-1)
  mu_Ky = R0/n/(n-1)

  p1 = t*(t-1)/n/(n-1)
  p2 = p1*(t-2)/(n-2)
  p3 = p2*(t-3)/(n-3)

  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = q1*(n-t-2)/(n-2)
  q3 = q2*(n-t-3)/(n-3)

  var_Kx = (2*R1*p1 + 4*R2*p2 + R3*p3)/t/t/(t-1)/(t-1) - mu_Kx^2
  var_Ky = (2*R1*q1 + 4*R2*q2 + R3*q3)/(n-t)/(n-t)/(n-t-1)/(n-t-1) - mu_Ky^2
  cov_Kx_Ky = R3/n/(n-1)/(n-2)/(n-3) - mu_Kx*mu_Ky

  u.D = t*(t-1)
  v.D = -(n-t)*(n-t-1)
  mean.D = mu_Kx*u.D + mu_Ky*v.D
  var.D = (u.D^2)*var_Kx + (v.D^2)*var_Ky + 2*u.D*v.D*cov_Kx_Ky

  u.W1 = r1*t*(t-1)*(n-t)/n
  v.W1 = (n-t)*(n-t-1)*t/n
  mean.W1 = mu_Kx*u.W1 + mu_Ky*v.W1
  var.W1 = var_Kx*u.W1^2 + var_Ky*v.W1^2 + 2*u.W1*v.W1*cov_Kx_Ky

  u.W2 = r2*t*(t-1)*(n-t)/n
  v.W2 = (n-t)*(n-t-1)*t/n
  mean.W2 = mu_Kx*u.W2 + mu_Ky*v.W2
  var.W2 = var_Kx*u.W2^2 + var_Ky*v.W2^2 + 2*u.W2*v.W2*cov_Kx_Ky

  rD =  (ED.3- 3*mean.D*var.D - mean.D^3)/(var.D^(3/2)) # E(Z.D^3)
  rW1 =  (EW1.3- 3*mean.W1*var.W1 - mean.W1^3)/(var.W1^(3/2)) # E(Z.W1^3)
  rW2 =  (EW2.3- 3*mean.W2*var.W2 - mean.W2^3)/(var.W2^(3/2)) # E(Z.W2^3)

  for(i in 1:length(rD)){
    if (is.na(rD[i])==TRUE){
      rD[i]=0
    }
    if (is.na(rW1[i])==TRUE){
      rW1[i]=0
    }
    if (is.na(rW2[i])==TRUE){
      rW2[i]=0
    }
  }
  if (rD[n/2]==0) {rD[n/2]=rD[n/2+1]}

  c1_D = rhoD(n, t)
  for(i in 1:length(c1_D)){
    if ((abs(c1_D[i]))=="Inf"){
      c1_D[i]=0
    }
  }
  c1_W1 = rhoW(n, t, r1, R0, R1, R2)
  for(i in 1:length(c1_W1)){
    if ((abs(c1_W1[i]))=="Inf"){
      c1_W1[i]=0
    }
  }
  c1_W2 = rhoW(n, t, r2, R0, R1, R2)
  for(i in 1:length(c1_W2)){
    if ((abs(c1_W2[i]))=="Inf"){
      c1_W2[i]=0
    }
  }

  bD = scanZ$Z_D$Zmax
  bW1 = scanZ$Z_W1$Zmax
  bW2 = scanZ$Z_W2$Zmax

  result_D = pval1_sub_1(n, bD, rD, c1_D, lower, upper)
  if (!is.numeric(result_D) || result_D==0 ) {
    if (result_D ==0) {
      message("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    message("Z_D: p-value approximation without skewness correction is reported.\n")
    bD = scanZ$Z_D$Zmax
    integrandD = function(t){
      c1 = rhoD(n, t)
      c1*Nu(sqrt(2*bD^2*c1))
    }
    pval.ZD = 2*dnorm(bD)*bD*integrate(integrandD, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    output$Z_D = min(pval.ZD,1)
  } else {
    output$Z_D = min(result_D,1)
  }

  result_W1 = pval1_sub_2(n, bW1, rW1, c1_W1, lower, upper)
  if (is.numeric(result_W1) && result_W1 > 0) {
    output$Z_W1 = min(result_W1, 1)
  } else {
    if (result_W1 ==0) {
      message("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    message("Z_W1: p-value approximation without skewness correction is reported.\n")
    bW1 = scanZ$Z_W1$Zmax
    if (bW1>0) {
      integrandW1 = function(t){
        c1 = rhoW(n, t, r1, R0, R1, R2)
        c1*Nu(sqrt(2*bW1^2*c1))
      }
      pval.ZW1 = dnorm(bW1)*bW1*integrate(integrandW1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW1 = 1
    }
    output$Z_W1 = min(pval.ZW1,1)
  }

  result_W2 = pval1_sub_2(n, bW2, rW2, c1_W2, lower, upper)
  if (is.numeric(result_W2) && result_W2 > 0) {
    output$Z_W2 = min(result_W2, 1)
  } else {
    if (result_W2 ==0) {
      message("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    message("Z_W2: p-value approximation without skewness correction is reported.\n")
    bW2 = scanZ$Z_W2$Zmax
    if (bW2>0) {
      integrandW2 = function(t){
        c1 = rhoW(n, t, r2, R0, R1, R2)
        c1*Nu(sqrt(2*bW2^2*c1))
      }
      pval.ZW2 = dnorm(bW2)*bW2*integrate(integrandW2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW2 = 1
    }
    output$Z_W2 = min(pval.ZW2,1)
  }

  return(output)
}


# p value approximation for single change-point
pval2 = function(n, K, scanZ, r1, r2, skew.corr=TRUE, lower=ceiling(0.05*n), upper=floor(0.95*n)) {
  output = list()

  if(lower<2){
    lower = 2
    #print(lower)
  }
  if(upper>(n-2)){
    upper = n-2
    #print(upper)
  }

  t = 1:n
  t = as.numeric(t)

  R0 = sum(K)
  R1 = sum(K^2)
  R2 = sum(rowSums(K)^2) - R1
  R3 = R0^2 - 2*R1 - 4*R2

  if (skew.corr==FALSE) {
    bD = scanZ$Z_D$Zmax
    if (bD>0) {
      integrandD = function(t){
        c1 = rhoD(n, t)
        (c1*Nu(sqrt(2*bD^2*c1)))^2*(n-t)
      }
      pval.ZD = 2*dnorm(bD)*bD^3*integrate(integrandD, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZD = 1
    }
    output$Z_D = min(pval.ZD,1)

    bW1 = scanZ$Z_W1$Zmax
    if (bW1>0) {
      integrandW1 = function(t){
        c1 = rhoW(n, t, r1, R0, R1, R2)
        (c1*Nu(sqrt(2*bW1^2*c1)))^2*(n-t)
      }
      pval.ZW1 = dnorm(bW1)*bW1^3*integrate(integrandW1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW1 = 1
    }
    output$Z_W1 = min(pval.ZW1,1)

    bW2 = scanZ$Z_W2$Zmax
    if (bW2>0) {
      integrandW2 = function(t){
        c1 = rhoW(n, t, r2, R0, R1, R2)
        (c1*Nu(sqrt(2*bW2^2*c1)))^2*(n-t)
      }
      pval.ZW2 = dnorm(bW2)*bW2^3*integrate(integrandW2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW2 = 1
    }
    output$Z_W2 = min(pval.ZW2,1)

    return(output)
  }


  temp = skewcorr(n, K, r1, r2)
  ED.3 = temp$ED.3
  EW1.3 = temp$EW1.3
  EW2.3 = temp$EW2.3

  mu_Kx = R0/n/(n-1)
  mu_Ky = R0/n/(n-1)

  p1 = t*(t-1)/n/(n-1)
  p2 = p1*(t-2)/(n-2)
  p3 = p2*(t-3)/(n-3)

  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = q1*(n-t-2)/(n-2)
  q3 = q2*(n-t-3)/(n-3)

  var_Kx = (2*R1*p1 + 4*R2*p2 + R3*p3)/t/t/(t-1)/(t-1) - mu_Kx^2
  var_Ky = (2*R1*q1 + 4*R2*q2 + R3*q3)/(n-t)/(n-t)/(n-t-1)/(n-t-1) - mu_Ky^2
  cov_Kx_Ky = R3/n/(n-1)/(n-2)/(n-3) - mu_Kx*mu_Ky

  u.D = t*(t-1)
  v.D = -(n-t)*(n-t-1)
  mean.D = mu_Kx*u.D + mu_Ky*v.D
  var.D = (u.D^2)*var_Kx + (v.D^2)*var_Ky + 2*u.D*v.D*cov_Kx_Ky

  u.W1 = r1*t*(t-1)*(n-t)/n
  v.W1 = (n-t)*(n-t-1)*t/n
  mean.W1 = mu_Kx*u.W1 + mu_Ky*v.W1
  var.W1 = var_Kx*u.W1^2 + var_Ky*v.W1^2 + 2*u.W1*v.W1*cov_Kx_Ky

  u.W2 = r2*t*(t-1)*(n-t)/n
  v.W2 = (n-t)*(n-t-1)*t/n
  mean.W2 = mu_Kx*u.W2 + mu_Ky*v.W2
  var.W2 = var_Kx*u.W2^2 + var_Ky*v.W2^2 + 2*u.W2*v.W2*cov_Kx_Ky

  rD =  (ED.3- 3*mean.D*var.D - mean.D^3)/(var.D^(3/2)) # E(Z.D^3)
  rW1 =  (EW1.3- 3*mean.W1*var.W1 - mean.W1^3)/(var.W1^(3/2)) # E(Z.W1^3)
  rW2 =  (EW2.3- 3*mean.W2*var.W2 - mean.W2^3)/(var.W2^(3/2)) # E(Z.W2^3)

  for(i in 1:length(rD)){
    if (is.na(rD[i])==TRUE){
      rD[i]=0
    }
    if (is.na(rW1[i])==TRUE){
      rW1[i]=0
    }
    if (is.na(rW2[i])==TRUE){
      rW2[i]=0
    }
  }
  if (rD[n/2]==0) {rD[n/2]=rD[n/2+1]}

  c1_D = rhoD(n, t)
  for(i in 1:length(c1_D)){
    if ((abs(c1_D[i]))=="Inf"){
      c1_D[i]=0
    }
  }
  c1_W1 = rhoW(n, t, r1, R0, R1, R2)
  for(i in 1:length(c1_W1)){
    if ((abs(c1_W1[i]))=="Inf"){
      c1_W1[i]=0
    }
  }
  c1_W2 = rhoW(n, t, r2, R0, R1, R2)
  for(i in 1:length(c1_W2)){
    if ((abs(c1_W2[i]))=="Inf"){
      c1_W2[i]=0
    }
  }

  bD = scanZ$Z_D$Zmax
  bW1 = scanZ$Z_W1$Zmax
  bW2 = scanZ$Z_W2$Zmax

  result_D = pval2_sub_1(n, bD, rD, c1_D, lower, upper)
  if (!is.numeric(result_D) || result_D==0 ) {
    if (result_D ==0) {
      message("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    message("Z_D: p-value approximation without skewness correction is reported.\n")
    bD = scanZ$Z_D$Zmax
    integrandD = function(t){
      c1 = rhoD(n, t)
      (c1*Nu(sqrt(2*bD^2*c1)))^2*(n-t)
    }
    pval.ZD = 2*dnorm(bD)*bD^3*integrate(integrandD, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    output$Z_D = min(pval.ZD,1)
  } else {
    output$Z_D = min(result_D,1)
  }

  result_W1 = pval2_sub_2(n, bW1, rW1, c1_W1, lower, upper)
  if (is.numeric(result_W1) && result_W1 > 0) {
    output$Z_W1 = min(result_W1, 1)
  } else {
    if (result_W1 ==0) {
      message("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    message("Z_W1: p-value approximation without skewness correction is reported.\n")
    bW1 = scanZ$Z_W1$Zmax
    if (bW1>0) {
      integrandW1 = function(t){
        c1 = rhoW(n, t, r1, R0, R1, R2)
        (c1*Nu(sqrt(2*bW1^2*c1)))^2*(n-t)
      }
      pval.ZW1 = dnorm(bW1)*bW1^3*integrate(integrandW1, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW1 = 1
    }
    output$Z_W1 = min(pval.ZW1,1)
  }

  result_W2 = pval2_sub_2(n, bW2, rW2, c1_W2, lower, upper)
  if (is.numeric(result_W2) && result_W2 > 0) {
    output$Z_W2 = min(result_W2, 1)
  } else {
    if (result_W2 ==0) {
      message("Extrapolation for skewness-corrected p-value approximation could not be performed. \n")
    }
    message("Z_W2: p-value approximation without skewness correction is reported.\n")
    bW2 = scanZ$Z_W2$Zmax
    if (bW2>0) {
      integrandW2 = function(t){
        c1 = rhoW(n, t, r2, R0, R1, R2)
        (c1*Nu(sqrt(2*bW2^2*c1)))^2*(n-t)
      }
      pval.ZW2 = dnorm(bW2)*bW2^3*integrate(integrandW2, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value
    } else {
      pval.ZW2 = 1
    }
    output$Z_W2 = min(pval.ZW2,1)
  }

  return(output)
}



rhoD = function(n, t) {
  n/2/t/(n-t)
}

# r = r1 or r2
rhoW = function(n, t, r, R0, R1, R2) {
  R3 = R0^2 - 2*R1 - 4*R2
  ft = t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  qhat = (n-t)/n
  phat = t/n

  temp1 = 2*t-n+2*n*t-3*t^2
  temp2 = 3*t^2-4*n*t+2*t+n^2-n

  r1 = ((n-t)/(n^3*(n-1)))*(r^2*temp1 + r*t + t*(n-t-1))
  r2 = ((n-t)/(n^3*(n-1)*(n-2)))*(r^2*temp1*(t-2) + r*t*(t^2-n*t+n-2) + t*(n-t-1)*(n-2*t-2))

  part1 = r^2*temp1*(t-2)*(t-3)
  part2 = r*temp1*(n-t-1)*t
  part3 = r*t*(3*t^3-4*n*t^2-5*t^2+n^2*t+7*n*t-2*t-n^2-3*n+6)
  part4 = t*(n-t-1)*(3*t^2-4*n*t+10*t+n^2-5*n+6)
  r3 = ((n-t)/(n^3*(n-1)*(n-2)*(n-3)))*(part1 + part2 + part3 + part4)

  r0 = (t*(n-t)/(n^4*(n-1)^2))*(r^2*temp1*(t-1) + r*temp1*(n-t-1) + r*temp2*(t-1) + temp2*(n-t-1))

  numerator_deriv = 2*R1*r1 + 4*R2*r2 + R3*r3 - R0^2*r0

  vart1 = ( (r*qhat+phat)^2*(2*R1-2*R0^2/n/(n-1)) + (r^2*qhat^2*(t-2)/(n-t-1)+phat^2*(n-t-2)/(t-1)-2*r*qhat*phat)*(4*R1+4*R2-4*R0^2/n) )
  vart = ft*vart1

  vart_deriv1 = (4*t^3-6*n*t^2+2*n^2*t+2*n*t-2*t-n^2+n)/n/(n-1)/(n-2)/(n-3)
  vart_deriv21 = ( 2*(1-r)*(r*(n-t)+t)/n/n )
  vart_deriv22 = r^2*(2*t^2-3*n*t+t+n^2+n-4)*(n-t)/(n-t-1)/(n-t-1)
  vart_deriv23 = t*(-2*t^2+n*t+t-2*n+4)/(t-1)/(t-1)

  vart_deriv = vart_deriv1*vart1 + ft*( vart_deriv21*(2*R1-2*R0^2/n/(n-1)) + (4*R1+4*R2-4*R0^2/n)*(vart_deriv22+vart_deriv23-2*r*(n-2*t))/n/n )

  res = (2*numerator_deriv - vart_deriv)/2/vart

  return(res)
}


skewcorr = function(n, K, r1=1.2, r2=0.8) {
  t = 1:n

  p1 = t*(t-1)/n/(n-1)
  p2 = p1*(t-2)/(n-2)
  p3 = p2*(t-3)/(n-3)
  p4 = p3*(t-4)/(n-4)
  p5 = p4*(t-5)/(n-5)

  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = q1*(n-t-2)/(n-2)
  q3 = q2*(n-t-3)/(n-3)
  q4 = q3*(n-t-4)/(n-4)
  q5 = q4*(n-t-5)/(n-5)

  f1 = t*(t-1)*(n-t)*(n-t-1)/n/(n-1)/(n-2)/(n-3)
  f2 = f1*(t-2)/(n-4)
  f3 = f2*(t-3)/(n-5)

  g1 = f1
  g2 = g1*(n-t-2)/(n-4)
  g3 = g2*(n-t-3)/(n-5)

  R0 = sum(K)
  R1 = sum(K^2)
  R_temp = rowSums(K)
  R_temp1 = rowSums(K^2)
  R_temp2 = R_temp^2 - R_temp1
  R2 = sum(R_temp2)
  temp1 = sum(K^3)
  temp2 = sum( R_temp*R_temp1 ) - temp1
  temptemp = skew(K, R_temp, R_temp2, R0, R2)
  temp3 = 2*temptemp[1]
  temp4 = 2*temptemp[2]
  temp5 = 2*temptemp[3]
  temp6 = 2*temptemp[4]
  temp7 = 2*temptemp[5]
  temp8 = sum(K)^3 - 4*temp1 - 24*temp2 - 8*temp3 - 6*temp4 - 8*temp5 - 24*temp6 - 12*temp7

  # alpha^3
  EKx3 = ( 4*temp1*p1+24*temp2*p2+8*temp3*p2+6*temp4*p3+8*temp5*p3+24*temp6*p3+12*temp7*p4+temp8*p5 )/t/t/t/(t-1)/(t-1)/(t-1)
  # beta^3
  EKy3 = ( 4*temp1*q1+24*temp2*q2+8*temp3*q2+6*temp4*q3+8*temp5*q3+24*temp6*q3+12*temp7*q4+temp8*q5 )/(n-t)/(n-t)/(n-t)/(n-t-1)/(n-t-1)/(n-t-1)
  # alpha^2beta
  EKx2Ky = ( 2*temp4*f1+4*temp7*f2+temp8*f3 )/t/t/(t-1)/(t-1)/(n-t)/(n-t-1)
  # beta^2alpha
  EKxKy2 = ( 2*temp4*g1+4*temp7*g2+temp8*g3 )/t/(t-1)/(n-t)/(n-t)/(n-t-1)/(n-t-1)

  u.D = t*(t-1)
  v.D = -(n-t)*(n-t-1)
  u.W1 = r1*t*(t-1)*(n-t)/n
  v.W1 = (n-t)*(n-t-1)*t/n
  u.W2 = r2*t*(t-1)*(n-t)/n
  v.W2 = (n-t)*(n-t-1)*t/n

  ED.3 = u.D^3*EKx3 + 3*u.D^2*v.D*EKx2Ky + 3*u.D*v.D^2*EKxKy2 + v.D^3*EKy3
  EW1.3 = u.W1^3*EKx3 + 3*u.W1^2*v.W1*EKx2Ky + 3*u.W1*v.W1^2*EKxKy2 + v.W1^3*EKy3
  EW2.3 = u.W2^3*EKx3 + 3*u.W2^2*v.W2*EKx2Ky + 3*u.W2*v.W2^2*EKxKy2 + v.W2^3*EKy3

  return( list(ED.3 = ED.3, EW1.3 = EW1.3, EW2.3 = EW2.3) )
}

# Support function(skewness correction) for approximated p-value in single change point setting
# approximated p-value with extrapolation for Z.D
pval1_sub_1 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = x*Nu(sqrt(2*b^2*x)) * ratio

  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):n]*b>0))

  if (nn.l>0.35*n){
    return(0)
  }

  if (nn.l>=lower){
    neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  if (nn.r>=(n-upper)){
    neg = which(1+2*r[ceiling(n/2):n]*b<=0)
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:n] = ratio[id2-1]+inc*((id2:n)-id2)
    ratio[ratio<0]=0
    a[(n/2):n] = (x*Nu(sqrt(2*b^2*x))*ratio)[(n/2):n] # update a after extrapolation
  }

  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]
  }

  result = try(2*dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}

# approximated p-value(skewness correction) with extrapolation for Z.W
pval1_sub_2 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b) # S(t) in integrand
  a = x*Nu(sqrt(2*b^2*x)) * ratio
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  if (nn>=(lower-1)+(n-upper)){
    neg = which(1+2*r*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0 | is.na(a)==TRUE)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]
  }

  result = try(dnorm(b)*b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T)
  return(result)
}


# Support function(skewness correction) for approximated p-value in changed interval setting
# approximated p-value with extrapolation for max-count statistic
pval2_sub_1 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  for(i in 1:length(theta_b[pos])){
    if (is.na(theta_b[pos][i])==TRUE){
      theta_b[pos][i]=0
    }
  }
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b)
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2*ratio  ## this is different from single change-point alternative

  nn.l = ceiling(n/2)-length(which(1+2*r[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*r[ceiling(n/2):n]*b>0))

  if (nn.l>0.35*n){
    return(0)
  }

  neg = which(1+2*r[1:ceiling(n/2)]*b<=0)
  if (nn.l>=lower){      ## this part also differs from single change-point
    dif = c(diff(neg),n/2-nn.l)
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
  }
  neg = which(1+2*r[ceiling(n/2):n]*b<=0)
  if (nn.r>=(n-upper)){        ## this part also differs from single change-point
    id1 = min(neg+ceiling(n/2)-1,ceiling(n/2)-1)
    id2 = id1 - ceiling(0.03*n)
    id3 = id2 - ceiling(0.09*n)
    inc = (ratio[id3]-ratio[id2])/(id3-id2)
    ratio[id2:n] = ratio[id2-1]+inc*((id2:n)-id2)
    ratio[ratio<0]=0
    a[(n/2):n] = ((b^2*x*Nu(sqrt(2*b^2*x)))^2*ratio)[(n/2):n]
  }

  neg2 = which(a<0)
  a[neg2] = 0
  integrand = function(s){
    a[s]*(n-s)       # different from single change point
  }

  result = try(2*dnorm(b)/b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T) # differs from single
  return(result)
}

# approximated p-value(skewness correction) with extrapolation for weighted statistic
pval2_sub_2 = function(n,b,r,x,lower,upper){
  theta_b = rep(0,n)
  pos = which(1+2*r*b>0)
  theta_b[pos] = (sqrt((1+2*r*b)[pos])-1)/r[pos]
  ratio = exp((b-theta_b)^2/2 + r*theta_b^3/6)/sqrt(1+r*theta_b) # S(t) in integrand
  a = (b^2*x*Nu(sqrt(2*b^2*x)))^2 * ratio  ## this is different from single change-point
  a_na = which(is.na(a)==TRUE )
  a[a_na] = 0
  nn = n-length(pos)
  if (nn>0.75*n){
    return(0)
  }
  neg = which(1+2*r*b<=0)  ## this is different from single change-point
  if (nn>=(lower-1)+(n-upper)){
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (a[id3]-a[id2])/(id3-id2)
    a[id2:1] = a[id2+1]-inc*(1:id2)
    a[(n/2+1):n] = a[(n/2):1]
    neg2 = which(a<0 | is.na(a)==TRUE)
    a[neg2] = 0
  }
  integrand = function(s){
    a[s]*(n-s)    # different from single change point
  }

  result = try(dnorm(b)/b*integrate(integrand, lower, upper, subdivisions=3000, stop.on.error=FALSE)$value, silent=T) ## this is different from single
  return(result)
}


# p value from permutation for single change point (permutation p-value)
permpval1 = function(n, K, scanZ, r1, r2, B=100, n0=ceiling(0.05*n), n1=floor(0.95*n)) {
  if(n0<2){
    n0=2
  }
  if(n1>(n-2)){
    n1=n-2
  }

  Z.D = Z.W1 = Z.W2 = gkcp = matrix(0,B,n)

  for(b in 1:B) {
    if(b%%1000 ==0) {
      message(b, "permutations completed.\n")
    }
    id0 = sample(1:n, replace = FALSE)    # permute the data

    K_id = K[id0, id0]
    gcpstar = computeSTAT1(n, K_id, r1, r2, n0, n1)

    Z.D[b,] = gcpstar$Z_D$seq
    Z.W1[b,] = gcpstar$Z_W1$seq
    Z.W2[b,] = gcpstar$Z_W2$seq
    gkcp[b,] = gcpstar$GKCP$seq
  }

  output = list()
  p=1-(0:(B-1))/B

  maxZ_D = apply(Z.D[,n0:n1],1,max)
  maxZs_D = sort(maxZ_D)
  ZD.pval.perm = min(1, 2*length(which(maxZs_D>=scanZ$Z_D$Zmax))/B )
  output$Z_D = list(pval=ZD.pval.perm, curve=cbind(maxZs_D,p), maxZs_D=maxZs_D, Z=Z.D)

  maxZ_W1 = apply(Z.W1[,n0:n1],1,max)
  maxZs_W1 = sort(maxZ_W1)
  ZW1.pval.perm = min(1, length(which(maxZs_W1>=scanZ$Z_W1$Zmax))/B )
  output$Z_W1 = list(pval=ZW1.pval.perm, curve=cbind(maxZs_W1,p), maxZs_W1=maxZs_W1, Z=Z.W1)

  maxZ_W2 = apply(Z.W2[,n0:n1],1,max)
  maxZs_W2 = sort(maxZ_W2)
  ZW2.pval.perm = min(1, length(which(maxZs_W2>=scanZ$Z_W2$Zmax))/B )
  output$Z_W2 = list(pval=ZW2.pval.perm, curve=cbind(maxZs_W2,p), maxZs_W2=maxZs_W2, Z=Z.W2)

  maxZ_D = apply(gkcp[,n0:n1],1,max)
  maxZs_D = sort(maxZ_D)
  GKCP.pval.perm = min(1, length(which(maxZs_D>=scanZ$GKCP$Zmax))/B )
  output$GKCP = list(pval=GKCP.pval.perm, curve=cbind(maxZs_D,p), maxZs_D=maxZs_D, Z=gkcp)

  return(output)
}


# p value from permutation for changed interval (permutation p-value)
permpval2 = function(n, K, scanZ, r1, r2, B=100, l0=ceiling(0.05*n), l1=floor(0.95*n)) {
  if(l0<2){
    l0=2
  }
  if(l1>(n-2)){
    l1=n-2
  }

  Z.D.max = Z.W1.max = Z.W2.max = gkcp.max = rep(0,n)

  for(b in 1:B) {
    if(b%%1000 ==0) {
      message(b, "permutations completed.\n")
    }
    id0 = sample(1:n, replace = FALSE)    # permute the data

    K_id = K[id0, id0]
    gcpstar = computeSTAT2(n, K_id, r1, r2, l0, l1)

    Z.D.max[b] = gcpstar$Z_D$Zmax
    Z.W1.max[b] = gcpstar$Z_W1$Zmax
    Z.W2.max[b] = gcpstar$Z_W2$Zmax
    gkcp.max[b] = gcpstar$GKCP$Zmax
  }

  output = list()
  p=1-(0:(B-1))/B

  maxZ_D = Z.D.max
  maxZs_D = sort(maxZ_D)
  ZD.pval.perm = min(1, 2*length(which(maxZs_D>=scanZ$Z_D$Zmax))/B )
  output$Z_D = list(pval=ZD.pval.perm, curve=cbind(maxZs_D,p), maxZs_D=maxZs_D, Z=Z.D.max)

  maxZ_W1 = Z.W1.max
  maxZs_W1 = sort(maxZ_W1)
  ZW1.pval.perm = min(1, length(which(maxZs_W1>=scanZ$Z_W1$Zmax))/B )
  output$Z_W1 = list(pval=ZW1.pval.perm, curve=cbind(maxZs_W1,p), maxZs_W1=maxZs_W1, Z=Z.W1.max)

  maxZ_W2 = Z.W2.max
  maxZs_W2 = sort(maxZ_W2)
  ZW2.pval.perm = min(1, length(which(maxZs_W2>=scanZ$Z_W2$Zmax))/B )
  output$Z_W2 = list(pval=ZW2.pval.perm, curve=cbind(maxZs_W2,p), maxZs_W2=maxZs_W2, Z=Z.W2.max)

  maxZ_D = gkcp.max
  maxZs_D = sort(maxZ_D)
  GKCP.pval.perm = min(1, length(which(maxZs_D>=scanZ$GKCP$Zmax))/B )
  output$GKCP = list(pval=GKCP.pval.perm, curve=cbind(maxZs_D,p), maxZs_D=maxZs_D, Z=gkcp.max)


  return(output)
}


