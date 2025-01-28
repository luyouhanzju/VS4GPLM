

functionG <- function(t,exampleIdx) # done
{
  
  if (exampleIdx == 1 || exampleIdx ==2){
  
  g = log(matrix(t))
  
  return(g)
  
  }
  
  else if (exampleIdx >= 3){
    
    g = matrix(log(t/(1-t)))
    
  }
  else{ g = matrix(t)}
  
  }

functionGInversePrime <- function(t,exampleIdx) #done
{
  if (exampleIdx == 1 || exampleIdx ==2){
    
    gInversePrime = matrix(exp(t))
    
  }
  
  else if (exampleIdx >= 3){
    
    gInversePrime = matrix(exp(t)/(1+exp(t))^2)
    
  }
  else{ gInversePrime = matrix(repmat(1,length(t),1)) }
  return(gInversePrime)
  
  
}

functionGPrime <- function(t,exampleIdx) # done
{
  if (exampleIdx == 1 || exampleIdx ==2){
    
    gPrime = matrix(1/t)
    
  }
  
  else if (exampleIdx >= 3){
    
    gPrime = matrix(1/(t*(1-t)))
    
  }
  else{gPrime = matrix(repmat(1,length(t),1))}
  return(gPrime)
  
}

functionGInverse <- function(t,exampleIdx) #done
{
  if (exampleIdx == 1 || exampleIdx ==2){
    
    gInverse = matrix(exp(t))
    
  }
  
  else if (exampleIdx >= 3){
    
    gInverse = matrix(exp(t)/(1+exp(t)))
    
  }
  else{gInverse = matrix(t)}
  
  return(gInverse)
}

functionVar <- function(etaAll,exampleIdx) #done
{
  
  if (exampleIdx == 1 || exampleIdx ==2){
    
    miuAll = functionGInverse(etaAll,exampleIdx)
    v = matrix(miuAll)
    
  }
  
  else if (exampleIdx >= 3){
    
    miuAll = functionGInverse(etaAll,exampleIdx)
    
    v = matrix(miuAll * (1-miuAll))
    
  }
  else{ v = matrix(repmat(2,length(etaAll),1))}
  return(v)
  
}




logLikelihood <- function (Yall,etaAll,exampleIdx)
{
  if (exampleIdx == 1 || exampleIdx ==2){
    
    miuAll = functionGInverse(etaAll,exampleIdx)
    Yakk_matrix = matrix(Yall)
    log_N_fac = rep(100000, length(Yall))

    for (junk in 1:length(Yall)){
      if (Yall[junk] == 0) {
        log_N_fac[junk] = 1
      }
      else {
      log_N_fac[junk] = sum(seq(1,Yall[junk],1))
      }
    }
    logLikelihoodRes = apply( (Yall * log(miuAll) - miuAll - log_N_fac),2,sum)
    
  }
  
  else if (exampleIdx >= 3){
    
    miuAll = functionGInverse(etaAll,exampleIdx)
    
    logLikelihoodRes = apply( (Yall * log(miuAll) + (1-Yall) * log(1-miuAll) ),2,sum)
    
  }
  else{    
    n = length(Yall)
    logLikelihoodRes = (4*pi)^(-n/2) * exp(-1/4 * sum((Yall - etaAll)^2))
    
    }
  
  return(logLikelihoodRes)
  
}

gendata <- function (exampleIdx){
  if (exampleIdx == 1){
  n = 100
  p = 5
  q = 5
  r = 0.5
  d = p+q
  rawpred = matrix(rnorm(n*d),n,d)%*%chol(r^(as.matrix(dist(1:d))))
  x = rawpred[,1:p]
  u = 2*pnorm(rawpred[,p+(1:q)])-1
  betacoef = rep(0,p)
  betacoef[2] = 1
  betacoef[4] = 1

  eta = x%*%betacoef + sqrt(2) * sin(pi*u[,1] * u[,3])
  miu = functionGInverse(eta,exampleIdx)
  y = matrix(rpois(n,miu))
  return(list(x = x, u = u, y = y,eta = eta,miu = miu))
  
  }
  else if (exampleIdx == 2) {  
    n = 100
    p = 5
    q = 4
    r = 0.5
    d = p+q
    rawpred = matrix(rnorm(n*d),n,d)%*%chol(r^(as.matrix(dist(1:d))))
    x = rawpred[,1:p]
    u = 2*pnorm(rawpred[,p+(1:q)])-1
    betacoef = rep(0,p)
    betacoef[2] = 1
    betacoef[4] = 1
    
    #eta = x%*%betacoef + 1.5 * sin(pi*u[,1]) + 1.5*(u[,3]-1)^2
    eta = x%*%betacoef + 2 * sin(pi*u[,1]) + 1.5*(u[,3]-1)^2
    miu = functionGInverse(eta,exampleIdx)
    y = matrix(rpois(n,miu))
    #y=f+sqrt(2)*rnorm(n)
    #y = f + rnorm(n)
    return(list(x = x, u = u, y = y,eta = eta,miu = miu))
  }
  else if (exampleIdx == 3) {
    n = 300
    p = 6
    q = 4
    r = 0.5
    d = p+q
    rawpred = matrix(rnorm(n*d),n,d)%*%chol(r^(as.matrix(dist(1:d))))
    x = rawpred[,1:p]
    u = 2*pnorm(rawpred[,p+(1:q)])-1
    betacoef = rep(0,p)
    betacoef[1] = 1
    betacoef[4] = 1
    
    #eta = 1.5*x%*%betacoef + 1.2*(u[,3]-1)^2
    #eta = 1.3*x%*%betacoef + 3* sin(pi*u[,1] * u[,3]) # TEST_bic_GPLM_Example3_%d_170_20231010_newBIC_tau1_refitting_n300_U13_1.Rdata
    eta = 1.6*x%*%betacoef + 4* (u[,1] * u[,3] -1)^2 # TEST_bic_GPLM_Example3_%d_170_20231010_newBIC_tau1_refitting_n300_U13_2.Rdata
    
    miu = functionGInverse(eta,exampleIdx)
    
    y = matrix(rbinom(n,1,miu))
    
    mea = mean(y)
    
    return(list(x = x, u = u, y = y,eta = eta, mea = mea))
  }
  
  else if (exampleIdx == 4)
  {
    data = read.table('~/VSGeneralizedPartialLinearModel/dataset/bankrupcy_1by1.csv', header = TRUE, sep = ',')
    
    x = as.matrix(data[, 1:3])
    u = as.matrix(data[, c(7,11:13,15)])
    y = data[, 16]
     
    return(list(x = x, u = u, y = y))
  }
  
  else if (exampleIdx == 5)
  {
    data = read.table('~/VSGeneralizedPartialLinearModel/dataset/bankrupcy_stratified.csv', header = TRUE, sep = ',')
    
    x = as.matrix(data[, 1:3])
    u = as.matrix(data[, c(7,11:13,15)])
    y = data[, 16]
    
    return(list(x = x, u = u, y = y))
  }
  
  else if (exampleIdx == 6)
  {
    data = read.table('~/VSGeneralizedPartialLinearModel/dataset/diabetes.csv', header = TRUE, sep = ',')
    
    x = as.matrix(data[, 1:3])
    u = as.matrix(data[, 4:8])
    y = data[, 9]
    
    return(list(x = x, u = u, y = y))
  }
  
  else if (exampleIdx == 7)
  {
    data = read.table('~/VSGeneralizedPartialLinearModel/dataset/diabetes_full.csv', header = TRUE, sep = ',', nrows = 768)
    
    x = as.matrix(data[, 1:3])
    u = as.matrix(data[, 4:8])
    y = data[, 9]
    
    return(list(x = x, u = u, y = y))
  }
  
  else if (exampleIdx == 8)
  {
    load("~/PartialLinearModelVS_Wu/TaiWanHousePrice/TaiwanHousePrice.Rdata")
    
    meanYall = mean(Yall)
    y = NULL
    for (junk in (1:length(Yall))){
      if (Yall[junk]>= meanYall){
        y[junk] = 1
      }
      else{
        y[junk] = 0
      }
    }
    
    return(list(x = as.matrix(Xall), u = Zall, y = y))
  }
  
  else if (exampleIdx == 9)
  {
    load("~/VSGeneralizedPartialLinearModel/dataset/framingham.rdata")
    
    # meanYall = mean(Yall)
    # y = NULL
    # for (junk in (1:length(Yall))){
    #   if (Yall[junk]>= meanYall){
    #     y[junk] = 1
    #   }
    #   else{
    #     y[junk] = 0
    #   }
    # }
    
    return(list(x = as.matrix(Xall), u = Zall, y = Yall))
  }
  
  
  else if (exampleIdx == 10)
  {
    load("~/VSGeneralizedPartialLinearModel/dataset/diabetes.rdata")
    
    # meanYall = mean(Yall)
    # y = NULL
    # for (junk in (1:length(Yall))){
    #   if (Yall[junk]>= meanYall){
    #     y[junk] = 1
    #   }
    #   else{
    #     y[junk] = 0
    #   }
    # }
    
    return(list(x = as.matrix(Xall), u = Zall, y = Yall))
  }
  
  
  
}
