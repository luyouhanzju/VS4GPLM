
source("examples.R")

getLCSmat<-function(DD, lambda,W_0) {
  n=dim(DD[[1]])[1]
  if(length(lambda)!=length(DD)) {
    stop('THERE IS A ERROR: DD and lambda does not match')
  }
  q=length(lambda)
  
  tempDUMP=sapply(1:q,function(t) DD[[t]] * lambda[t], simplify=FALSE)
  temp=Reduce("+",tempDUMP)
  W_mat = matrix( rep( W_0  , n ) , ncol =  n , byrow = FALSE )
  K = exp(-temp^2) * (W_mat)
  #K = exp(-temp^2)
  InvKrowSum=1/apply(K,1,sum)
  Smat=t(sapply(1:n, function(j) K[j,]*InvKrowSum[j], simplify=TRUE))
  return(Smat)
}



getHmat<-function(x, lambdaX, Smat,W_0) {
  # Here we use WLS (20220406)
  x=as.matrix(x)
  n=nrow(Smat)
  p=NCOL(x)
  #XlamHalf=x%*%diag(sqrt(lambdaX))
  
  # diag_W is a N * N matrix 
  diag_W = diag(as.numeric(sqrt(W_0)),names = TRUE) 
  
  RecipN=1/n
  
  Imat=eye(n)
  
  ImS=Imat-Smat
  ImSx=ImS%*%x
  if (length(lambdaX) > 1)
    {
  DsqrtL=diag(sqrt(lambdaX))}
  
  else{
    
    DsqrtL = sqrt(lambdaX)
      
  }
  
  Bmat=RecipN*(t(ImSx)%*% diag_W %*%ImS) # Dim: 5 * 100
  #Bmat=RecipN*(t(ImSx) %*%ImS)
  temp=ImSx%*%DsqrtL
  Amat=DsqrtL%*%solve(RecipN*( t(temp) %*% diag_W %*% temp)+eye(p))%*%DsqrtL
  #Amat=DsqrtL%*%solve(RecipN*( t(temp) %*% temp)+eye(p))%*%DsqrtL
  Hmat=ImS%*%(Imat-x%*%Amat%*%Bmat) # This Hmart corresponds to Chapter 3, R(gamma, lambta)
  # IMPORTANT 2022/04/03
  # BetaMat = Amat %*% Bmat %*%Y
  # The df orginally = n - tr(Hmat)
  # df=n-tr(Hmat)
  df = sum(lambdaX > 0) + tr(Smat)
  return(list(Hmat=Hmat, df=df,Amat = Amat,Bmat = Bmat))
}





objfun<-function(X, DD, Z_sups, lambda, W_0, exampleIdx,Yall) {
  
  p=NCOL(X)
  q=length(DD)
  if(length(lambda)!=(p+q)) {
    stop('Somehting wrong: dimension does not match!!!')
  }  
  lambdaX=lambda[1:p]
  lambdaZ=lambda[p+(1:q)]
  
  Smat=getLCSmat(DD, lambdaZ,W_0)
  temp = getHmat(X, lambdaX, Smat,W_0)
  Hmat = temp$Hmat
  betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
  ## Use likelihood in obj instead; 
  etaAll = X%*%betaMat + Smat%*% (Z_sups - X%*%betaMat)
  res = -logLikelihood(Yall,etaAll,exampleIdx)

  #return(sum(res^2))
  return(res)
}



getHmat_refit<-function(x, lambdaX, Smat,W_0) {
  # Here we use WLS (20220406)
  x=as.matrix(x)
  n=nrow(Smat)
  p=NCOL(x)
  #XlamHalf=x%*%diag(sqrt(lambdaX))
  
  # diag_W is a N * N matrix 
  diag_W = diag(as.numeric(sqrt(W_0)),names = TRUE) 
  
  RecipN=1/n
  
  Imat=eye(n)
  
  ImS=Imat-Smat
  ImSx=ImS%*%x
  # if (length(lambdaX) > 1)
  # {
  #   DsqrtL=diag(sqrt(lambdaX))}
  # 
  # else{
  #   
  #   DsqrtL = sqrt(lambdaX)
  #   
  # }
  
  Bmat=(t(ImSx)%*% diag_W %*%ImS) # Dim: 5 * 100
  #Bmat=RecipN*(t(ImSx) %*%ImS)
  temp=ImSx#%*%DsqrtL
  Amat=solve((t(ImSx) %*% diag_W %*% ImSx))
  #Amat=DsqrtL%*%solve(RecipN*( t(temp) %*% temp)+eye(p))%*%DsqrtL
  Hmat=ImS%*%(Imat-x%*%Amat%*%Bmat) # This Hmart corresponds to Chapter 3, R(gamma, lambta)
  
  
  
  # IMPORTANT 2022/04/03
  # BetaMat = Amat %*% Bmat %*%Y
  # The df orginally = n - tr(Hmat)
  # df=n-tr(Hmat)
  df = sum(lambdaX > 0) + tr(Smat)
  return(list(Hmat=Hmat, df=df,Amat = Amat,Bmat = Bmat))
}

objfun_refit<-function(x_ind, dd_ind, Z_sups, lambda, W_0, exampleIdx,Yall) {
  
  p=NCOL(x_ind)
  q=length(dd_ind)
  if(length(lambda)!=(p+q)) {
    stop('Somehting wrong: dimension does not match!!!')
  }  
  lambdaX=lambda[1:p]
  lambdaZ=lambda[p+(1:q)]
  
  Smat=getLCSmat(dd_ind, lambdaZ,W_0)
  temp = getHmat_refit(x_ind, lambdaX, Smat,W_0)
  Hmat = temp$Hmat
  betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
  ## Use likelihood in obj instead; 
  etaAll = x_ind%*%betaMat + Smat%*% (Z_sups - x_ind%*%betaMat)
  res = -logLikelihood(Yall,etaAll,exampleIdx)

  #return(sum(res^2))
  return(res)
}


## LLobjfun(X, U, Z_sups, lambdaCUR,W_sups_std,exampleIdx，Yall)$MinusLoglikelihood
LLobjfun<-function(x, Uall, Z_sups, lambdaZ,W_sups_std,exampleIdx,Yall) {
  # this version does not apply Ridge penalty
  # if want to use ridge penalty, adjust the above codes for local constant smoothing
  x = as.matrix(x)
  Uall = as.matrix(Uall)
  
  n=nrow(x)
  
  p=NCOL(x)
  q=NCOL(Uall)
  
  ##
  # p=NCOL(X)
  # q=length(DD)
  # if(length(lambda)!=(p+q)) {
  #   stop('Somehting wrong: dimension does not match!!!')
  # }  
  # lambdaX=lambda[1:p]
  # lambdaZ=lambda[p+(1:q)]
  # 
  # Smat=getLCSmat(DD, lambdaZ,W_0)
  # temp = getHmat(X, lambdaX, Smat,W_0)
  # Hmat = temp$Hmat
  # betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
  # ## Use likelihood in obj instead; 
  # etaAll = X%*%betaMat + Smat%*% (Z_sups - X%*%betaMat)
  # res = -logLikelihood(Yall,etaAll,exampleIdx)
  Smat=LLsmoothingFAST(Uall,Uall,Z_sups,lambdaZ,W_sups_std)$LLSmat
  diag_W = diag(as.numeric(sqrt(W_sups_std)),names = TRUE) 
  
  # betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
  # etaAll = X%*%betaMat + Smat%*% (Z_sups - X%*%betaMat)
  # res = -logLikelihood(Yall,etaAll,exampleIdx)
  if(p>0.5){
    Imat=eye(n)
    temp=((Imat-Smat)%*%x)
    Amat = solve(t(temp)%*% diag_W %*%temp)
    Bmat = t(temp)%*% diag_W %*%(Imat-Smat)
    Hmat=(Imat-Smat)%*%(Imat- x%*%Amat%*%Bmat)
    
    res=Hmat%*%Z_sups
    RSS=sum(res^2)
    df=n-tr(Hmat)
  }else{
    df=tr(Smat)
    res=Z_sups-Smat%*%Z_sups
    RSS=sum(res^2)
  }
  
  betaMat = Amat %*% Bmat %*% Z_sups
  ## Use likelihood in obj instead; 
  etaAll = x%*%betaMat + Smat %*% (Z_sups - x%*%betaMat)
  res = -logLikelihood(Yall,etaAll,exampleIdx)
  
  return(list(MinusLoglikelihood = res, df=df,RSS = RSS))
}



PredictLLobjfun<-function(x, z, y, lambdaZ, xtest, ztest, ytest) {
  # this version does not apply Ridge penalty
  # if want to use ridge penalty, adjust the above codes for local constant smoothing
  x=as.matrix(x)
  z=as.matrix(z)
  n=nrow(x)
  p=NCOL(x)
  q=NCOL(z)
  
  xtest=as.matrix(xtest)
  ztest=as.matrix(ztest)
  
  #  Smat=LLsmoothing(z,z,y,lambdaZ)$LLSmat
  Smat = LLsmoothingFAST(z,z,y,lambdaZ)$LLSmat
  testSmat=LLsmoothingFAST(z,ztest,y,lambdaZ)$LLSmat
  
  if(p>0.5){
    Imat=eye(n)
    temp=((Imat-Smat)%*%x)
    Hmat=(Imat-Smat)%*%(Imat- x%*%solve(t(temp)%*%temp)%*%(t(temp)%*%(Imat-Smat)))
    res=Hmat%*%y
    RSS=sum(res^2)
    df=n-tr(Hmat)
    
    betahat=(solve(t(temp)%*%temp)%*%(t(temp)%*%(Imat-Smat)))%*%y
    
    ## RSS is same as  sum((y-x%*%betahat-Smat%*%(y-x%*%betahat))^2)
    
    testyhat=xtest%*%betahat+testSmat%*%(y-x%*%betahat)
  }else{
    df=tr(Smat)
    res=y-Smat%*%y
    RSS=sum(res^2)
    
    testyhat=testSmat%*%y
  }
  testres=ytest-testyhat
  testRSS=mean(testres^2)
  
  return(list(RSS=RSS, df=df, testyhat=testyhat, 
              testRSS=testRSS, testSmat=testSmat, Smat=Smat, betahat=betahat))
  
}






#   Smat=LLsmoothingFAST(Uall,Uall,Z_sups,lambdaZ,W_sups_std)$LLSmat

# tempDUMP=sapply(1:q,function(t) DD[[t]] * lambda[t], simplify=FALSE)
# temp=Reduce("+",tempDUMP)
# W_mat = matrix( rep( W_0  , n ) , ncol =  n , byrow = FALSE )
# K = exp(-temp^2) * (W_mat)
# #K = exp(-temp^2)
# InvKrowSum=1/apply(K,1,sum)
# Smat=t(sapply(1:n, function(j) K[j,]*InvKrowSum[j], simplify=TRUE))

LLsmoothing<-function(z,ztest,Z_sups,lam) {
  
  z=as.matrix(z)
  ztest=as.matrix(ztest)
  n=nrow(z)
  q=NCOL(z)
  LLSmat=NULL
  nt=nrow(ztest)
  for (i in 1:nt) {
    
    temp=z-repmat(ztest[i,],n,1)
    XX=cbind(rep(1,n), temp)
    
    Kmat=exp(-apply((temp*repmat(lam, n, 1))^2,1,sum))
    # Kmat2=exp(-apply((temp%*%diag(lam))^2,1,sum))
    
    
    tempDUMP=sapply(1:n,function(s) (XX[s,]%*%t(XX[s,]))*Kmat[s], simplify=FALSE)
    temp=Reduce("+",tempDUMP)  ##checking#  t(XX)%*%diag(Kmat)%*%((XX))
    
    #  temp2=t(XX)%*%diag(Kmat2)%*%XX
    #  s2=(solve(temp2)%*%(t(XX)%*%diag(Kmat2)))[1,]
    
    s=(solve(temp)%*%t(XX*repmat(Kmat,1,q+1)))[1,]
    
    LLSmat=rbind(LLSmat, s)
  }
  
  yhat=LLSmat %*% Z_sups
  
  return(list(yhat=yhat, LLSmat=LLSmat))
}
#plot(y, yhat,type='p')


#   Smat=LLsmoothingFAST(Uall,Uall,Z_sups,lambdaZ,W_sups_std)$LLSmat

# tempDUMP=sapply(1:q,function(t) DD[[t]] * lambda[t], simplify=FALSE)
# temp=Reduce("+",tempDUMP)
# W_mat = matrix( rep( W_0  , n ) , ncol =  n , byrow = FALSE )
# K = exp(-temp^2) * (W_mat)
# #K = exp(-temp^2)
# InvKrowSum=1/apply(K,1,sum)
# Smat=t(sapply(1:n, function(j) K[j,]*InvKrowSum[j], simplify=TRUE))
# Uall,Uall,Z_sups,lambdaZ,W_sups_std
LLsmoothingFAST<-function(z,ztest,Z_sups,lam,W_sups_std) {
  
  z=as.matrix(z)
  ztest=as.matrix(ztest)
  n=nrow(z)
  q=NCOL(z)
  LLSmat=NULL
  nt=nrow(ztest)
  
  Zdiffztest=sapply(1:nt, function(j) sapply(1:n, function(s) z[s,]-ztest[j,]), simplify=FALSE)
  #W_mat = matrix( rep( W_sups_std  , nt ) , ncol =  nt , byrow = FALSE )
  
  if(q>1.5) {
    getLLsvec<-function(i){
      temp=t(Zdiffztest[[i]])
      XX=cbind(1, temp)
      Kmat= exp(-apply((sapply(1:q, function(s) temp[,s]*lam[s])^2),1,sum)) * W_sups_std
      tempDUMP=sapply(1:n,function(s) (XX[s,]%*%t(XX[s,]))*Kmat[s], simplify=FALSE)
      temp=Reduce("+",tempDUMP)  ##checking#  t(XX)%*%diag(Kmat)%*%((XX))
      s=(solve(temp)%*%t(XX*repmat(Kmat,1,q+1)))[1,]
      return(s)
    }
  }else{
    getLLsvec<-function(i){
      temp=as.matrix(Zdiffztest[[i]])
      XX=cbind(1, temp)
      Kmat= exp(-apply((sapply(1:q, function(s) temp[,s]*lam[s])^2),1,sum))
      Kmat=exp(-(temp*lam)^2) * W_sups_std
      tempDUMP=sapply(1:n,function(s) (XX[s,]%*%t(XX[s,]))*Kmat[s], simplify=FALSE)
      temp=Reduce("+",tempDUMP)  ##checking#  t(XX)%*%diag(Kmat)%*%((XX))
      s=(solve(temp)%*%t(XX*repmat(Kmat,1,q+1)))[1,]
      return(s)
    }
  }
  
  LLSmat =   t(sapply(1:nt, function(i) getLLsvec(i)))
  for (i in (1:nt)){
    getLLsvec(i)
  }
  yhat=LLSmat %*% Z_sups
  
  return(list(yhat=yhat, LLSmat=LLSmat))
}
#plot(y, yhat,type='p')



#####################################################################################################

