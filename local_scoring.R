source("main_2tau_deviance.R")

local_scoring <- function(Xall,DDall,Yall,Uall,exampleIdx,lambta){
  
  local_scoring = 1
  nloop = 0
  ## g(E(Y|X,U)) = X^t * beta + f(U)
  # m == 1 
  p = NCOL(Xall)
  q = NCOL(Uall)
  
  lam1 = lambta[1:p]
  lam2 = lambta[(p+1):(p+q)]
  
  Beta_0 = matrix(zeros(1,p)) # Beta_0 = [0,0,0...]
  Ybar = repmat(matrix(sum(Yall)/length(Yall)),NROW(Xall),1)
  fU_0 = functionG(Ybar,exampleIdx) # fU_0 = g(Y_bar)
  eta_0 = Xall%*%Beta_0 + fU_0
  miu_0 = Ybar
  # W_0 is a matrix of N * 1
  W_0 = (1/(functionVar(eta_0,exampleIdx)))*(functionGInversePrime(eta_0,exampleIdx))^2
  W_0_std = W_0/sum(W_0)
  
  Z_sups = eta_0 + (Yall - miu_0) * functionGPrime(miu_0,exampleIdx) # replace y every time
  W_pre_std = W_0_std
  W_sups_std = W_0_std
  epsilon = 1e-5
  
  while (local_scoring) {
    Smat = getLCSmat(DDall, lam2, W_pre_std)
    temp = getHmat(Xall, lam1, Smat, W_pre_std)
    myBeta = temp$Amat %*% temp$Bmat %*% Z_sups
    myfU = Smat %*% (Z_sups - Xall%*%myBeta)
 
    eta_0 = Xall%*%myBeta + myfU
    miu_0 = functionGInverse(eta_0,exampleIdx)
    W_sups = 1/(functionVar(eta_0,exampleIdx)) * (functionGInversePrime(eta_0,exampleIdx))^2
    W_sups_std = W_sups/sum(W_sups)
    Z_sups = eta_0 + (Yall - miu_0) * functionGPrime(miu_0,exampleIdx)
    
    # ### calculate Z_sups_refit
    # Xind=which(lam1>1e-2)
    # Uind=which(lam2>1e-2)
    # 
    # x_ind = Xall[,Xind]
    # dd_ind = DDall[Uind]
    # u_ind = Uall[,Uind] 
    # myBeta_refit = rep(0,p)
    # temp_refit = getHmat_refit(Xall[,Xind], lam1, Smat, W_sups_std)
    # 
    # eta_0_refit = Xall%*%myBeta_refit + Smat %*% (Z_sups_refit - Xall%*%myBeta_refit)
    # miu_0_refit = functionGInverse(eta_0_refit,exampleIdx)
    # Z_sups_refit = eta_0_refit + (Yall - miu_0_refit) * functionGPrime(miu_0_refit,exampleIdx)
    
    
    
    
    ### compare W convergence
    if (any(is.na(W_sups_std))){
      # print(paste("irep:", irep, "W_sups_std NA in LS"))
      print(paste("lambda:", lambta))
    }
    if (any(is.na(W_pre_std))){
      # print(paste("irep:", irep, "W_pre_std NA in LS"))
      print(paste("lambda:", lambta))
    }
    if (max(W_sups_std - W_pre_std )< epsilon){
      Z_sups_final = Z_sups
      W_sups_std_final = W_sups_std
      
      Xind=which(lam1>1e-2)
      Uind=which(lam2>1e-2)
      
      x_ind = Xall[,Xind]
      dd_ind = DDall[Uind]
      u_ind = Uall[,Uind] 
      myBeta_refit = myBeta
      myBeta_refit[which(lam1<1e-2),] = 0
      # temp_refit = getHmat_refit(Xall[,Xind], lam1, Smat, W_sups_std)
      # myBeta_refit[Xind] = temp_refit$Amat %*% temp_refit$Bmat %*% Z_sups
      
      eta_0_refit = Xall%*%myBeta_refit + Smat %*% (Z_sups - Xall%*%myBeta_refit)
      miu_0_refit = functionGInverse(eta_0_refit,exampleIdx)
      Z_sups_refit = eta_0_refit + (Yall - miu_0_refit) * functionGPrime(miu_0_refit,exampleIdx)
      Z_sups_refit_final = Z_sups_refit
      break
    }
    
    if (nloop >= 15){
      Z_sups_final = Z_sups
      W_sups_std_final = W_sups_std
      
      Xind=which(lam1>1e-2)
      Uind=which(lam2>1e-2)
      
      x_ind = Xall[,Xind]
      dd_ind = DDall[Uind]
      u_ind = Uall[,Uind] 
      myBeta_refit = myBeta
      myBeta_refit[which(lam1<1e-2),] = 0
      # myBeta_refit = rep(0,p)
      # temp_refit = getHmat_refit(Xall[,Xind], lam1, Smat, W_sups_std)
      # myBeta_refit[Xind] = temp_refit$Amat %*% temp_refit$Bmat %*% Z_sups
      
      eta_0_refit = Xall%*%myBeta_refit + Smat %*% (Z_sups - Xall%*%myBeta_refit)
      miu_0_refit = functionGInverse(eta_0_refit,exampleIdx)
      
      Z_sups_refit = eta_0_refit + (Yall - miu_0_refit) * functionGPrime(miu_0_refit,exampleIdx)
      Z_sups_refit_final = Z_sups_refit
      print('Break in local scoring while loop after more than 15 iterations')
      break
    }
    
    W_pre_std = W_sups_std
    nloop = nloop + 1
  }
  
  return (list(Z_sups = Z_sups_final,W_0 = W_sups_std_final,Z_sups_refit = Z_sups_refit_final))

}



local_scoring_refit <- function(Xall,DDall,Yall,Uall,exampleIdx,lambta){
  
  local_scoring = 1
  nloop = 0
  ## g(E(Y|X,U)) = X^t * beta + f(U)
  # m == 1 
  p = NCOL(Xall)
  q = NCOL(Uall)
  
  lam1 = lambta[1:p]
  lam2 = lambta[(p+1):(p+q)]
  
  Beta_0 = matrix(zeros(1,p)) # Beta_0 = [0,0,0...]
  Ybar = repmat(matrix(sum(Yall)/length(Yall)),NROW(Xall),1)
  fU_0 = functionG(Ybar,exampleIdx) # fU_0 = g(Y_bar)
  eta_0 = Xall%*%Beta_0 + fU_0
  miu_0 = Ybar
  # W_0 is a matrix of N * 1
  W_0 = (1/(functionVar(eta_0,exampleIdx)))*(functionGInversePrime(eta_0,exampleIdx))^2
  W_0_std = W_0/sum(W_0)
  
  Z_sups = eta_0 + (Yall - miu_0) * functionGPrime(miu_0,exampleIdx) # replace y every time
  W_pre_std = W_0_std
  W_sups_std = W_0_std
  epsilon = 1e-5
  
  while (local_scoring) {
    Smat = getLCSmat(DDall, lam2, W_pre_std)
    temp = getHmat_refit(Xall, lam1, Smat, W_pre_std)
    myBeta = temp$Amat %*% temp$Bmat %*% Z_sups
    myfU = Smat %*% (Z_sups - Xall%*%myBeta)
    
    eta_0 = Xall%*%myBeta + myfU
    miu_0 = functionGInverse(eta_0,exampleIdx)
    W_sups = 1/(functionVar(eta_0,exampleIdx)) * (functionGInversePrime(eta_0,exampleIdx))^2
    W_sups_std = W_sups/sum(W_sups)
    Z_sups = eta_0 + (Yall - miu_0) * functionGPrime(miu_0,exampleIdx)
    ### compare W convergence
    if (any(is.na(W_sups_std))){
      # print(paste("irep:", irep, "W_sups_std NA in LS"))
      print(paste("lambda:", lambta))
    }
    if (any(is.na(W_pre_std))){
      # print(paste("irep:", irep, "W_pre_std NA in LS"))
      print(paste("lambda:", lambta))
    }
    if (max(W_sups_std - W_pre_std )< epsilon){
      Z_sups_final = Z_sups
      W_sups_std_final = W_sups_std
      myBeta_final = myBeta
      mydf_final = temp$df
      break
    }
    
    if (nloop >= 15){
      Z_sups_final = Z_sups
      W_sups_std_final = W_sups_std
      myBeta_final = myBeta
      mydf_final = temp$df
      print('Break in local scoring while loop after more than 15 iterations')
      break
    }
    
    W_pre_std = W_sups_std
    nloop = nloop + 1
  }
  
  return (list(Z_sups = Z_sups_final,W_0 = W_sups_std_final,myBeta = myBeta_final,mydf = mydf_final))
  
}



local_scoring_LL <- function(X,Yall,U,exampleIdx,lam2){
  
  local_scoring = 1
  nloop = 0
  
  X = as.matrix(X)
  U = as.matrix(U)
  
  ## g(E(Y|X,U)) = X^t * beta + f(U)
  # m == 1 
  n = nrow(X)

  p = NCOL(X)
  q = NCOL(U)
  
  Beta_0 = matrix(zeros(1,p)) # Beta_0 = [0,0,0...]
  Ybar = repmat(matrix(sum(Yall)/length(Yall)),nrow(X),1)
  fU_0 = functionG(Ybar,exampleIdx) # fU_0 = g(Y_bar)
  eta_0 = X%*%Beta_0 + fU_0
  miu_0 = Ybar
  # W_0 is a matrix of N * 1
  W_0 = (1/(functionVar(eta_0,exampleIdx)))*(functionGInversePrime(eta_0,exampleIdx))^2
  W_0_std = W_0/sum(W_0)
  
  Z_sups = eta_0 + (Yall - miu_0) * functionGPrime(miu_0,exampleIdx) # replace y every time
  W_pre_std = W_0_std
  W_sups_std = W_0_std
  epsilon = 1e-5
  
  while (local_scoring) {
    Smat=LLsmoothingFAST(U,U,Z_sups,lam2,W_pre_std)$LLSmat
    diag_W = diag(as.numeric(sqrt(W_pre_std)),names = TRUE) 
    
    Imat=eye(n)
    temp=((Imat-Smat)%*%X)
    Amat = solve(t(temp)%*% diag_W %*%temp)
    Bmat = t(temp)%*% diag_W %*%(Imat-Smat)
    Hmat=(Imat-Smat)%*%(Imat- X%*%Amat%*%Bmat)
    myBeta = Amat %*% Bmat %*% Z_sups
    myfU = Smat %*% (Z_sups - X%*%myBeta)
    
    
    
    
    
    eta_0 = X%*%myBeta + myfU
    miu_0 = functionGInverse(eta_0,exampleIdx)
    W_sups = 1/(functionVar(eta_0,exampleIdx)) * (functionGInversePrime(eta_0,exampleIdx))^2
    W_sups_std = W_sups/sum(W_sups)
    Z_sups = eta_0 + (Yall - miu_0) * functionGPrime(miu_0,exampleIdx)
    ### compare W convergence
    # if (any(is.na(W_sups_std))){
    #   # print(paste("irep:", irep, "W_sups_std NA in LSLL"))
    #   print(paste("lambda:", lambta))
    # }
    # if (any(is.na(W_pre_std))){
    #   print(paste("irep:", irep, "W_pre_std NA in LSLL"))
    #   print(paste("lambda:", lambta))
    # }
    if (max(W_sups_std - W_pre_std )< epsilon){
      Z_sups_final = Z_sups
      W_sups_std_final = W_sups_std
      break
    }
    
    if (nloop >= 15){
      Z_sups_final = Z_sups
      W_sups_std_final = W_sups_std
      print('Break in LL local scoring while loop after more than 15 iterations')
      break
    }
    
    W_pre_std = W_sups_std
    nloop = nloop + 1
  }
  
  return (list(Z_sups = Z_sups_final,W_0 = W_sups_std_final))
  
}

