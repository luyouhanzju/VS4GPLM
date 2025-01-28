library("matlab")
library(psych)
source("main_2tau_deviance.R")
source("MCD.R")

naiveBICFULL<-function(tauall_1,tauall_2, Xall, DDall,Uall,exampleIdx,Yall) {
  n = length(Xall)
  p = NCOL(Xall)
  q = NCOL(Uall)
  mydf = rep(-1000,length(tauall_1)*length(tauall_2))
  myRSS = rep(-10000, length(tauall_1)*length(tauall_2))
  res = rep(-10000, length(tauall_1)*length(tauall_2))
  # myBeta: 5 * size(tau) 
  myBeta = matrix(rep(-1000,length(tauall_1)*length(tauall_2)*p),ncol = p) 
  # myfU: N * size(tau)
  myfU = matrix(rep(-1000,n * length(tauall_1)*length(tauall_2) ),ncol = n) 
    
  lamhatALL = NULL
  W_sups_std = NULL
  Z_sups_all = NULL
  
  for (j in 1:(length(tauall_1)*length(tauall_2))) {
    
    tau_1 = tauall_1[(j-1)%/%length(tauall_2)+1]
    tau_2 = tauall_2[(j-1)%%length(tauall_2)+1]
    
    result_MCD = MCD(Xall, DDall, tau_1,tau_2,exampleIdx,Yall,Uall)
    
    lamhat = result_MCD$lamhat
    Z_sups = result_MCD$Z_sups
    W_0 = result_MCD$W_0
    
    lamhatALL = rbind(lamhatALL,lamhat)
    W_sups_std = cbind(W_sups_std,W_0)
    Z_sups_all = cbind(Z_sups_all,Z_sups)
    
    lam1 = lamhat[1:p]
    lam2 = lamhat[p+(1:q)]
    
    if(sum(lam2)>1e-3) {
      Smat = getLCSmat(DDall, lam2, W_0)
      if(sum(lam1)>1e-3) {
        temp = getHmat(Xall, lam1, Smat, W_0)
        mydf[j] = temp$df
        myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
        myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
        
        Hmat = temp$Hmat
        betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
        ## Use likelihood in obj instead; 
        etaAll = Xall%*%betaMat + Smat%*% (Z_sups - Xall%*%betaMat)
        res[j] = -logLikelihood(Yall,etaAll,exampleIdx)

        
      }else {
        mydf[j] = tr(Smat)
        myRSS[j] = sum((Z_sups-Smat%*%Z_sups)^2)
        myBeta[j,] = repmat(0,p)
        res[j] = - myRSS[j]
      } #end if(sum(lam1)>1e-3)
    }else {
      temp = getHmat(Xall, lam1, matrix(0,n,n), W_0)
      mydf[j] = temp$df
      myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
      myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
      
      ## res
      Hmat = temp$Hmat
      betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
      ## Use likelihood in obj instead; 
      etaAll = Xall%*%betaMat + Smat%*% (Z_sups - Xall%*%betaMat)
      res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
      
    }# end if(sum(lam2)>1e-3)
    myfU[j,] = Smat %*% (Z_sups - Xall%*%myBeta[j,])
  } # end for (j)
  # BIC = n*log(myRSS)+log(n)*mydf
  #BIC = (res)+log(n)*mydf
  BIC = (res)+log(n)*mydf
  #print('now print sum(res)')
  #print(sum(res))
  #print('now print log(n) * mydf')
  #print(log(n)*mydf)
  return(list(BIC = BIC, df = mydf, RSS = myRSS, lamhatALL = lamhatALL,myfU = myfU, 
              myBeta = myBeta, W_sups_std = W_sups_std, Z_sups_all = Z_sups_all))
}  # end function


naiveBICFULL_latest<-function(tauall_1,tauall_2, Xall, DDall,Uall,exampleIdx,Yall) {
  n = NROW(Xall)
  p = NCOL(Xall)
  q = NCOL(Uall)
  mydf = rep(-1000,length(tauall_1)*length(tauall_2))
  mydf_refit = rep(-1000,length(tauall_1)*length(tauall_2))
  myRSS = rep(-10000, length(tauall_1)*length(tauall_2))
  myRSS_refit = rep(-10000, length(tauall_1)*length(tauall_2))
  res = rep(-10000, length(tauall_1)*length(tauall_2))
  res_refit = rep(-10000, length(tauall_1)*length(tauall_2))
  # myBeta: 5 * size(tau) 
  myBeta = matrix(rep(-1000,length(tauall_1)*length(tauall_2)*p),ncol = p) 
  myBeta_refit = matrix(rep(0,length(tauall_1)*length(tauall_2)*p),ncol = p)
  # myfU: N * size(tau)
  myfU = matrix(rep(-1000,n * length(tauall_1)*length(tauall_2) ),ncol = n) 
  myfU_refit = matrix(rep(-1000,n * length(tauall_1)*length(tauall_2) ),ncol = n)
  
  lamhatALL = NULL
  lamhatALL_refit = NULL
  W_sups_std = NULL
  W_sups_std_refit = NULL
  Z_sups_all = NULL
  Z_sups_all_refit = NULL
  
  for (j in 1:(length(tauall_1)*length(tauall_2))) {
  
  #for (j in 1:(6*length(tauall_2))) {
    
    tau_1 = tauall_1[(j-1)%/%length(tauall_2)+1]
    tau_2 = tauall_2[(j-1)%%length(tauall_2)+1]
    
    result_MCD = MCD(Xall, DDall, tau_1,tau_2,exampleIdx,Yall,Uall)
    
    lamhat = result_MCD$lamhat
    Z_sups = result_MCD$Z_sups
    
    #### modified on 2023/09/21, instead of using Z_sups in refitting, we recalulcate the Z_sups_refit using Xall;
    #Z_sups_refit = result_MCD$Z_sups_refit
    Z_sups_refit = Z_sups
    
    W_0 = result_MCD$W_0
    
    lamhatALL = rbind(lamhatALL,lamhat)
    W_sups_std = cbind(W_sups_std,W_0)
    Z_sups_all = cbind(Z_sups_all,Z_sups)
    
    lam1 = lamhat[1:p]
    lam2 = lamhat[p+(1:q)]
    
    if(sum(lam2)>1e-3) {
      Smat = getLCSmat(DDall, lam2, W_0)
      if(sum(lam1)>1e-3) {
        temp = getHmat(Xall, lam1, Smat, W_0)
        mydf[j] = temp$df
        myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
        myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
        
        Hmat = temp$Hmat
        betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
        ## Use likelihood in obj instead; 
        etaAll = Xall%*%betaMat + Smat%*% (Z_sups - Xall%*%betaMat)
        res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
        
        ############ refitting #############
        lam1 = lamhat[1:p]
        lam2 = lamhat[p+(1:q)]
        Xind=which(lam1>1e-2)
        Uind=which(lam2>1e-2)
        
        x_ind = Xall[,Xind]
        dd_ind = DDall[Uind]
        u_ind = Uall[,Uind] 
        
        lam2_U = lam2[Uind]
        tmpLam = repmat(sum(lam2)/length(lam2_U),1,length(lam2_U))
        
        local_scoring_refit_res = local_scoring_refit(Xall[,Xind],DDall[Uind],Yall,Uall[,Uind],exampleIdx,c(lam1[Xind],tmpLam))
        Z_sups_refit = local_scoring_refit_res$Z_sups
        mydf_refit[j] = local_scoring_refit_res$mydf
        myBeta_refit[j,Xind] = local_scoring_refit_res$myBeta
        
        betaMat_refit = rep(0,p)
        betaMat_refit[Xind] = local_scoring_refit_res$myBeta
        
        ## Use likelihood in obj instead; 
        etaAll_refit = Xall%*%betaMat_refit + Smat%*% (Z_sups_refit - Xall%*%betaMat_refit)
        res_refit[j] = -logLikelihood(Yall,etaAll_refit,exampleIdx)
        
        
        
      }else {
        mydf[j] = tr(Smat)
        myRSS[j] = sum((Z_sups-Smat%*%Z_sups)^2)
        myBeta[j,] = repmat(0,p)
        res[j] = - myRSS[j]
        
        mydf_refit[j] = tr(Smat)
        myRSS_refit[j] = sum((Z_sups_refit-Smat%*%Z_sups_refit)^2)
        myBeta_refit[j,] = repmat(0,p)
        res_refit[j] = - myRSS[j]
      } #end if(sum(lam1)>1e-3)
    }else {
      temp = getHmat(Xall, lam1, matrix(0,n,n), W_0)
      mydf[j] = temp$df
      myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
      myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
      
      ## res
      Hmat = temp$Hmat
      betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
      ## Use likelihood in obj instead; 
      etaAll = Xall%*%betaMat + Smat%*% (Z_sups - Xall%*%betaMat)
      res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
      
      ############ refitting #############
      ############ refitting #############
      lam1 = lamhat[1:p]
      lam2 = lamhat[p+(1:q)]
      Xind=which(lam1>1e-2)
      Uind=which(lam2>1e-2)
      
      x_ind = Xall[,Xind]
      dd_ind = DDall[Uind]
      u_ind = Uall[,Uind] 
      
      lam2_U = lam2[Uind]
      tmpLam = repmat(sum(lam2)/length(lam2_U),1,length(lam2_U))
      
      local_scoring_refit_res = local_scoring_refit(Xall[,Xind],DDall[Uind],Yall,Uall[,Uind],exampleIdx,c(lam1[Xind],tmpLam))
      Z_sups_refit = local_scoring_refit_res$Z_sups
      mydf_refit[j] = local_scoring_refit_res$mydf
      myBeta_refit[j,Xind] = local_scoring_refit_res$myBeta
      
      betaMat_refit = rep(0,p)
      betaMat_refit[Xind] = local_scoring_refit_res$myBeta
      
      ## Use likelihood in obj instead; 
      etaAll_refit = Xall%*%betaMat_refit + Smat%*% (Z_sups_refit - Xall%*%betaMat_refit)
      res_refit[j] = -logLikelihood(Yall,etaAll_refit,exampleIdx)
      
      
    }# end if(sum(lam2)>1e-3)
    myfU[j,] = Smat %*% (Z_sups - Xall%*%myBeta[j,])
    myfU_refit[j,] = Smat %*% (Z_sups_refit - Xall%*%myBeta_refit[j,])
    
  } # end for (j)
  # BIC = n*log(myRSS)+log(n)*mydf
  BIC = (res)+log(n)*mydf
  BIC_refit = res_refit + log(n) * mydf_refit
  #BIC = (res)+mydf
  #BIC_refit = res_refit + mydf_refit
  return(list(BIC = BIC, df = mydf, RSS = myRSS, lamhatALL = lamhatALL,myfU = myfU, 
              myBeta = myBeta, W_sups_std = W_sups_std, Z_sups_all = Z_sups_all,BIC_refit = BIC_refit))
}  # end function


naiveBICFULL_postRefitting<-function(tauall_1,tauall_2, Xall, DDall,Uall,exampleIdx,Yall) {
  n = length(Xall)
  p = NCOL(Xall)
  q = NCOL(Uall)
  mydf = rep(-1000,length(tauall_1)*length(tauall_2))
  mydf_refit = rep(-1000,length(tauall_1)*length(tauall_2))
  myRSS = rep(-10000, length(tauall_1)*length(tauall_2))
  myRSS_refit = rep(-10000, length(tauall_1)*length(tauall_2))
  res = rep(-10000, length(tauall_1)*length(tauall_2))
  res_refit = rep(-10000, length(tauall_1)*length(tauall_2))
  # myBeta: 5 * size(tau) 
  myBeta = matrix(rep(-1000,length(tauall_1)*length(tauall_2)*p),ncol = p) 
  myBeta_refit = matrix(rep(-1000,length(tauall_1)*length(tauall_2)*p),ncol = p)
  # myfU: N * size(tau)
  myfU = matrix(rep(-1000,n * length(tauall_1)*length(tauall_2) ),ncol = n) 
  myfU_refit = matrix(rep(-1000,n * length(tauall_1)*length(tauall_2) ),ncol = n)
  
  lamhatALL = NULL
  lamhatALL_refit = NULL
  W_sups_std = NULL
  W_sups_std_refit = NULL
  Z_sups_all = NULL
  Z_sups_all_refit = NULL
  
  for (j in 1:(length(tauall_1)*length(tauall_2))) {
    
    tau_1 = tauall_1[(j-1)%/%length(tauall_2)+1]
    tau_2 = tauall_2[(j-1)%%length(tauall_2)+1]
    
    result_MCD = MCD(Xall, DDall, tau_1,tau_2,exampleIdx,Yall,Uall)
    
    lamhat = result_MCD$lamhat
    Z_sups = result_MCD$Z_sups
    W_0 = result_MCD$W_0

    lam1 = lamhat[1:p]
    lam2 = lamhat[p+(1:q)] 
    Xind=which(lam1>1e-2)
    Uind=which(lam2>1e-2)
    
    x_ind = Xall[,Xind]
    dd_ind = DDall[Uind]
    u_ind = Uall[,Uind]
    
    
    #result_MCD = MCD_refit(x_ind, dd_ind, tau_1,tau_2,exampleIdx,Yall,u_ind)
    result_MCD = MCD_refit(x_ind, dd_ind, tau_1,tau_2,exampleIdx,Yall,u_ind,Xall, Uall,DDall,Xind,Uind,lamhat,Z_sups,W_0)
    
    #lamhat_refit_temp = result_MCD$lamhat
    lamhat_refit = result_MCD$lamhat
    
    # idx_old = which(lamhat < 1e-2)
    # idx_new = which(lamhat >= 1e-2)
    # lamhat_refit = rep(-100,length(lamhat))
    # lamhat_refit[idx_old] = lamhat[idx_old]
    # lamhat_refit[idx_new] = lamhat_refit_temp
    
    Z_sups_refit = result_MCD$Z_sups
    W_0_refit = result_MCD$W_0
    
    lamhatALL = rbind(lamhatALL,lamhat)
    lamhatALL_refit = rbind(lamhatALL_refit,lamhat_refit)
    W_sups_std = cbind(W_sups_std,W_0)
    W_sups_std_refit = cbind(W_sups_std_refit,W_0_refit)
    Z_sups_all = cbind(Z_sups_all,Z_sups)
    Z_sups_all_refit = cbind(Z_sups_all_refit,Z_sups_refit)
    
    
    lam1 = lamhat[1:p]
    lam2 = lamhat[p+(1:q)]

    lam1_refit = lamhat_refit[1:p]
    lam2_refit = lamhat_refit[p+(1:q)]
    
    if(sum(lam2)>1e-3) {
      Smat = getLCSmat(DDall, lam2, W_0)
      if(sum(lam1)>1e-3) {
        temp = getHmat(Xall, lam1, Smat, W_0)
        mydf[j] = temp$df
        myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
        myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
        
        Hmat = temp$Hmat
        betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
        
        ################ Refitting ##################
        
        
        ## Use likelihood in obj instead; 
        etaAll = Xall%*%betaMat + Smat%*% (Z_sups - Xall%*%betaMat)
        res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
        
        
      }else {
        mydf[j] = tr(Smat)
        myRSS[j] = sum((Z_sups-Smat%*%Z_sups)^2)
        myBeta[j,] = repmat(0,p)
        res[j] = - myRSS[j]
      } #end if(sum(lam1)>1e-3)
    }else {
      temp = getHmat(Xall, lam1, matrix(0,n,n), W_0)
      mydf[j] = temp$df
      myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
      myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
      
      ## res
      Hmat = temp$Hmat
      betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
      ## Use likelihood in obj instead; 
      etaAll = Xall%*%betaMat + Smat%*% (Z_sups - Xall%*%betaMat)
      res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
      
    }# end if(sum(lam2)>1e-3)
    myfU[j,] = Smat %*% (Z_sups - Xall%*%myBeta[j,])
    
    ############################# post tuning ############################
    
    if(sum(lam2_refit)>1e-3) {
      Smat_refit = getLCSmat(DDall, lam2_refit, W_0_refit)
      if(sum(lam1_refit)>1e-3) {
        temp = getHmat_refit(Xall, lam1_refit, Smat_refit, W_0_refit)
        mydf_refit[j] = temp$df
        myRSS[j] = sum((temp$Hmat%*%Z_sups_refit)^2)
        #myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
        
        Hmat = temp$Hmat
        betaMat = temp$Amat %*% temp$Bmat %*% Z_sups_refit
        ## Use likelihood in obj instead; 
        etaAll = Xall%*%betaMat + Smat_refit%*% (Z_sups - Xall%*%betaMat)
        res_refit[j] = -logLikelihood(Yall,etaAll,exampleIdx)
        
        
      }else {
        mydf[j] = tr(Smat_refit)
        myRSS[j] = sum((Z_sups_refit-Smat_refit%*%Z_sups_refit)^2)
        #myBeta[j,] = repmat(0,p)
        res_refit[j] = - myRSS[j]
      } #end if(sum(lam1)>1e-3)
    }else {
      temp = getHmat_refit(Xall, lam1_refit, matrix(0,n,n), W_0_refit)
      mydf_refit[j] = temp$df
      myRSS[j] = sum((temp$Hmat%*%Z_sups_refit)^2)
      #myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
      ## res
      Hmat = temp$Hmat
      betaMat = temp$Amat %*% temp$Bmat %*% Z_sups_refit
      ## Use likelihood in obj instead; 
      etaAll = Xall%*%betaMat + Smat%*% (Z_sups_refit - Xall%*%betaMat)
      res_refit[j] = -logLikelihood(Yall,etaAll,exampleIdx)
      
    }# end if(sum(lam2)>1e-3)
    
    
    myfU_refit[j,] = Smat_refit %*% (Z_sups_refit - Xall%*%myBeta[j,])
    
    
  } # end for (j)
  # BIC = n*log(myRSS)+log(n)*mydf
  BIC = (res)+log(n)*mydf
  BIC_refit = (res_refit) + log(n) * mydf_refit

  return(list(BIC = BIC, df = mydf, RSS = myRSS, lamhatALL = lamhatALL,myfU = myfU, 
              myBeta = myBeta, W_sups_std = W_sups_std, Z_sups_all = Z_sups_all))
}  # end function




postRefitting <- function (tau_1,tau_2,exampleIdx,Yall,DDall,lamhat,Xall,Uall)
{
  # post refitting is performed by taking all candiate models and perform weighted least squared with smoothing without penalties. 
  
  p = NCOL(Xall)
  q = NCOL(Uall)
  n = length(Xall)
  
  lam1 = lamhat[1:p]
  lam2 = lamhat[p+(1:q)] 
  Xind=which(lam1>0)
  Uind=which(lam2>0)

  x_ind = Xall[,Xind]
  dd_ind = DDall[Uind]
  u_ind = Uall[,Uind]
  
  result_MCD = MCD_refit(x_ind, dd_ind, tau_1,tau_2,exampleIdx,Yall,u_ind) # NEED TO CHANGE 
  #a=LLMCD(Xall[,Xind], Uall[,Uind], tau2,exampleIdx,Yall,W_sups_std_min,Z_sups_all_min) # NEED TO CHANGE
  
  lamhat_refit = result_MCD$lamhat
  # Z_sups = result_MCD$Z_sups
  # W_0 = result_MCD$W_0
  # 
  # lamhatALL = rbind(lamhatALL,lamhat)
  # W_sups_std = cbind(W_sups_std,W_0)
  # Z_sups_all = cbind(Z_sups_all,Z_sups)
  # 
  # 
  # 
  # Smat = getLCSmat(DDall[Uind], lam2[Uind], W_sups_std[,j])
  # temp = getHmat(Xall[,Xind], lam1[Xind], Smat, W_sups_std[,j])
  # mydf[j] = temp$df
  # #myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
  # #myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
  # 
  # Hmat = temp$Hmat
  # Z_sups = Z_sups_all[,j]
  # betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
  ## Use likelihood in obj instead; 
  #etaAll = Xall[,Xind]%*%betaMat + Smat%*% (Z_sups - Xall[,Xind]%*%betaMat)

  #res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
  
#}
#pBIC = (res)+log(n)*mydf
  
  return(list(lamhat_refit = lamhat_refit))
}


postRefitting_tuning <- function (tauall_1,tauall_2,exampleIdx,Yall,DDall,SolutionPath,Xall,Uall,candMODEL)
{
  # post refitting is performed by taking all candiate models and perform weighted least squared with smoothing without penalties. 
  
  p0 = NCOL(Xall)
  q0 = NCOL(Uall)
  n = length(Xall)
  pBIC = NULL
  pRes = NULL
  pDF=NULL
  pRSS=NULL
  for (k in 1:nrow(candMODEL)) {
    model=candMODEL[k,]
    lam1_j=model[1:p0]
    lam2_j=model[p0+(1:q0)]
    Xind=which(lam1_j>0)
    Uind=which(lam2_j>0)
    
    x_ind = Xall[,Xind]
    dd_ind = DDall[Uind]
    u_ind = Uall[,Uind]
    
    mydf = rep(-1000,length(tauall_1)*length(tauall_2))
    myRSS = rep(-10000, length(tauall_1)*length(tauall_2))
    res = rep(-10000, length(tauall_1)*length(tauall_2))
    # myBeta: 5 * size(tau) 
    #myBeta = matrix(rep(-1000,length(tauall_1)*length(tauall_2)*p),ncol = p0) 
    # myfU: N * size(tau)
    myfU = matrix(rep(-1000,n * length(tauall_1)*length(tauall_2) ),ncol = n) 
    
    lamhatALL = NULL
    W_sups_std = NULL
    Z_sups_all = NULL
    
    
    for (j in 1:(length(tauall_1)*length(tauall_2))) {
      
      tau_1 = tauall_1[(j-1)%/%length(tauall_2)+1]
      tau_2 = tauall_2[(j-1)%%length(tauall_2)+1]
      
      result_MCD = MCD_refit(x_ind, dd_ind, tau_1,tau_2,exampleIdx,Yall,u_ind) # NEED TO CHANGE 

      lamhat = result_MCD$lamhat
      Z_sups = result_MCD$Z_sups
      W_0 = result_MCD$W_0
      
      lamhatALL = rbind(lamhatALL,lamhat)
      W_sups_std = cbind(W_sups_std,W_0)
      Z_sups_all = cbind(Z_sups_all,Z_sups)
      
      p = length(Xind)
      q = length(Uind)
      
      lam1 = lamhat[1:p]
      lam2 = lamhat[p+(1:q)]
      
      if(sum(lam2)>1e-3) {
        Smat = getLCSmat(dd_ind, lam2, W_0)
        if(sum(lam1)>1e-3) {
          temp = getHmat_refit(x_ind, lam1, Smat, W_0)
          mydf[j] = temp$df
          myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
          #myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
          
          Hmat = temp$Hmat
          betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
          ## Use likelihood in obj instead; 
          etaAll = x_ind%*%betaMat + Smat%*% (Z_sups - x_ind%*%betaMat)
          res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
          
          
        }else {
          mydf[j] = tr(Smat)
          myRSS[j] = sum((Z_sups-Smat%*%Z_sups)^2)
          #myBeta[j,] = repmat(0,p)
          res[j] = - myRSS[j]
        } #end if(sum(lam1)>1e-3)
      }else {
        temp = getHmat_refit(x_ind, lam1, matrix(0,n,n), W_0)
        mydf[j] = temp$df
        myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
        #myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
        ## res
        Hmat = temp$Hmat
        betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
        ## Use likelihood in obj instead; 
        etaAll = x_ind%*%betaMat + Smat%*% (Z_sups - x_ind%*%betaMat)
        res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
        
      }# end if(sum(lam2)>1e-3)
      #myfU[j,] = Smat %*% (Z_sups - x_ind%*%myBeta[j,])
    }
    
    BIC = (res) + log(n)*mydf
    ind=which.min(BIC)
    pRes[k]=res[ind]
    pDF[k]=mydf[ind]
  }
    
  pBIC = pRes + log(n)*pDF
    
  # Z_sups = result_MCD$Z_sups
  # W_0 = result_MCD$W_0
  # 
  # lamhatALL = rbind(lamhatALL,lamhat)
  # W_sups_std = cbind(W_sups_std,W_0)
  # Z_sups_all = cbind(Z_sups_all,Z_sups)
  # 
  # 
  # 
  # Smat = getLCSmat(DDall[Uind], lam2[Uind], W_sups_std[,j])
  # temp = getHmat(Xall[,Xind], lam1[Xind], Smat, W_sups_std[,j])
  # mydf[j] = temp$df
  # #myRSS[j] = sum((temp$Hmat%*%Z_sups)^2)
  # #myBeta[j,] = temp$Amat %*% temp$Bmat %*% Z_sups
  # 
  # Hmat = temp$Hmat
  # Z_sups = Z_sups_all[,j]
  # betaMat = temp$Amat %*% temp$Bmat %*% Z_sups
  ## Use likelihood in obj instead; 
  #etaAll = Xall[,Xind]%*%betaMat + Smat%*% (Z_sups - Xall[,Xind]%*%betaMat)
  
  #res[j] = -logLikelihood(Yall,etaAll,exampleIdx)
  
  #}
  #pBIC = (res)+log(n)*mydf
  
  return(list(pBIC = pBIC))
}

profileBICFULL<-function(Xall, Uall, candMODEL, tau2allmain,exampleIdx, Yall,W_sups_std_min,Z_sups_all_min) {
  p=NCOL(Xall)
  q=NCOL(Uall)
  n=length(Xall)
  
  pBIC=NULL
  bic2ALL=NULL
  pDF=NULL
  pRSS=NULL
  pLambdaZ = list()
  
  for (j in 1:nrow(candMODEL)) {
    model=candMODEL[j,]
    lam1=model[1:p]
    lam2=model[p+(1:q)]
    Xind=which(lam1)
    Uind=which(lam2)
    
    if(length(Uind)<0.5){
      
      
      pRSS[j]=deviance(lm(Z_sups~Xall[,Xind]))
      pDF[j]=length(Xind)+1
      pLambdaZ[j] = NULL
        
    }else {
      mydf2=rep(-1000,length(tau2allmain))
      myRSS2=rep(-10000, length(tau2allmain))
      myLambdaZ = c()
      for (k in 1:length(tau2allmain)){
        tau2=tau2allmain[k]
        a=LLMCD(Xall[,Xind], Uall[,Uind], tau2,exampleIdx,Yall,W_sups_std_min,Z_sups_all_min) # NEED TO CHANGE
        Z_sups = a$Z_sups
        W_sups_std = a$W_0
        lambdaZ = a$Flambda
        myLambdaZ = rbind(myLambdaZ,lambdaZ)
        # x, Uall, Z_sups, lambdaZ,W_sups_std,exampleIdx,Yall
        temp=LLobjfun(Xall[,Xind], Uall[,Uind], Z_sups,lambdaZ,W_sups_std,exampleIdx,Yall)
        mydf2[k]=temp$df
        myRSS2[k]=temp$RSS
      } # end for (k in 1:length(tau2all))
      BIC2=n*log(myRSS2)+log(n)*mydf2
      bic2ALL=rbind(bic2ALL,BIC2)
      ind=which.min(BIC2)
      pRSS[j]=myRSS2[ind]
      pDF[j]=mydf2[ind]
      print(pLambdaZ)
      pLambdaZ[j] = list(myLambdaZ[ind,])
    } 
    
    #myfU[j,] = Smat %*% (Z_sups - Xall%*%myBeta[j,])
    
    # end if(length(Zind)<0.5){
  } # end for (j in 1:nrow(candMODEL))
  
  pBIC=n*log(pRSS)+log(n)*pDF
  
  return(list(pBIC=pBIC, bic2ALL=bic2ALL,pLambdaZ = pLambdaZ))
} # end function profileBIC

#############################################################

naiveBIC<-function(tauall, Xall, DDall, Yall,Zall) {
  n=length(Yall)
  p=NCOL(Xall)
  q=NCOL(Zall)
  mydf=NULL # rep(-1000,length(tauall))
  myRSS=NULL # rep(-10000, length(tauall))
  lamhatALL=NULL
  
  for (j in 1:length(tauall)) {
    tau=tauall[j]
    
    lamhat=MCD(Xall, DDall, Yall, tau)
    lamhatALL=rbind(lamhatALL,lamhat)
    
    lam1=lamhat[1:p]
    lam2=lamhat[p+(1:q)]
    
    if(sum(lam2)>1e-3) {
      Smat=getLCSmat(DDall, lam2)
      if(sum(lam1)>1e-3) {
        temp=getHmat(Xall, lam1, Smat)
        #  mydf[j]=temp$df
        #  myRSS[j]=sum((temp$Hmat%*%Yall)^2)
        mydf=c(mydf,temp$df)
        myRSS=c(myRSS,sum((temp$Hmat%*%Yall)^2))
      }else {
        #  mydf[j]=tr(Smat)
        #  myRSS[j]=sum((Yall-Smat%*%Yall)^2)
        mydf=c(mydf,tr(Smat))
        myRSS=c(myRSS,sum((Yall-Smat%*%Yall)^2))
      } #end ifif(sum(lam1)>1e-3)
    }else {
      temp=getHmat(Xall, lam1, matrix(0,n,n))
      #  mydf[j]=temp$df
      #  myRSS[j]=sum((temp$Hmat%*%Yall)^2)
      mydf=c(mydf,temp$df)
      myRSS=c(myRSS,sum((temp$Hmat%*%Yall)^2))
    }# end if(sum(lam2)>1e-3)
    BIC=n*log(myRSS)+log(n)*mydf
    if(j>1.5) {
      if(BIC[j]-BIC[j-1]>1e-2) {
        break;
      }
      
    }
  } # end for (j)
  BIC=n*log(myRSS)+log(n)*mydf
  return(list(BIC=BIC, df=mydf, RSS=myRSS, lamhatALL=lamhatALL))
}  # end function


profileBIC<-function(Xall, Zall, Yall, candMODEL, tau2all) {
  p=NCOL(Xall)
  q=NCOL(Zall)
  n=length(Yall)
  pBIC=NULL
  bic2ALL=NULL
  
  pDF=NULL
  pRSS=NULL
  
  for (j in 1:nrow(candMODEL)) {
    model=candMODEL[j,]
    lam1=model[1:p]
    lam2=model[p+(1:q)]
    Xind=which(lam1)
    Zind=which(lam2)
    
    if(length(Zind)<0.5){
      #  pRSS[j]=deviance(lm(Yall~Xall[,Xind]))
      #  pDF[j]=length(Xind)+1
      pRSS=c(pRSS, deviance(lm(Yall~Xall[,Xind])))
      pDF=c(pDF, length(Xind)+1)
    }else {
      #mydf2=rep(-1000,length(tau2all))
      #myRSS2=rep(-10000, length(tau2all))
      mydf2=NULL
      myRSS2=NULL
      for (k in 1:length(tau2all)){
        tau2=tau2all[k]
        a=LLMCD(Xall[,Xind], Zall[,Zind], Yall, tau2)
        temp=LLobjfun(Xall[,Xind], Zall[,Zind], Yall, a)
        #  mydf2[k]=temp$df
        #  myRSS2[k]=temp$RSS
        mydf2=c(mydf2,temp$df)
        myRSS2=c(myRSS2,temp$RSS)
        BIC2=n*log(myRSS2)+log(n)*mydf2
        
        if(k>1.5) {
          if(BIC2[k]-BIC2[k-1]>1e-2) {
            break;
          }
        }
        
      } # end for (k in 1:length(tau2all))
      BIC2=n*log(myRSS2)+log(n)*mydf2
      
      bic2ALL=rbind(bic2ALL,c(BIC2, rep(999999999, length(tau2all)-length(BIC2))))
      ind=which.min(BIC2)
      
      # pRSS[j]=myRSS2[ind]
      # pDF[j]=mydf2[ind]
      
      pRSS=c(pRSS, myRSS2[ind])
      pDF=c(pDF,mydf2[ind])
      
    } # end if(length(Zind)<0.5){
    
    pBIC=n*log(pRSS)+log(n)*pDF
    
    if(j>1.5) {
      if(pBIC[j]-pBIC[j-1]>1e-2) {
        break;
      }
    }
    
  } # end for (j in 1:nrow(candMODEL))
  
  pBIC=n*log(pRSS)+log(n)*pDF
  
  return(list(pBIC=pBIC, bic2ALL=bic2ALL))
} # end function profileBIC