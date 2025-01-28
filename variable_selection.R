
source("main_2tau_deviance.R")
source("examples.R")
source("getVS.R")

variable_selection <-function(tauall_1,tauall_2,Xall,DDall,Yall,Uall,exampleIdx,tau2allmain) {
  
  if_variable_selection = 1
  p = NCOL(Xall)
  #wuNaive = naiveBICFULL(tauall_1,tauall_2, Xall, DDall,Uall,exampleIdx,Yall)
  wuNaive = naiveBICFULL_latest(tauall_1,tauall_2, Xall, DDall,Uall,exampleIdx,Yall)
  #print(paste("irep:", irep, "BIC done"))
  #W_sups_std = wuNaive$W_sups_std
  Beta_sups = wuNaive$myBeta
  Beta_sups_min = matrix(Beta_sups[which.min(wuNaive$BIC),])
  fU_sups = wuNaive$myfU
  fU_sups_min = matrix(fU_sups[which.min(wuNaive$BIC),])
  SolutionPath = wuNaive$lamhatALL[which.min(wuNaive$BIC),]
  W_sups_std = wuNaive$W_sups_std
  Z_sups_all = wuNaive$Z_sups_all
  W_sups_std_min = W_sups_std[,which.min(wuNaive$BIC)]
  Z_sups_all_min = Z_sups_all[,which.min(wuNaive$BIC)]
  lamhatALL = wuNaive$lamhatALL
  candMODELmain = unique(lamhatALL>1e-2)
  print('pBIC start....')
  j = which.min(wuNaive$BIC)
  tau_1 = tauall_1[(j-1)%/%length(tauall_2)+1]
  tau_2 = tauall_2[(j-1)%%length(tauall_2)+1]
  # Youhan Lu add refitting step on 2023.07.14. 
  # The refitting is essentially taking all candidate variables and do vs again. 
  # lamhat = SolutionPath
  # wuRefitting_lambta = postRefitting(tau_1,tau_2,exampleIdx,Yall,DDall,SolutionPath,Xall,Uall)
  # wuRefitting_pBIC = postRefitting_tuning(tauall_1,tauall_2,exampleIdx,Yall,DDall,SolutionPath,Xall,Uall,candMODELmain)
  
  #wuRefitting = naiveBICFULL(tauall_1,tauall_2, Xall_refit, DDall_refit,Uall_refit,exampleIdx,Yall)
  #wuProfile = profileBICFULL(Xall, Uall, candMODELmain, tau2allmain,exampleIdx,Yall,W_sups_std_min,Z_sups_all_min)
  #result=list(wuNaive = wuNaive,tauall_1 = tauall_1, tauall_2 = tauall_2,
  #            SolutionPath = SolutionPath,wuProfile = wuProfile)
  # idx_old = which(lamhat < 1e-2)
  # idx_new = which(lamhat >= 1e-2)
  # wuRefitting = rep(-100,length(lamhat))
  # wuRefitting[idx_old] = lamhat[idx_old]
  # wuRefitting[idx_new] = wuRefitting_lambta$lamhat_refit
  # 
  # wuRefitting_tuning = candMODELmain[which.min(wuRefitting_pBIC$pBIC),]
  
  result=list(wuNaive = wuNaive,tauall_1 = tauall_1, tauall_2 = tauall_2,SolutionPath = SolutionPath)
  
  print('VS done')
  return(result)
  
}
  # if (max(W_sups_std - W_sups_all[,nloop] )< epsilon){
  #   break
  # }
  # while (if_variable_selection) {
    
    # eta_sups = Xall%*%Beta_sups_min + fU_sups_min
    # miu_sups = functionGInverse(eta_sups,exampleIdx)
    # W_sups = 1/(functionVar(eta_sups,exampleIdx)) * (functionGInversePrime(eta_sups,exampleIdx))^2
    # W_sups_std = W_sups/sum(W_sups)
    # Z_sups = eta_sups + (Yall - miu_sups) * functionGPrime(miu_sups,exampleIdx)
    # epsilon = 1e-3
    # wuNaive = naiveBICFULL(tauall_1,tauall_2, Xall, DDall,Uall,exampleIdx,Yall)
    
    # Beta_sups = wuNaive$myBeta
    # Beta_sups_min = matrix(Beta_sups[which.min(wuNaive$BIC),])
    # fU_sups = wuNaive$myfU
    # fU_sups_min = matrix(fU_sups[which.min(wuNaive$BIC),])
    # SolutionPath = wuNaive$lamhatALL[which.min(wuNaive$BIC),]
    # 
    # nloop = nloop + 1
    # W_sups_all = cbind(W_sups_all,W_sups_std)
    # SolutionPath_all = rbind(SolutionPath_all,SolutionPath)
    
    # if (max(W_sups_std - W_sups_all[,nloop] )< epsilon){
    #   break
    # }
    
  #   if (nloop >= 15){
  #     break
  #   }
  # }
