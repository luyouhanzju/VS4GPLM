source("local_scoring.R")
source("main_2tau_deviance.R")

MCD<-function(Xall, DDall, tau_1, tau_2,exampleIdx,Yall,Uall) {
  
  eps1=1e-4 # convergence criterion
  
  eps0=2*eps1
  # tau=10
  p=NCOL(Xall)
  q=length(DDall)
  d=p+q
  gammaCUR = rep(tau_1/p,p)
  #gammaCUR=c(0,tau1/3, 0, tau1/3, 0, 0, 0, tau1/3,rep(0,p-8)) # see formula 12
  lambdaCUR=rep(tau_2/q, q)
  #lambdaCUR=c(tau2/2, rep(0,(q-3)),tau2/2,0)
  
  
  #objCUR=objfun(X, DD, y, lambdaCUR)
  mygrid=c(0,0.2,0.5,0.8, 1)
  # if (tau2 < 3){
  #   mygrid=c(0,0.25, 1)
  # }
  
  #mygrid=c(0, 0.15,0.3,0.45,0.6,0.75,0.9,1)
  #mygrid=c(0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
  
  test=1
  
  nloop=0
  gammaLambdaALL=NULL
  objVall=NULL
  W_0_all = NULL
  Z_sups_all = NULL
  Z_sups_refit_all = NULL
  
  while(test) {
    gammaOLD = gammaCUR
    lambdaOLD = lambdaCUR
    
    for (i in 1:p){
      gamma1 = gammaCUR
      gamma1[i] = 0
      # if (any(is.na(gamma1))){
      #   print(paste("irep:", irep, "gamma1 NA"))
      # }
      if (sum(gamma1)>eps0){
        gamma1 = gamma1/sum(gamma1)
        gamma2 = rep(0,p)
        gamma2[i] = 1
        objobj2 = rep(100000,length(mygrid))
        for (m in 1:length(mygrid)){
          gam_gamma = mygrid[m]
          tempGAMMA=(gam_gamma*gamma1+(1-gam_gamma)*gamma2)*tau_1
          
          local_scoring_res = local_scoring(Xall,DDall,Yall,Uall,exampleIdx,c(tempGAMMA,lambdaCUR))
          
          Z_sups = local_scoring_res$Z_sups
          W_0 = local_scoring_res$W_0
          
          objobj2[m]=objfun(Xall, DDall, Z_sups , c(tempGAMMA,lambdaCUR),W_0,exampleIdx, Yall)
          
        }
        gam_gamma = mygrid[which.min(objobj2)]
        gammaCUR=(gam_gamma*gamma1+(1-gam_gamma)*gamma2)*tau_1
      } # end for over j
      
    }
    
    for(j in 1:q) {
      
      lam1 = lambdaCUR
      lam1[j] = 0
      # if (any(is.na(lam1))){
      #   print(paste("irep:", irep, "lam1 NA"))
      # }
      if(sum(lam1)>eps0) {
        lam1 = lam1/sum(lam1)
        
        lam2 = rep(0, q)
        lam2[j] = 1
        objobj = rep(100000, length(mygrid))
        for (k in 1:length(mygrid)){
          gam_lam = mygrid[k]
          tempLAM = (gam_lam*lam1+(1-gam_lam)*lam2)*tau_2
          
          local_scoring_res = local_scoring(Xall,DDall,Yall,Uall,exampleIdx,c(gammaCUR,tempLAM))
          Z_sups = local_scoring_res$Z_sups
          W_0 = local_scoring_res$W_0
          
          objobj[k] = objfun(Xall, DDall, Z_sups, c(gammaCUR,tempLAM),W_0, exampleIdx,Yall)
        }  # end for over k
        gam_lam = mygrid[which.min(objobj)]
        lambdaCUR = (gam_lam*lam1+(1-gam_lam)*lam2)*tau_2
        #print(lambdaCUR)
      } # end if
    }
    
    
    gammaLambdaALL=rbind(gammaLambdaALL,c(gammaCUR,lambdaCUR))
    
    local_scoring_res = local_scoring(Xall,DDall,Yall,Uall,exampleIdx,c(gammaCUR,lambdaCUR))
    
    W_0 = local_scoring_res$W_0
    Z_sups = local_scoring_res$Z_sups
    Z_sups_refit = local_scoring_res$Z_sups_refit
    
    objVall = c(objVall, objfun(Xall, DDall, Z_sups, c(gammaCUR,lambdaCUR), W_0, exampleIdx, Yall))
    nloop = nloop+1
    
    W_0_all = cbind(W_0_all,W_0)
    Z_sups_all = cbind(Z_sups_all,Z_sups)
    Z_sups_refit_all = cbind(Z_sups_refit_all,Z_sups_refit)
    
    if(nloop>8){
      Flambda = gammaLambdaALL[which.min(objVall),]
      Z_sups = Z_sups_all[,which.min(objVall)] # 100 * 1 matrix
      Z_sups_refit = Z_sups_refit_all[,which.min(objVall)]
      W_0 = W_0_all[,which.min(objVall)] # 100 * 1 matrix
      #print(c(gammaCUR,lambdaCUR))
      #print(min(objVall))
      #print(which.min(objVall))
      print('Break while loop after more than 8 iterations')
      break
    }
    
    
    #  print(lambdaCUR)
    if(max(c(abs(lambdaOLD-lambdaCUR),abs(gammaOLD-gammaCUR)))<eps1) {
      #  test=0
      Flambda=c(gammaCUR,lambdaCUR)
      
      #print(objfun(X, DD, y, c(gammaCUR,lambdaCUR)))
      #print ('Break because convergence')
      #print(nloop)
      #print(c(gammaCUR,lambdaCUR))
      #print(c(abs(lambdaOLD-lambdaCUR),abs(gammaOLD-gammaCUR)))
      break
      
    } # end if 
    
  } # end while
  #print(Flambda)
  
  return(list(Flambda = Flambda, lamhat = Flambda,Z_sups = Z_sups, W_0 = W_0, Z_sups_refit = Z_sups_refit))
}# end of function



############### POST REFITTING MCD


MCD_refit<-function(x_ind, dd_ind, tau_1, tau_2,exampleIdx,Yall,u_ind,Xall, Uall,DDall,Xind,Uind,lamhat,Z_sups,W_0) {
  
  eps1=1e-4 # convergence criterion
  
  eps0=2*eps1
  p=NCOL(x_ind)
  q=length(dd_ind)
  d=p+q
  gammaCUR = rep(tau_1/p,p)
  lambdaCUR=rep(tau_2/q, q)

  
  mygrid=c(0,0.2,0.5,0.8, 1)
  
  test=1
  
  nloop=0
  gammaLambdaALL=NULL
  objVall=NULL
  W_0_all = NULL
  Z_sups_all = NULL
  
  
  while(test) {
    gammaOLD = gammaCUR
    lambdaOLD = lambdaCUR
    
    for (i in 1:p){
      gamma1 = gammaCUR
      gamma1[i] = 0

      if (sum(gamma1)>eps0){
        gamma1 = gamma1/sum(gamma1)
        gamma2 = rep(0,p)
        gamma2[i] = 1
        objobj2 = rep(100000,length(mygrid))
        for (m in 1:length(mygrid)){
          gam_gamma = mygrid[m]
          tempGAMMA=(gam_gamma*gamma1+(1-gam_gamma)*gamma2)*tau_1
          
          FullLambta = c(tempGAMMA,lambdaCUR)
          
          #local_scoring_res = local_scoring_refit(x_ind,dd_ind,Yall,u_ind,exampleIdx,FullLambta)
          
          FullLambta_refit = rep(0,length(lamhat))
          
          #idx_old = which(lamhat < 1e-2)
          idx_new = which(lamhat >= 1e-2)
          FullLambta_refit[idx_new] = FullLambta
          print(FullLambta_refit)
          #local_scoring_res = local_scoring_refit(Xall,DDall,Yall,Uall,exampleIdx,FullLambta_refit)
          #local_scoring_res = local_scoring(Xall,DDall,Yall,Uall,exampleIdx,FullLambta_refit)
          
          # Z_sups = local_scoring_res$Z_sups
          # W_0 = local_scoring_res$W_0
          
          #objobj2[m]=objfun_refit(x_ind, dd_ind, Z_sups , FullLambta,W_0,exampleIdx, Yall)
          #objobj2[m]=objfun(Xall, DDall, Z_sups , FullLambta_refit,W_0,exampleIdx, Yall)
          lam1_temp = FullLambta_refit[1:NCOL(Xall)]
          lam2_temp = FullLambta_refit[(NCOL(Xall)+1):(NCOL(Xall) + length(DDall))]
          objobj2[m] = objfun_refit(Xall[,which(lam1_temp >= 1e-2)], DDall[lam2_temp >= 1e-2], Z_sups , FullLambta_refit[FullLambta_refit>=1e-2],W_0,exampleIdx, Yall)
        }
        gam_gamma = mygrid[which.min(objobj2)]
        gammaCUR=(gam_gamma*gamma1+(1-gam_gamma)*gamma2)*tau_1
      } # end for over j
      
    }
    
    # for(j in 1:q) {
    #   
    #   lam1 = lambdaCUR
    #   lam1[j] = 0
    # 
    #   if(sum(lam1)>eps0) {
    #     lam1 = lam1/sum(lam1)
    #     
    #     lam2 = rep(0, q)
    #     lam2[j] = 1
    #     objobj = rep(100000, length(mygrid))
    #     for (k in 1:length(mygrid)){
    #       gam_lam = mygrid[k]
    #       tempLAM = (gam_lam*lam1+(1-gam_lam)*lam2)*tau_2
    #       FullLambta = c(gammaCUR,tempLAM)
    #       
    #       FullLambta_refit = rep(0,length(lamhat))
    #       idx_new = which(lamhat >= 1e-2)
    #       FullLambta_refit[idx_new] = FullLambta
    #       
    #       #local_scoring_res = local_scoring_refit(x_ind,dd_ind,Yall,u_ind,exampleIdx,FullLambta)
    #       #local_scoring_res = local_scoring_refit(Xall,DDall,Yall,Uall,exampleIdx,FullLambta_refit)
    #       #local_scoring_res = local_scoring(Xall,DDall,Yall,Uall,exampleIdx,FullLambta_refit)
    #       
    #       #Z_sups = local_scoring_res$Z_sups
    #       #W_0 = local_scoring_res$W_0
    #       
    #       #objobj[k] = objfun_refit(x_ind, dd_ind, Z_sups , FullLambta, W_0, exampleIdx, Yall)
    #       #objobj[k] = objfun(Xall, DDall, Z_sups , FullLambta_refit, W_0, exampleIdx, Yall)
    #       lam1_temp = FullLambta_refit[1:NCOL(Xall)]
    #       lam2_temp = FullLambta_refit[(NCOL(Xall)+1):(NCOL(Xall) + length(DDall))]
    #       objobj[k]=objfun_refit(Xall[,which(lam1_temp >= 1e-2)], DDall[lam2_temp >= 1e-2], Z_sups , FullLambta_refit[FullLambta_refit>=1e-2],W_0,exampleIdx, Yall)
    #       
    #     }  # end for over k
    #     gam_lam = mygrid[which.min(objobj)]
    #     lambdaCUR = (gam_lam*lam1+(1-gam_lam)*lam2)*tau_2
    #     #print(lambdaCUR)
    #   } # end if
    # }
    
    
    gammaLambdaALL=rbind(gammaLambdaALL,c(gammaCUR,lambdaCUR))
    
    FullLambta = c(gammaCUR,lambdaCUR)
    FullLambta_refit = rep(0,length(lamhat))
    
    #idx_old = which(lamhat < 1e-2)
    idx_new = which(lamhat >= 1e-2)
    FullLambta_refit[idx_new] = FullLambta
    
#    local_scoring_res = local_scoring_refit(x_ind,dd_ind,Yall,u_ind,exampleIdx,FullLambta)
#    local_scoring_res = local_scoring_refit(Xall,DDall,Yall,Uall,exampleIdx,FullLambta_refit)
    
    # local_scoring_res = local_scoring(Xall,DDall,Yall,Uall,exampleIdx,FullLambta_refit)
    # W_0 = local_scoring_res$W_0
    # Z_sups = local_scoring_res$Z_sups

    #objVall = c(objVall, objfun_refit(x_ind, dd_ind, Z_sups, FullLambta, W_0, exampleIdx, Yall))
    #objVall = c(objVall, objfun(Xall, DDall, Z_sups, FullLambta_refit, W_0, exampleIdx, Yall))
    lam1_temp = FullLambta_refit[1:NCOL(Xall)]
    lam2_temp = FullLambta_refit[(NCOL(Xall)+1):(NCOL(Xall) + length(DDall))]
    objVall = c(objVall, objfun_refit(Xall[,which(lam1_temp >= 1e-2)], DDall[lam2_temp >= 1e-2], Z_sups, FullLambta_refit[FullLambta_refit>=1e-2], W_0, exampleIdx, Yall))
    #objobj[k]=objfun_refit(Xall[,which(lam1_temp >= 1e-2)], DDall[lam2_temp >= 1e-2], Z_sups , FullLambta_refit[FullLambta_refit>=1e-2],W_0,exampleIdx, Yall)
    
    nloop = nloop+1
    
    W_0_all = cbind(W_0_all,W_0)
    Z_sups_all = cbind(Z_sups_all,Z_sups)
    
    if(nloop>30){
      Flambda = gammaLambdaALL[which.min(objVall),]
      FullLambta_refit = rep(0,length(lamhat))
      
      #idx_old = which(lamhat < 1e-2)
      idx_new = which(lamhat >= 1e-2)
      FullLambta_refit[idx_new] = Flambda
      
      Z_sups = Z_sups_all[,which.min(objVall)] # 100 * 1 matrix
      W_0 = W_0_all[,which.min(objVall)] # 100 * 1 matrix
      print('Break while loop after more than 30 iterations')
      break
    }
    
    
    #  print(lambdaCUR)
    if(max(c(abs(gammaOLD-gammaCUR)))<eps1) {
      #  test=0
      Flambda=c(gammaCUR,lambdaCUR)
      FullLambta_refit = rep(0,length(lamhat))
      
      #idx_old = which(lamhat < 1e-2)
      idx_new = which(lamhat >= 1e-2)
      FullLambta_refit[idx_new] = Flambda
      
      break
      
    } # end if 
    
  } # end while
  #print(Flambda)
  
  return(list(Flambda = FullLambta_refit, lamhat = FullLambta_refit,Z_sups = Z_sups, W_0 = W_0))
}# end of function


###############


###############
#
#   local linear smoothing for refitting
#
############### X, DD, Z_sups, tau1, tau2,W_0,exampleIdx,Yall

# a=LLMCD(Xall[,Xind], Uall[,Uind], Z_sups, tau2)
# temp=LLobjfun(Xall[,Xind], Uall[,Uind], Z_sups, a)

#        a=LLMCD(Xall[,Xind], Uall[,Uind], tau2,exampleIdx,Yall) # NEED TO CHANGE
#Z_sups = a$Z_sups
#W_sups_std = a$W_0
LLMCD<-function(X, U, tau,exampleIdx,Yall,W_sups_std_min,Z_sups_all_min) {
  # this version only do variable selection for nonparametric part
  #print(X)
  X = as.matrix(X)
  U = as.matrix(U)
  eps1=1e-4 # convergence criterion
  
  eps0=2*eps1
  #  p=NCOL(X)
  q=NCOL(U)
  #  d=p+q
  lambdaCUR=rep(tau/q, q)
  #lambdaCUR=c(tau/2, rep(0,(q-3)),tau/2,0)
  #  objCUR=objfun(X, DD, y, lambdaCUR)
  
  mygrid=c(0, 0.2,0.5,0.8, 1)
  
  test=1
  
  nloop=0
  lambdaALL=NULL
  objVall=NULL
  W_0_all = NULL
  Z_sups_all = NULL
  
  while(test) {
    
    lambdaOLD=lambdaCUR
    
    for(j in 1:q) {
      lam1=lambdaCUR
      lam1[j]=0
      # if (any(is.na(lam1))){
      #   print(paste("irep:", irep, "lam1 NA"))
      # }
      if(sum(lam1)>eps0) {
        lam1 = lam1/sum(lam1)
        
        lam2=rep(0, q)
        lam2[j]=1
        objobj=rep(100000, length(mygrid))
        for (k in 1:length(mygrid)){
          gam = mygrid[k]
          tempLAM = (gam*lam1+(1-gam)*lam2)*tau

          # local_scoring_res = local_scoring_LL(X,Yall,U,exampleIdx,tempLAM) # NEED TO CHANGE HERE
          # W_sups_std = local_scoring_res$W_0
          # Z_sups = local_scoring_res$Z_sups
          W_sups_std = W_sups_std_min
          Z_sups =  Z_sups_all_min
          objobj[k] = LLobjfun(X,U,Z_sups,tempLAM,W_sups_std,exampleIdx,Yall)$MinusLoglikelihood 
        }  # end for over k
        
        gam = mygrid[which.min(objobj)]
        lambdaCUR = (gam*lam1+(1-gam)*lam2)*tau
      } # end if
    } # end for over j
    
    
    lambdaALL = rbind(lambdaALL,lambdaCUR)
    # local_scoring_res = local_scoring_LL(X,Yall,U,exampleIdx,lambdaCUR) # NEED TO CHANGE HERE
    # W_0 = local_scoring_res$W_0
    # Z_sups = local_scoring_res$Z_sups
    W_0 = W_sups_std_min
    Z_sups = Z_sups_all_min
    objVall = c(objVall, LLobjfun(X,U,Z_sups,lambdaCUR,W_0,exampleIdx,Yall)$MinusLoglikelihood)
    nloop = nloop+1
    
    W_0_all = cbind(W_0_all,W_0)
    Z_sups_all = cbind(Z_sups_all,Z_sups)
    
    if(nloop>15){
      
      Flambda = lambdaALL[which.min(objVall),]
      Z_sups = Z_sups_all[,which.min(objVall)] # 100 * 1 matrix
      W_0 = W_0_all[,which.min(objVall)] # 100 * 1 matrix
      
      print('Break while loop after more than 15 iterations')
      break
    }
    
    
    #  print(lambdaCUR)
    if(max(abs(lambdaOLD-lambdaCUR))<eps1) {
      #test=0
      #print(max(abs(lambdaOLD-lambdaCUR)))
      Flambda=lambdaCUR
      break
      
    } # end if 
    
  } # end while
  # print(Flambda)
  return(list(Flambda = Flambda, Z_sups = Z_sups,W_0 = W_0))
}# end of function