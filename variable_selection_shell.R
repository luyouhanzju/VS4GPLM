setwd("~/VSGeneralizedPartialLinearModel/latest codes")
source("local_scoring.R")
source("examples.R")
source("variable_selection.R")
source("main_2tau_deviance.R")
source("getVS.R")
mymain <- function(irep){
  
  tryCatch(
    expr = {
      #tauall_1= seq(1,45,5)
      # Ex1/Ex2: tauall_1= seq(1,45,5), tauall_2 = seq(0.1,5,1) -------
      # Ex1: tauall_1= seq(1,45,5), tauall_2 = seq(0.1,3,0.5)
      # Ex2: tauall_1= seq(1,45,5), tauall_2 = seq(0.1,3,0.5)
      # Ex3: tauall_1= seq(1,200,5), tauall_2 = seq(0.1,5,1)
      # Ex4: and Ex 5 seq(1,350,10); seq(0.1,8,1); Ex 6 and Ex7 seq(1,600,10)
      tauall_1 = seq(1,300,30)
      tauall_2 = seq(0.1,5,1)
      
      SelModNaive=NULL
      SelModProfile=NULL
      samp=NULL
      set.seed(irep)
      #print("seed")
      exampleIdx = 1
      traindata=gendata(exampleIdx)

      Xall = traindata$x
      Uall = traindata$u
      Yall = traindata$y

      if (exampleIdx == 3)
      {mea = traindata$mea
      }
      p = ncol(Xall)
      q = ncol(Uall)
      DDall = sapply(1:q, function(j) as.matrix(dist(Uall[,j])), simplify=FALSE)
      tau2allmain = seq(0.1,5,1)
      res = variable_selection(tauall_1,tauall_2,Xall,DDall,Yall,Uall,exampleIdx,tau2allmain)
      print('res completed!')
      wuNaive = res$wuNaive
      print('wuNaive!')
      #wuProfile = res$wuProfile
      wuRefitting = res$wuRefitting
      SelModNaive = rbind(SelModNaive,wuNaive$lamhatALL[which.min(wuNaive$BIC),]>1e-2)
      print('SelModNaive!')
      print(wuNaive$lamhatALL[which.min(wuNaive$BIC),]>1e-2)
      #SelModNaive
      SelModNaive_list = list(SelModNaive = SelModNaive)
      tau2allmain_list = list(tau2allmain = tau2allmain)
      print(SelModProfile)
      result = append(res,SelModNaive_list)
      result = append(result,SelModNaive_list)
      result = append(result,SelModNaive_list)
      # if (exampleIdx == 3)
      # {
      #   result = append(result,traindata)
      # }

      return(result)
    },
    error = function(e){
      result <- vector('list', 7)
      names(result) <- c("wuNaive","tauall_1","tauall_2","SelModNaive","SelModNaive","SelModNaive")
      print('error')
      return(result)
    },
    finally = {
      print(paste('irep:', irep, 'done'))
    }

  )
}

seed.start = 170

result = list()

library(foreach)
library(quadprog)
library(iterators)
library(doParallel)
registerDoParallel(20)

samp=foreach(irep=1:170, .verbose=TRUE, .combine=c, .packages=c("matlab", "psych")) %dopar% mymain(irep + seed.start)
#save(file=sprintf("diabetes_tau2allmain_01_21_02%d.Rdata",seed.start), list = ls())
##### after to 2023/01/25, we use bic as log likelihood
#save(file=sprintf("TEST_bic_GPLM_Example3_%d_170_20231010_newBIC_tau1_refitting_n300_U13_2.Rdata",seed.start), list = ls())
save(file=sprintf("GPLM_Example1_20241208.Rdata",seed.start), list = ls())
stopImplicitCluster()
