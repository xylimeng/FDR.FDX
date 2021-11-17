source('OptSol.R')

library(MASS)
library(purrr)
library(nleqslv)
library(mixtools)
library(scalreg)
library(natural)
library(svd)
library(selectiveInference)
library(mvtnorm)

#helper func
indOfMaxCalpha=function(calpha.seq, numRejs, numRej.target){
  result=which.max(numRejs<numRej.target)-1
  if(result==0){
    return(length(calpha.seq))
  }
  return(result)
}


#helper func
EabsNormal=function(mu,sd,medianBetahat){
  return(sd*sqrt(2/pi)*exp(-mu^2/2/sd^2)+mu*(1-2*pnorm(-mu/sd))-medianBetahat)
}

upBound_FDPcov=function(lambdas,c.uniq,numRej.uniq,deltas.uniq,beta.hat,covMat,orders,p,fdr.target){
  #pre-screening step: reduce the times of conducting GMM
  nullBias_est=nleqslv(0.001, partial(EabsNormal, sd=mean(sqrt(diag(covMat))),medianBetahat=median(abs(beta.hat))))$x
  superset_Null=which((beta.hat<(nullBias_est+qnorm(fdr.target/2/p,lower.tail = F)*sqrt(diag(covMat)))) & (beta.hat>(nullBias_est-qnorm(fdr.target/2/p,lower.tail = F)*sqrt(diag(covMat)))))
  
  max_nrej_allowed=(p-length(superset_Null))*(1+fdr.target)
  start_ind=which(numRej.uniq<=max_nrej_allowed)[1]
  if(is.na(start_ind)){return(0)}
  for (i in start_ind:length(c.uniq)) {
    numRej.tp=numRej.uniq[i]
    c.tp=c.uniq[i]
    delta.tp=deltas.uniq[i]
    maxFDallowed=floor(numRej.tp*fdr.target)
    if(length(superset_Null)==p){next}
    est_FP=sum(abs(beta.hat[superset_Null])>lambdas[numRej.tp]*delta.tp)
    if(est_FP<maxFDallowed |est_FP==0){
      #refine the superset
      superset_Null_update=which((beta.hat<(nullBias_est+qnorm((fdr.target-est_FP/numRej.tp)/2/p,lower.tail = F)*sqrt(diag(covMat)))) & 
                                   (beta.hat>(nullBias_est-qnorm((fdr.target-est_FP/numRej.tp)/2/p,lower.tail = F)*sqrt(diag(covMat)))))
      set.seed(123)
      threshold_list=matrix(nrow = 2,ncol=2)
      k=1
      while(k<=2) {
        tryCatch(
          {mixmdl=normalmixEM(beta.hat/sqrt(diag(covMat)), k = 2,verb = F)
          mus.tp=mixmdl$mu
          sigmas.tp=mixmdl$sigma
          id=which.min(abs(mus.tp))
          threshold_list[k,]=c(mus.tp[id]-qnorm((fdr.target-est_FP/numRej.tp)/2/p,lower.tail = F)*sigmas.tp[id],
                               mus.tp[id]+qnorm((fdr.target-est_FP/numRej.tp)/2/p,lower.tail = F)*sigmas.tp[id])
          k=k+1
          },
          error=function(cond){
            #message(cond)
            threshold_list[k,]=c(NA,NA)
            return(NA)}
        )
      }
      superset_Null_update2=which((beta.hat>threshold_list[1,1]*sqrt(diag(covMat)))&(beta.hat<threshold_list[1,2]*sqrt(diag(covMat))))
      superset_Null_update3=which((beta.hat>threshold_list[2,1]*sqrt(diag(covMat)))&(beta.hat<threshold_list[2,2]*sqrt(diag(covMat))))
      if(length(superset_Null_update2)>=length(superset_Null_update3)){
        superset_Null_update=union(superset_Null_update,superset_Null_update2)
      }else{
        superset_Null_update=union(superset_Null_update,superset_Null_update3)
      }
      
      est_FP=sum(abs(beta.hat[superset_Null_update])>lambdas[numRej.tp]*delta.tp)
      if(est_FP<=maxFDallowed){
        return(numRej.tp)
      }
    }
    
  }
  return(0)
}


#helper function
get_numFD_missby_superset=function(aid_rvs,superset_tail){
  corrBernous=rowSums(aid_rvs>qnorm(p=superset_tail,lower.tail = F))
  return(corrBernous)
}

upBound_FDPcov_new=function(lambdas,c.uniq,numRej.uniq,deltas.uniq,beta.hat,covMat,orders,p,fdr.target,FDP_target=0){
  nullBias_est=nleqslv(0.001, partial(EabsNormal, sd=mean(sqrt(diag(covMat))),medianBetahat=median(abs(beta.hat))))$x
  tailprobs=seq(1/p,fdr.target/2/p,length.out=20)
  power_record=rep(0,length(tailprobs))
  fdp_controlFreq_record=rep(1,length(tailprobs))
  
  #the corr matrix of beta.hat
  corrMat=diag(1/sqrt(diag(covMat)))%*%covMat%*%diag(1/sqrt(diag(covMat)))
  set.seed(123)
  aid_rvs=abs(rmvnorm(n=2000,sigma =corrMat))
  
  for (t in 1:length(tailprobs)) {
    superset_Null=which((beta.hat<(nullBias_est+qnorm(tailprobs[t],lower.tail = F)*sqrt(diag(covMat)))) & 
                          (beta.hat>(nullBias_est-qnorm(tailprobs[t],lower.tail = F)*sqrt(diag(covMat)))))
    set.seed(123)
    threshold_list=matrix(nrow = 2,ncol=2)
    k=1
    while(k<=2) {
      tryCatch(
        {mixmdl=normalmixEM(beta.hat/sqrt(diag(covMat)), k = 2,verb = F)
        mus.tp=mixmdl$mu
        sigmas.tp=mixmdl$sigma
        id=which.min(abs(mus.tp))
        threshold_list[k,]=c(mus.tp[id]-qnorm(tailprobs[t],lower.tail = F)*sigmas.tp[id],
                             mus.tp[id]+qnorm(tailprobs[t],lower.tail = F)*sigmas.tp[id])
        k=k+1
        },
        error=function(cond){
          threshold_list[k,]=c(NA,NA)
          return(NA)}
      )
    }
    superset_Null2=which((beta.hat>threshold_list[1,1]*sqrt(diag(covMat)))&(beta.hat<threshold_list[1,2]*sqrt(diag(covMat))))
    superset_Null3=which((beta.hat>threshold_list[2,1]*sqrt(diag(covMat)))&(beta.hat<threshold_list[2,2]*sqrt(diag(covMat))))
    if(length(superset_Null2)>=length(superset_Null3)){
      superset_Null=union(superset_Null,superset_Null2)
    }else{
      superset_Null=union(superset_Null,superset_Null3)
    }
    
    max_nrej_allowed=(p-length(superset_Null))*(1+fdr.target)
    if(max_nrej_allowed<max(power_record)){
      break
    }
    start_ind=which(numRej.uniq<=max_nrej_allowed)[1]
    if(is.na(start_ind)){next}
    
    #--an upper bound on the number of FDs that are missed by the superset
    prob_cutpoint=seq(0,1,length.out=50)
    numFD_missby_super_draws=get_numFD_missby_superset(aid_rvs,superset_tail=tailprobs[t])
    
    for (i in start_ind:length(c.uniq)) {
      numRej.tp=numRej.uniq[i]
      c.tp=c.uniq[i]
      delta.tp=deltas.uniq[i]
      numFD_in_super=sum(abs(beta.hat[superset_Null])>lambdas[numRej.tp]*delta.tp)
      discrete_integrals=rep(NA,length(prob_cutpoint))
      for (j in 1:length(prob_cutpoint)) {
        discrete_integrals[j]=mean(numFD_missby_super_draws>=(prob_cutpoint[j]*numRej.tp-numFD_in_super))
      }
      fdr_est=mean(discrete_integrals)
      fdp_controlFreq_est=1-discrete_integrals[which(prob_cutpoint>fdr.target)[1]-1]
      
      if(fdr_est<=fdr.target & fdp_controlFreq_est>=FDP_target){
        power_record[t]=numRej.tp
        fdp_controlFreq_record[t]=fdp_controlFreq_est
        break
      }
    }
  }
  
  #identify the one with the highest power
  numRej.chosen=max(power_record[which(fdp_controlFreq_record>=FDP_target)])
  return(numRej.chosen)
}

#main func: FDP control fit func
FDPcontrolFit=function(x,y,lambdas,fdr.target,FDP_target=0,sigma=NA,version='fdpc+',isOrthogonal=F,standardize=F){
  
  # ----Input Parameter:
  #     x: the design matrix containing predictors
  #     y: the response variable
  #     fdr.target: target FDR level
  #     FDP_target: FDP exceedance control target, i.e. Pr(FDP<=fdr.target)>=FDP_target
  #     version: either 'fdpc' or 'fdpc+'
  #     isOrthogonal: whether X is orthogonal design
  #     standardize: whether columns in X should be standardized 
  #     lambdas: the decresing regulaziation sequence
  #     sigma: the noise level of the normal error term; default unknown(NA).
  
  # ----Output:
  #     index.rej: index of the selected predictors
  
  p=ncol(x)
  n=nrow(x)
  if(standardize){
    x=scale(x)/sqrt(n)
    y=scale(y,scale = F)
  }
  #estimate sigma if not provided
  if(!is.numeric(sigma)){
    set.seed(123)
    #use organic lasso estimate for sigma
    sigma=tryCatch(
      {olasso_cv(x = x, y =y)$sig_obj
      },error=function(cond){
        estimateSigma(x,y,intercept = F,standardize = F)$sigmahat
      }
    )
  }
  
  x=x*sqrt(n)
  if(isOrthogonal){
    tau=0.01
    beta.hat=1/(n+tau)*t(x)%*%y
  }else{
    tryCatch(
      {svdxx=svd(t(x)%*%x)
      },
      error=function(cond){
        svdxx=propack.svd(t(x)%*%x)}
    )
    #choose tau
    tau_list=n^seq(max(0.2,min(log(20,n),0.9)),0.95,length.out = 40)
    Diff_list=rep(0,length(tau_list))
    for (i in 1:length(tau_list)) {
      betahat.tp=abs(svdxx$u%*%diag(1/(svdxx$d+tau_list[i]))%*%t(svdxx$u)%*%t(x)%*%y)
      nullSize=median(betahat.tp)
      Vsize=max(sqrt(diag(svdxx$u%*%diag(svdxx$d/(svdxx$d+tau_list[i])^2)%*%t(svdxx$u))))*sigma
      Esize=median(betahat.tp[which(betahat.tp>(nullSize+1.645*Vsize))])
      Diff_list[i]=(Esize-nullSize)/Vsize
    }
    
    #choose tau that maximize the difference btw posterior mean signals and nulls
    if(all(is.na(Diff_list))){
      return(c("index.rej"=0))
      }
    tau=tau_list[which(Diff_list==max(Diff_list,na.rm = T))]
    if(is.na(tau)){
      return(c("index.rej"=0))
      }
    beta.hat=svdxx$u%*%diag(1/(svdxx$d+tau))%*%t(svdxx$u)%*%t(x)%*%y
    covMat=svdxx$u%*%diag(svdxx$d/(svdxx$d+tau)^2)%*%t(svdxx$u)*sigma^2
  }
  
  #First, solve the constrained optimization
  orders=order(abs(beta.hat),decreasing = TRUE)
  by=0.1
  calpha.seq=seq(from=0.01, to=sum(beta.hat^2)*(n+tau)/sigma^2-1,by=by)
  results=numRej_forAllC(abs(beta.hat[orders]),lambdas,n,sigma,tau,endpoint=calpha.seq[length(calpha.seq)], by)
  numRejs=results$nRejs
  deltas=results$deltas
  
  numRej.uniq=unique(numRejs)
  calpha.ids=sapply(numRej.uniq, indOfMaxCalpha, calpha.seq=calpha.seq,numRejs=numRejs)
  calpha.uniq=calpha.seq[calpha.ids]
  deltas.uniq=deltas[calpha.ids]
  
  c.uniq=calpha.uniq*sigma^2/(n+tau)
  if(version=='fdpc'){
    numRej.chosen=upBound_FDPcov(lambdas,c.uniq,numRej.uniq,deltas.uniq,beta.hat,covMat,orders,p,fdr.target)
  }else{
    numRej.chosen=upBound_FDPcov_new(lambdas,c.uniq,numRej.uniq,deltas.uniq,beta.hat,covMat,orders,p,
                                     fdr.target,FDP_target)
  }
  
  if(numRej.chosen!=0){
    index.rej=c(1:p)[orders[1:numRej.chosen]]
  }else{
    index.rej=0
  }
  return(c("index.rej"=index.rej))
}





