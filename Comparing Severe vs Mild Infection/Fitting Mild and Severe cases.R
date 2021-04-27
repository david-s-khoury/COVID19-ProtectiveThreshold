library('ggplot2')
library('plyr')
library('car')

#### Set directory to same directory as the r-script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Import Tables For fitting to severe cases only
SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe<-read.csv("SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe.csv", fileEncoding="UTF-8-BOM") 


### The Neutralisation ratio of vaccine to convalescence using reported neut titres.
SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutConv)


### Import Tables For fitting Models (mild infection cases)
SummaryTable_Efficacy_NeutRatio_SD_SEM<-read.csv("SummaryTable_Efficacy_NeutRatio_SD_SEM.csv", fileEncoding="UTF-8-BOM") 


### The Neutralisation ratio of vaccine to convalescence using reported neut titres.
SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)



##########################################################################################
##################### The Models

####Logistic Model
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  Output<-NULL
  for (i in 1:length(mu_titre)) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk,C50)*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}

#############################################################################################################################
#######  Logistic model for Raw Efficacy Counts (Different parameters for both severe and mild infection) ###################

FittingLogistic_Raw_Combined_Diff_All<-function(logRisk0_Severe,logk_Severe,C50_Severe,logRisk0,logk,C50,N_C_Severe,N_V_Severe,Inf_C_Severe,Inf_V_Severe,MeanVector_Severe,SDVector_Severe,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector)
{
  
  Risk0=exp(logRisk0_Severe)
  
  if (length(SDVector_Severe)==1) {
    SDVector_Severe=rep(SDVector_Severe,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C_Severe) | is.na(MeanVector_Severe) | is.na(SDVector_Severe))
  N_C_SEVERE=N_C_Severe[!IndexNA]
  N_V_SEVERE=N_V_Severe[!IndexNA]
  Inf_V_Severe=Inf_V_Severe[!IndexNA]
  Inf_C_Severe=Inf_C_Severe[!IndexNA]
  MeanVector_Severe=MeanVector_Severe[!IndexNA]
  SDVector_Severe=SDVector_Severe[!IndexNA]
  
  LL_Severe=0
  for (i in 1:length(N_C_SEVERE)) {
    
    LL_Severe=LL_Severe-log(dbinom(Inf_C_Severe[i],N_C_Severe[i],Risk0[i]))-log(dbinom(Inf_V_Severe[i],N_V_Severe[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector_Severe[i],SDVector_Severe[i],logk_Severe,C50_Severe))))
  }
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk,C50))))
  }
  
  
  LL_TOT = LL+LL_Severe
  
}


#Initial values (just pick random initial values)
LogisticEstimate=c("logk_Severe"=log(runif(1,0,1)),"C50_Severe"=log10(runif(1,0,0.05)),"logk"=log(runif(1,0,1)),"C50"=log10(runif(1,0,0.5)) )##InitialGues


#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different parameter for each. 
#Will use reported mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_All(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+1],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+2],
                                        p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+8],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+8],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+8],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)


#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different parameter for each. 
#Will use estimated mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_All(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+1],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+2],
                                        p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+8],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+8],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+8],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)



####################################################################################################################################
#######  Logistic model for Raw Efficacy Counts (Different n50 for both severe and mild infection; shared slope) ###################

FittingLogistic_Raw_Combined_Diff_EC50_Only<-function(logRisk0_Severe,C50_Severe,logRisk0,logk,C50,N_C_Severe,N_V_Severe,Inf_C_Severe,Inf_V_Severe,MeanVector_Severe,SDVector_Severe,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector)
{
  
  Risk0=exp(logRisk0_Severe)
  
  if (length(SDVector_Severe)==1) {
    SDVector_Severe=rep(SDVector_Severe,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C_Severe) | is.na(MeanVector_Severe) | is.na(SDVector_Severe))
  N_C_SEVERE=N_C_Severe[!IndexNA]
  N_V_SEVERE=N_V_Severe[!IndexNA]
  Inf_V_Severe=Inf_V_Severe[!IndexNA]
  Inf_C_Severe=Inf_C_Severe[!IndexNA]
  MeanVector_Severe=MeanVector_Severe[!IndexNA]
  SDVector_Severe=SDVector_Severe[!IndexNA]
  
  
  LL_Severe=0
  for (i in 1:length(N_C_SEVERE)) {
    
    LL_Severe=LL_Severe-log(dbinom(Inf_C_Severe[i],N_C_Severe[i],Risk0[i]))-log(dbinom(Inf_V_Severe[i],N_V_Severe[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector_Severe[i],SDVector_Severe[i],logk,C50_Severe))))
  }
  
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk,C50))))
  }
  
  
  LL_TOT = LL+LL_Severe
  
}


#Initial values (just pick random initial values)
LogisticEstimate=c("C50_Severe"=log10(runif(1,0,0.001)),"logk"=log(runif(1,0,1)),"C50"=log10(runif(1,0,0.5)) )

#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different n50 for each; but with a shared slope. 
#Will use reported mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_EC50_Only(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                              p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+1],
                                              p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+7],
                                              p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+7],
                                              p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+7],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)


#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different n50 for each; but with a shared slope. 
#Will use estimated mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_EC50_Only(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                              p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+1],
                                              p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+7],
                                              p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+7],
                                              p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+7],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                              SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)





####################################################################################################################################
#######  Logistic model for Raw Efficacy Counts (Different slope for both severe and mild infection; shared n50) ###################

FittingLogistic_Raw_Combined_Diff_Slope_Only<-function(logRisk0_Severe,logk_Severe,logRisk0,logk,C50,N_C_Severe,N_V_Severe,Inf_C_Severe,Inf_V_Severe,MeanVector_Severe,SDVector_Severe,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector)
{
  Risk0=exp(logRisk0_Severe)
  
  if (length(SDVector_Severe)==1) {
    SDVector_Severe=rep(SDVector_Severe,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C_Severe) | is.na(MeanVector_Severe) | is.na(SDVector_Severe))
  N_C_SEVERE=N_C_Severe[!IndexNA]
  N_V_SEVERE=N_V_Severe[!IndexNA]
  Inf_V_Severe=Inf_V_Severe[!IndexNA]
  Inf_C_Severe=Inf_C_Severe[!IndexNA]
  MeanVector_Severe=MeanVector_Severe[!IndexNA]
  SDVector_Severe=SDVector_Severe[!IndexNA]
  
  LL_Severe=0
  for (i in 1:length(N_C_SEVERE)) {
    
    LL_Severe=LL_Severe-log(dbinom(Inf_C_Severe[i],N_C_Severe[i],Risk0[i]))-log(dbinom(Inf_V_Severe[i],N_V_Severe[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector_Severe[i],SDVector_Severe[i],logk_Severe,C50))))
  }
  
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk,C50))))
  }
  
  
  LL_TOT = LL+LL_Severe
  
}

#Initial values (just pick random initial values)
LogisticEstimate=c("logk_Severe"=log(runif(1,0,0.001)),"logk"=log(runif(1,0,1)),"C50"=log10(runif(1,0,0.5)) )

#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different slope for each; but with a shared n50. 
#Will use reported mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_Slope_Only(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                               p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+1],
                                               p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+7],
                                               p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+7],
                                               p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+7],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)


#Minimize the negative of the log-likelihood value to fit both severe and mild cases with different slope for each; but with a shared n50. 
#Will use estimated mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only<-nlm(function(p){
  FittingLogistic_Raw_Combined_Diff_Slope_Only(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                               p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))+1],
                                               p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+7],
                                               p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+7],
                                               p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+7],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                               SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)




####################################################################################################################################
#######  Logistic model for Raw Efficacy Counts (shared n50 and slope for both) ####################################################

FittingLogistic_Raw_Combined_Same_All<-function(logRisk0_Severe,logRisk0,logk,C50,N_C_Severe,N_V_Severe,Inf_C_Severe,Inf_V_Severe,MeanVector_Severe,SDVector_Severe,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector)
{
  
  Risk0=exp(logRisk0_Severe)
  
  if (length(SDVector_Severe)==1) {
    SDVector_Severe=rep(SDVector_Severe,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C_Severe) | is.na(MeanVector_Severe) | is.na(SDVector_Severe))
  N_C_SEVERE=N_C_Severe[!IndexNA]
  N_V_SEVERE=N_V_Severe[!IndexNA]
  Inf_V_Severe=Inf_V_Severe[!IndexNA]
  Inf_C_Severe=Inf_C_Severe[!IndexNA]
  MeanVector_Severe=MeanVector_Severe[!IndexNA]
  SDVector_Severe=SDVector_Severe[!IndexNA]
  
  
  LL_Severe=0
  for (i in 1:length(N_C_SEVERE)) {
    
    LL_Severe=LL_Severe-log(dbinom(Inf_C_Severe[i],N_C_Severe[i],Risk0[i]))-log(dbinom(Inf_V_Severe[i],N_V_Severe[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector_Severe[i],SDVector_Severe[i],logk,C50))))
  }
  
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(Efficacy))
  }
  
  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk,C50))))
  }
  
  
  LL_TOT = LL+LL_Severe
  
}

#Initial values (just pick random initial values)
LogisticEstimate=c("logk"=log(runif(1,0,0.5)),"C50"=log10(runif(1,0,0.5)) )

#Minimize the negative of the log-likelihood value to fit both severe and mild cases with shared parameters (both n50 and the slope) for each. 
#Will use reported mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All<-nlm(function(p){
  FittingLogistic_Raw_Combined_Same_All(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                        p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+6],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+6],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+6],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)

#Minimize the negative of the log-likelihood value to fit both severe and mild cases with shared parameters (both n50 and the slope) for each. 
#Will use estimated mean of titre from the original studies.
FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All<-nlm(function(p){
  FittingLogistic_Raw_Combined_Same_All(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont))],
                                        p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+6],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1+6],
                                        p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2+6],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)], 
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                        SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},
  c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe$NumCont)]),log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate)
  ,hessian=TRUE,iterlim=1000)


################### Table for AIC and p-value from the likelihood ratio test###################################################
npar_thres<-c(4,3,3,2)

LL_FittedLogistic_RawEfficacy_MeanCens_SDPool<-c(-FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$minimum,-FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$minimum,-FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$minimum,-FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$minimum)
LL_FittedLogistic_RawEfficacy_MeanRept_SDPool<-c(-FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$minimum,-FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$minimum,-FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$minimum,-FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$minimum)

AIC_FittedLogistic_RawEfficacy_MeanCens_SDPool<-2*npar_thres - 2*LL_FittedLogistic_RawEfficacy_MeanCens_SDPool
AIC_FittedLogistic_RawEfficacy_MeanRept_SDPool<-2*npar_thres - 2*LL_FittedLogistic_RawEfficacy_MeanRept_SDPool

Delta_AIC_FittedLogistic_RawEfficacy_MeanCens_SDPool<-AIC_FittedLogistic_RawEfficacy_MeanCens_SDPool - min(AIC_FittedLogistic_RawEfficacy_MeanCens_SDPool)
Delta_AIC_FittedLogistic_RawEfficacy_MeanRept_SDPool<-AIC_FittedLogistic_RawEfficacy_MeanRept_SDPool - min(AIC_FittedLogistic_RawEfficacy_MeanRept_SDPool)

p_val_LLR_test_FittedLogistic_RawEfficacy_MeanCens_SDPool<-NULL
p_val_LLR_test_FittedLogistic_RawEfficacy_MeanRept_SDPool<-NULL

llr = 2*(LL_FittedLogistic_RawEfficacy_MeanCens_SDPool[1] - LL_FittedLogistic_RawEfficacy_MeanCens_SDPool[2])
p_val_LLR_test_FittedLogistic_RawEfficacy_MeanCens_SDPool[1]<-pchisq(llr, df=1, lower.tail=FALSE)
llr = 2*(LL_FittedLogistic_RawEfficacy_MeanCens_SDPool[2] - LL_FittedLogistic_RawEfficacy_MeanCens_SDPool[4])
p_val_LLR_test_FittedLogistic_RawEfficacy_MeanCens_SDPool[4]<-pchisq(llr, df=1, lower.tail=FALSE)

llr = 2*(LL_FittedLogistic_RawEfficacy_MeanRept_SDPool[1] - LL_FittedLogistic_RawEfficacy_MeanRept_SDPool[2])
p_val_LLR_test_FittedLogistic_RawEfficacy_MeanRept_SDPool[1]<-pchisq(llr, df=1, lower.tail=FALSE)
llr = 2*(LL_FittedLogistic_RawEfficacy_MeanRept_SDPool[2] - LL_FittedLogistic_RawEfficacy_MeanRept_SDPool[4])
p_val_LLR_test_FittedLogistic_RawEfficacy_MeanRept_SDPool[4]<-pchisq(llr, df=1, lower.tail=FALSE)

Logistic_RawEfficacy_MeanRept_SDPool<-cbind(LL_FittedLogistic_RawEfficacy_MeanRept_SDPool,AIC_FittedLogistic_RawEfficacy_MeanRept_SDPool,Delta_AIC_FittedLogistic_RawEfficacy_MeanRept_SDPool,p_val_LLR_test_FittedLogistic_RawEfficacy_MeanRept_SDPool)

table_llr<-Logistic_RawEfficacy_MeanRept_SDPool
rownames(table_llr)<-c("All Diff","EC50 Only","Slope Only","All Same")
colnames(table_llr) = c("Likelihood","AIC","Delta_AIC", "p_val")

write.csv(table_llr,file="Table_AIC_Logistic.csv")



################### Table for the estimated parameters from each model ###################################################
parLogistic_MeanCens_SDPool_Diff_All <- c(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$estimate[7],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$estimate[8],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$estimate[17],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$estimate[18])
parLogistic_MeanCens_SDPool_Diff_EC50_Only <-c(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$estimate[16],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$estimate[16],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$estimate[7],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$estimate[17])
parLogistic_MeanCens_SDPool_Diff_Slope_Only <-c(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$estimate[7],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$estimate[16],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$estimate[17],FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$estimate[17])
parLogistic_MeanCens_SDPool_Same_All <-c(FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$estimate[15],FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$estimate[16])

parLogistic_MeanRept_SDPool_Diff_All <- c(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$estimate[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$estimate[8],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$estimate[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$estimate[18])
parLogistic_MeanRept_SDPool_Diff_EC50_Only <-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[17])
parLogistic_MeanRept_SDPool_Diff_Slope_Only <-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$estimate[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$estimate[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$estimate[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$estimate[17])
parLogistic_MeanRept_SDPool_Same_All <-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$estimate[15],FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$estimate[16])

table_par_Logistic<-rbind(parLogistic_MeanRept_SDPool_Diff_All,parLogistic_MeanRept_SDPool_Diff_EC50_Only,parLogistic_MeanRept_SDPool_Diff_Slope_Only,parLogistic_MeanRept_SDPool_Same_All
)
colnames(table_par_Logistic)=c("Severe Slope","Severe EC50", "Mild Slope", "Mild EC50")

write.csv(table_par_Logistic,file="Table_Par_Logistic.csv")


################### Table for the Lower and upper 95% CI for the estimated parameters ###################################################
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_All$estimate)
FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_Diff_All<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[8],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[8],
                                                          FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[18],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[18])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_EC50_Only$estimate)
FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_Diff_EC50_Only<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16],
                                                                FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16],
                                                                FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[7],
                                                                FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Diff_Slope_Only$estimate)
FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_Diff_Slope_Only<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16],
                                                                 FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanCens_SDPool_Same_All$estimate)
FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_Same_All<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[15],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[15],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_All$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Diff_All<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[8],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[8],
                                                          FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[18],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[18])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Diff_EC50_Only<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16],
                                                                FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16],
                                                                FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[7],
                                                                FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_Slope_Only$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Diff_Slope_Only<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[7],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16],
                                                                 FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[17],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[17])

FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper<- (sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower<- -(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$hessian))))*1.96 +(FittedLogistic_RawEfficacy_MeanRept_SDPool_Same_All$estimate)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Same_All<-c(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[15],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[15],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Lower[16],FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Upper[16])

table_par_CI_Logistic<-rbind(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Diff_All,FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Diff_EC50_Only, FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Diff_Slope_Only, FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_Same_All)
colnames(table_par_CI_Logistic) = c("Slope Severe Lower","Slope Severe Upper", "EC50 Severe Lower","EC50 Severe Upper","Slope Mild Lower", "Slope Mild Upper", "EC50 Mild Lower", "EC50 Mild Upper")

write.csv(table_par_CI_Logistic,file="Table_CI_Logistic.csv")