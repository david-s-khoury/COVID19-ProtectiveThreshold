library('ggplot2')
library('plyr')
library('car')
library('maxLik') #To compute numerical gradient as a function of the parameters 

#### Set directory to same directory as the r-script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Import Tables For fitting Models
SummaryTable_Efficacy_NeutRatio_SD_SEM<-read.csv("SummaryTable_Efficacy_NeutRatio_SD_SEM.csv") ##For threshold and logistic
IndividualNeutData_NormalisedbyConv<-read.csv("IndividualNeutData_NormalisedbyConv.csv") ##For PNC mdoel


### The Neutralisation ratio of vaccine to convalescence using reported neut titres.
SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutMean/SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)



##########################################################################################
##################### The Models

####Logistic Model
ProbRemainUninfected=function(logTitre,logk,C50){1/(1+exp(-exp(logk)*(logTitre-C50)))}

LogisticModel_PercentUninfected=function(mu_titre,sig_titre,logk,C50){
  Output<-NULL
  
  if (length(C50)==1) {
    C50=rep(C50,length(mu_titre))
  }
  
  if (length(logk)==1) {
    logk=rep(logk,length(mu_titre))
  }
  
  for (i in 1:length(mu_titre)) {
    Step=sig_titre[i]*0.001
    IntegralVector=seq(mu_titre[i]-5*sig_titre[i],mu_titre[i]+5*sig_titre[i],by=Step)
    Output[i]=sum(ProbRemainUninfected(IntegralVector,logk[i],C50[i])*dnorm(IntegralVector,mu_titre[i],sig_titre[i]))*Step
  }
  Output
}


### Logistic model for Raw Efficacy Counts
FittingLogistic_Raw<-function(logRisk0,logk,C50,N_C,N_V,Inf_C,Inf_V,MeanVector,SDVector){
  
  Risk0=exp(logRisk0)
  
  if (length(SDVector)==1) {
    SDVector=rep(SDVector,length(N_C))
  }

  IndexNA=(is.na(N_C) | is.na(MeanVector) | is.na(SDVector))
  N_C=N_C[!IndexNA]
  N_V=N_V[!IndexNA]
  Inf_V=Inf_V[!IndexNA]
  Inf_C=Inf_C[!IndexNA]
  MeanVector=MeanVector[!IndexNA]
  SDVector=SDVector[!IndexNA]
  
  if (length(C50)==1) {
    C50=rep(C50,length(N_C))
  }
  
  if (length(logk)==1) {
    logk=rep(logk,length(N_C))
  }
  
  LL=0
  for (i in 1:length(N_C)) {
    
    LL=LL-log(dbinom(Inf_C[i],N_C[i],Risk0[i]))-log(dbinom(Inf_V[i],N_V[i],Risk0[i]*(1-LogisticModel_PercentUninfected(MeanVector[i],SDVector[i],logk[i],C50[i]))))
  }
  LL
}



### Protective neut classification model
PNC_ProbObservingAboveThreshold=function(EfficacyForStudy,NumberInStudy,Values,Threshold,cens){
  if (cens<Threshold) {
    dbinom(sum(Values>Threshold),NumberInStudy,EfficacyForStudy)  
  } else if (sum(Values>cens)==length(Values)) {
    dbinom(sum(Values>Threshold),NumberInStudy,EfficacyForStudy)  
  } else {
    pbinom(sum(Values<=cens+abs(cens)*0.05),NumberInStudy,1-EfficacyForStudy) 
  }
  
}

PNC_LikelihoodObservingAbove=function(Threshold,Study,Efficacy,TableData,LOD){
  LL<-NULL
  for (i in 1:length(Study)){
    
    LL[i]=log(PNC_ProbObservingAboveThreshold(Efficacy[i],length(TableData[TableData$Study==Study[i],2]),TableData[TableData$Study==Study[i],2],Threshold,LOD[i]))
  }
  sum(LL)
  
}




##########################################################################################
##Fitting Models
LogisticEstimate=c("logk"=log(2.7),"C50"=log10(0.5))##InitialGuess

### Raw Efficacy Data Models - Fitting Logistic Model to Different combinations of Mean and SD
FittedLogistic_RawEfficacy_MeanCens_SDCens<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],    SummaryTable_Efficacy_NeutRatio_SD_SEM$SD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_RawEfficacy_MeanCens_SDMelb<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],    SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_RawEfficacy_MeanCens_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],    SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_RawEfficacy_MeanRept_SDCens<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$SD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_RawEfficacy_MeanRept_SDMelb<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)
FittedLogistic_RawEfficacy_MeanRept_SDPool<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+1],p[sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))+2],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                            SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},c(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),LogisticEstimate),hessian=TRUE)


## CI raw efficacy logistic model data
FittedLogistic_RawEfficacy_MeanCens_SDCens_CI<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDCens$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanCens_SDCens$estimate,1)
FittedLogistic_RawEfficacy_MeanCens_SDMelb_CI<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDMelb$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanCens_SDMelb$estimate,1)
FittedLogistic_RawEfficacy_MeanCens_SDPool_CI<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanCens_SDPool$estimate,1)
FittedLogistic_RawEfficacy_MeanRept_SDCens_CI<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDCens$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanRept_SDCens$estimate,1)
FittedLogistic_RawEfficacy_MeanRept_SDMelb_CI<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDMelb$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanRept_SDMelb$estimate,1)
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian))),1)*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)

FittedLogistic_RawEfficacy_MeanCens_SDCens_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDCens$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanCens_SDCens$estimate,2)[1]
FittedLogistic_RawEfficacy_MeanCens_SDMelb_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDMelb$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanCens_SDMelb$estimate,2)[1]
FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanCens_SDPool$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanCens_SDPool$estimate,2)[1]
FittedLogistic_RawEfficacy_MeanRept_SDCens_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDCens$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanRept_SDCens$estimate,2)[1]
FittedLogistic_RawEfficacy_MeanRept_SDMelb_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDMelb$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanRept_SDMelb$estimate,2)[1]
FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_logk<- tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian))),2)[1]*1.96*c(-1,1)+tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1]


### PNC (Protective Neutralisation Classification Model)

##Adjust the LOD by the mean of convalsecnce for each group
#(use the reported neut mean for convalescence group rather than censored estimate
#of mean for each group from extracted data because the censored model assumes normal
#distribution and trying to keep this "distribution free" model)
SummaryTable_Efficacy_NeutRatio_SD_SEM$LOD_adj=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$LOD)-log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)


### Create a table of the individuals log10(ratio neut titres) where each 
#individuals neut titres are scaled by the mean neut titre for convalescent 
#individuals reported in the study. (Not my estimated neut titres from censored model
# again trying to avoid assumptions about normal distribution in this model)
FittingTable<-IndividualNeutData_NormalisedbyConv[,c("Study","RatioTitres")]

###Scan Through Possible Thresholds and calculate log likelihood
PossibleThresholds=FittingTable$RatioTitres

LL<-NULL
for (i in 1:length(PossibleThresholds)) {
  LL[i]=PNC_LikelihoodObservingAbove(PossibleThresholds[i],SummaryTable_Efficacy_NeutRatio_SD_SEM$Study,SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittingTable,SummaryTable_Efficacy_NeutRatio_SD_SEM$LOD_adj)
}

##Find Best Fit
MaxLikelihood=max(LL)
PNC_Threshold=FittingTable[which.max(LL),2]


#####THIS WILL BE SLOW - Boostrap for CI
##PNC CI (Bootstrap - resample from table preserving number of data points per study)
Iteration=1000
PNC_Threshold_CIsampling<-NULL
for (j in 1:Iteration) {
  tempindexlist<-NULL
  for (jj in 1:length(SummaryTable_Efficacy_NeutRatio_SD_SEM$Study)) {
    tempindexlist[(length(tempindexlist)+1):(length(tempindexlist)+sum(FittingTable$Study==SummaryTable_Efficacy_NeutRatio_SD_SEM$Study[jj]))]=sample(which(FittingTable$Study==SummaryTable_Efficacy_NeutRatio_SD_SEM$Study[jj]),sum(FittingTable$Study==SummaryTable_Efficacy_NeutRatio_SD_SEM$Study[jj]),replace=TRUE)
  }
  
  tempFittingTable<-FittingTable[tempindexlist,]
  tempLL<-NULL
  for (i in 1:length(tempFittingTable$RatioTitres)) {
    tempLL[i]=PNC_LikelihoodObservingAbove(tempFittingTable$RatioTitres[i],SummaryTable_Efficacy_NeutRatio_SD_SEM$Study,SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,tempFittingTable,SummaryTable_Efficacy_NeutRatio_SD_SEM$LOD_adj)
  }
  PNC_Threshold_CIsampling[j]=tempFittingTable[which.max(tempLL),2]
}
CIthreshold=0.95

PNC_Threshold_CI=quantile(PNC_Threshold_CIsampling,c((1-CIthreshold)*0.5,CIthreshold+(1-CIthreshold)*0.5))



####################################################################################################################################
######################Create Summary Table Of Fitted Estimate Titres
TableOfEstimatedTitres=data.frame("Method"=c("PNT",rep("Logistic_Raw",6)),
                                  "Mean"=c("PNT","CensoredEstimate","CensoredEstimate","CensoredEstimate","Reported","Reported","Reported"),
                                  "SD"=c("PNT",rep(c("byStudy_cens","MelbStudy_cens","Pooled_cens"),2)),
                                  "Estimate_EC50"=NA,
                                  "L_CI"=NA,
                                  "U_CI"=NA,
                                  "Estimate_k"=NA,
                                  "L_CI_k"=NA,
                                  "U_CI_k"=NA)


TableOfEstimatedTitres$Estimate_EC50[1]=10^PNC_Threshold
TableOfEstimatedTitres$L_CI[1]=10^PNC_Threshold_CI[1]
TableOfEstimatedTitres$U_CI[1]=10^PNC_Threshold_CI[2]

TableOfEstimatedTitres$Estimate_EC50[2]=10^tail(FittedLogistic_RawEfficacy_MeanCens_SDCens$estimate,1)
TableOfEstimatedTitres$L_CI[2]=10^FittedLogistic_RawEfficacy_MeanCens_SDCens_CI[1]
TableOfEstimatedTitres$U_CI[2]=10^FittedLogistic_RawEfficacy_MeanCens_SDCens_CI[2]
TableOfEstimatedTitres$Estimate_k[2]=exp(tail(FittedLogistic_RawEfficacy_MeanCens_SDCens$estimate,2)[1])
TableOfEstimatedTitres$L_CI_k[2]=exp(FittedLogistic_RawEfficacy_MeanCens_SDCens_CI_logk[1])
TableOfEstimatedTitres$U_CI_k[2]=exp(FittedLogistic_RawEfficacy_MeanCens_SDCens_CI_logk[2])

TableOfEstimatedTitres$Estimate_EC50[3]=10^tail(FittedLogistic_RawEfficacy_MeanCens_SDMelb$estimate,1)
TableOfEstimatedTitres$L_CI[3]=10^FittedLogistic_RawEfficacy_MeanCens_SDMelb_CI[1]
TableOfEstimatedTitres$U_CI[3]=10^FittedLogistic_RawEfficacy_MeanCens_SDMelb_CI[2]
TableOfEstimatedTitres$Estimate_k[3]=exp(tail(FittedLogistic_RawEfficacy_MeanCens_SDMelb$estimate,2)[1])
TableOfEstimatedTitres$L_CI_k[3]=exp(FittedLogistic_RawEfficacy_MeanCens_SDMelb_CI_logk[1])
TableOfEstimatedTitres$U_CI_k[3]=exp(FittedLogistic_RawEfficacy_MeanCens_SDMelb_CI_logk[2])

TableOfEstimatedTitres$Estimate_EC50[4]=10^tail(FittedLogistic_RawEfficacy_MeanCens_SDPool$estimate,1)
TableOfEstimatedTitres$L_CI[4]=10^FittedLogistic_RawEfficacy_MeanCens_SDPool_CI[1]
TableOfEstimatedTitres$U_CI[4]=10^FittedLogistic_RawEfficacy_MeanCens_SDPool_CI[2]
TableOfEstimatedTitres$Estimate_k[4]=exp(tail(FittedLogistic_RawEfficacy_MeanCens_SDPool$estimate,2)[1])
TableOfEstimatedTitres$L_CI_k[4]=exp(FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_logk[1])
TableOfEstimatedTitres$U_CI_k[4]=exp(FittedLogistic_RawEfficacy_MeanCens_SDPool_CI_logk[2])

TableOfEstimatedTitres$Estimate_EC50[5]=10^tail(FittedLogistic_RawEfficacy_MeanRept_SDCens$estimate,1)
TableOfEstimatedTitres$L_CI[5]=10^FittedLogistic_RawEfficacy_MeanRept_SDCens_CI[1]
TableOfEstimatedTitres$U_CI[5]=10^FittedLogistic_RawEfficacy_MeanRept_SDCens_CI[2]
TableOfEstimatedTitres$Estimate_k[5]=exp(tail(FittedLogistic_RawEfficacy_MeanRept_SDCens$estimate,2)[1])
TableOfEstimatedTitres$L_CI_k[5]=exp(FittedLogistic_RawEfficacy_MeanRept_SDCens_CI_logk[1])
TableOfEstimatedTitres$U_CI_k[5]=exp(FittedLogistic_RawEfficacy_MeanRept_SDCens_CI_logk[2])

TableOfEstimatedTitres$Estimate_EC50[6]=10^tail(FittedLogistic_RawEfficacy_MeanRept_SDMelb$estimate,1)
TableOfEstimatedTitres$L_CI[6]=10^FittedLogistic_RawEfficacy_MeanRept_SDMelb_CI[1]
TableOfEstimatedTitres$U_CI[6]=10^FittedLogistic_RawEfficacy_MeanRept_SDMelb_CI[2]
TableOfEstimatedTitres$Estimate_k[6]=exp(tail(FittedLogistic_RawEfficacy_MeanRept_SDMelb$estimate,2)[1])
TableOfEstimatedTitres$L_CI_k[6]=exp(FittedLogistic_RawEfficacy_MeanRept_SDMelb_CI_logk[1])
TableOfEstimatedTitres$U_CI_k[6]=exp(FittedLogistic_RawEfficacy_MeanRept_SDMelb_CI_logk[2])

TableOfEstimatedTitres$Estimate_EC50[7]=10^tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)
TableOfEstimatedTitres$L_CI[7]=10^FittedLogistic_RawEfficacy_MeanRept_SDPool_CI[1]
TableOfEstimatedTitres$U_CI[7]=10^FittedLogistic_RawEfficacy_MeanRept_SDPool_CI[2]
TableOfEstimatedTitres$Estimate_k[7]=exp(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1])
TableOfEstimatedTitres$L_CI_k[7]=exp(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_logk[1])
TableOfEstimatedTitres$U_CI_k[7]=exp(FittedLogistic_RawEfficacy_MeanRept_SDPool_CI_logk[2])


write.csv(TableOfEstimatedTitres,file="TableOfEstimatedTitres.csv")



############################################################################################
####################Plotting Validation/Quality Of Fit Plots
###Logistic Model Raw Efficacy Data
FittedLogistic_RawEfficacy_MeanCens_SDCens_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens   ,SummaryTable_Efficacy_NeutRatio_SD_SEM$SD       ,tail(FittedLogistic_RawEfficacy_MeanCens_SDCens$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanCens_SDCens$estimate,1))
FittedLogistic_RawEfficacy_MeanCens_SDMelb_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens    ,SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD  ,tail(FittedLogistic_RawEfficacy_MeanCens_SDMelb$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanCens_SDMelb$estimate,1))
FittedLogistic_RawEfficacy_MeanCens_SDPool_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens    ,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD,tail(FittedLogistic_RawEfficacy_MeanCens_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanCens_SDPool$estimate,1))
FittedLogistic_RawEfficacy_MeanRept_SDCens_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM$SD      ,tail(FittedLogistic_RawEfficacy_MeanRept_SDCens$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDCens$estimate,1))
FittedLogistic_RawEfficacy_MeanRept_SDMelb_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD  ,tail(FittedLogistic_RawEfficacy_MeanRept_SDMelb$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDMelb$estimate,1))
FittedLogistic_RawEfficacy_MeanRept_SDPool_ModelOutput<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD,tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1))


pdf("LogisticModel_RawEfficacy_QualityOfFit.pdf",height=7,width=5)
par(mfrow=c(3,2))
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittedLogistic_RawEfficacy_MeanCens_SDCens_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Cens. mean, SD by study" )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittedLogistic_RawEfficacy_MeanCens_SDMelb_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Cens. mean, SD Melb." )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittedLogistic_RawEfficacy_MeanCens_SDPool_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Cens. mean, SD Pool" )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittedLogistic_RawEfficacy_MeanRept_SDCens_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD by study" )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittedLogistic_RawEfficacy_MeanRept_SDMelb_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD Melb." )
abline(0,1,lty=2)
plot(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,FittedLogistic_RawEfficacy_MeanRept_SDPool_ModelOutput,xlab="Observed efficacy",ylab="Estimated efficacy",main="Reported mean, SD Pool" )
abline(0,1,lty=2)
dev.off()



#### PNC model

pdf("PNCModelQualityOfFit.pdf",height=3.5,width=6)
###Diagnostics plotting log likelihood
par(mfrow=c(1,2))
plot(FittingTable[,2],LL,xlab="ProtectiveTitre",ylab="log-likelihood")

###Plotting fit.
FittingTable$Protected_PNC<-FittingTable$RatioTitres>PNC_Threshold
EstimatedEfficacy_PNC=aggregate(Protected_PNC~Study,data=FittingTable,FUN=function(x)sum(x)/length(x))
PNC_PlottingTable<-join(SummaryTable_Efficacy_NeutRatio_SD_SEM[,c("Study","Efficacy")],EstimatedEfficacy_PNC,by="Study")
plot(PNC_PlottingTable$Efficacy,PNC_PlottingTable$Protected_PNC,xlab="Observed efficacy",ylab="Estimated efficacy")
abline(0,1,lty=2)

dev.off()




###### Manuscript Figures

NeutValue=seq(0.1,11,by=0.001)

SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_LB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)-1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)
SummaryTable_Efficacy_NeutRatio_SD_SEM$RatioReported_UB=10^((SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported)+1.96*SummaryTable_Efficacy_NeutRatio_SD_SEM$SEM)

Efficacy_Logistic<-NULL
Efficacy_Logistic_Raw<-NULL
for (i in 1:length(NeutValue)) {
  # Efficacy_Logistic[i]=LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],coef(FittedLogistic_MeanRept_SDMelb)[1],coef(FittedLogistic_MeanRept_SDMelb)[2])  
  Efficacy_Logistic_Raw[i]=LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1))  
}

LogisticModel_withPoolSD=data.frame("NeutRatio_Reported"=log10(NeutValue),"Efficacy"=Efficacy_Logistic_Raw)
LogisticModel_withPoolSD$Study<-rep("LogisticModel",length(NeutValue))

##### Adding Covaxin
CovaxinTable<-SummaryTable_Efficacy_NeutRatio_SD_SEM[1,]
CovaxinTable[1,]<-NA
CovaxinTable$Study<-"Covaxin"
CovaxinTable$TechnicalName<-"Covaxin"
CovaxinTable$Age<-"all"
CovaxinTable$NeutMean<-127.6773881
CovaxinTable$NeutConv<-161.1587484
CovaxinTable$NeutRatio_Reported<-log10(CovaxinTable$NeutMean/CovaxinTable$NeutConv)
CovaxinTable$RatioReported_LB<-CovaxinTable$NeutRatio_Reported
CovaxinTable$RatioReported_UB<-CovaxinTable$NeutRatio_Reported
CovaxinTable$Efficacy<-0.806
Covaxin_PointEstimate=LogisticModel_PercentUninfected(CovaxinTable$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1))  

### Confidnece bounds:
Cov<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]

Covaxin_temp <- function(p_temp) LogisticModel_PercentUninfected(CovaxinTable$NeutRatio_Reported,SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],p_temp[1],p_temp[2])
grad1_covaxin<-numericGradient(Covaxin_temp, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)))[1]
grad2_covaxin<-numericGradient(Covaxin_temp, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)))[2]
G_covaxin<-cbind(grad1_covaxin,grad2_covaxin)
Covaxin_LowerB=Covaxin_PointEstimate-1.96*sqrt(G_covaxin%*%Cov%*%t(G_covaxin))
Covaxin_UpperB=Covaxin_PointEstimate+1.96*sqrt(G_covaxin%*%Cov%*%t(G_covaxin))

#################################################### Adding 95% Prediction Intervals #########################################################
grad1<-NULL
grad2<-NULL
Lower_Pred<-NULL
Upper_Pred<-NULL
G<-NULL

for (i in 1:length(NeutValue)) 
{
  f_temp <- function(p_temp) LogisticModel_PercentUninfected(log10(NeutValue[i]),SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],p_temp[1],p_temp[2])
  grad1[i]<-numericGradient(f_temp, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)))[1]
  grad2[i]<-numericGradient(f_temp, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)))[2]
  G<-cbind(grad1[i],grad2[i])
  Lower_Pred[i]=Efficacy_Logistic_Raw[i]-1.96*sqrt(G%*%Cov%*%t(G))
  Upper_Pred[i]=Efficacy_Logistic_Raw[i]+1.96*sqrt(G%*%Cov%*%t(G))
}

LogisticModel_withPoolSD$Lower<-100*c(Lower_Pred)
LogisticModel_withPoolSD$Upper<-100*c(Upper_Pred) 
##########################################################################################################################################


Figure1<-ggplot(data=SummaryTable_Efficacy_NeutRatio_SD_SEM, aes(y=100*Efficacy,x=(10^NeutRatio_Reported),group=Study)) +
  #adding the bands
  geom_ribbon(data=LogisticModel_withPoolSD,aes(ymin=Lower, ymax=Upper, fill = "95% Prediction Intervals"), alpha = 0.3)+
  geom_point(shape=1) +
  geom_errorbar(aes(ymin=Lower,ymax=Upper)) +
  geom_errorbarh(aes(xmin=RatioReported_LB,xmax=RatioReported_UB)) +
  scale_x_log10(lim=c(0.1,11),breaks=c(0.125,0.25,0.5,1,2,4,8),labels=c(0.125,0.25,0.5,1,2,4,8)) +
  scale_y_continuous(lim=c(25,100)) +
  
  theme_linedraw() +
  geom_line(data=LogisticModel_withPoolSD,aes(color=Study)) +

  
  geom_point(data=CovaxinTable,aes(color=Study),size=4,shape=1) +
  geom_text(aes(label=TechnicalName),vjust=0,hjust=0, nudge_x=0.05, nudge_y = 0.0) +
  geom_text(data=CovaxinTable,aes(label=TechnicalName),vjust=0,hjust=0, nudge_x=0.05, nudge_y = 0.0) +
  xlab("Neutralisation titre (/convalescent plasma)") +
  ylab("Protective efficacy")

pdf("Figure1.pdf",height=5,width=7)
Figure1
dev.off()
Figure1




####Figure 1a - this creates EC50 estimates for each cohort 
# for Deborah's figure with distributions
fixlogk=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1]
FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model<-nlm(function(p){FittingLogistic_Raw(log(SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]/SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)]),fixlogk,p[(1:sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)))],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$NumVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfCont[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$InfVac[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],
                                                                                          SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)])},rep(LogisticEstimate[2],sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))),hessian=TRUE)

FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model_CI<- matrix(c(-tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model$hessian))),sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)))*1.96+tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont))),tail(sqrt(diag(solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model$hessian))),sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)))*1.96+tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,sum(!is.na(SummaryTable_Efficacy_NeutRatio_SD_SEM$NumCont)))),ncol=2)

FittingEachStudySep_LogisticRaw_ReptMean_PooledSD=data.frame("Study"=SummaryTable_Efficacy_NeutRatio_SD_SEM$Study,"EC50"=10^tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model$estimate,length(SummaryTable_Efficacy_NeutRatio_SD_SEM$Study)),"CI_L"=10^FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model_CI[,1],"CI_U"=10^FittedLogistic_RawEfficacy_MeanRept_SDPool_EC50Model_CI[,2])

write.csv(FittingEachStudySep_LogisticRaw_ReptMean_PooledSD,file="FittingEachStudySep_LogisticRaw_ReptMean_PooledSD.csv")

#####Protective Titres by study
SuppTable_TableOfSDperStudy<-read.csv("SuppTable_TableOfSDperStudy.csv")
SuppTable_TableOfSDperStudy$ProtectiveTitre<-SuppTable_TableOfSDperStudy$NeutConv*TableOfEstimatedTitres$Estimate_EC50[TableOfEstimatedTitres$Method=="Logistic_Raw" & TableOfEstimatedTitres$Mean=="Reported" & TableOfEstimatedTitres$SD=="Pooled_cens"]
SuppTable_TableOfSDperStudy$ProtectiveTitre<-round(SuppTable_TableOfSDperStudy$ProtectiveTitre,0)
write.csv(SuppTable_TableOfSDperStudy,"SuppTable_TableOfSDperStudy_withProtectiveTitres.csv")
