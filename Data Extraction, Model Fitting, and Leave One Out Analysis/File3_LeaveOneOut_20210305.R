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
SummaryTable_Efficacy_NeutRatio_SD_SEM$LOD_adj=log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$LOD)-log10(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutConv)
FittingTable<-IndividualNeutData_NormalisedbyConv[,c("Study","RatioTitres")]


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
##Fitting Models with Leave one out.

### Fitting Logistic Model to Different combinations of Mean and SD
LogisticEstimate=c("logk"=log(2.7),"C50"=log10(0.5))##InitialGuess

StudyList<-SummaryTable_Efficacy_NeutRatio_SD_SEM$Study

LeaveOneOut_LogisticTable<-data.frame("Study"=rep(StudyList,each=7),
                                      "TechnicalName"=rep(SummaryTable_Efficacy_NeutRatio_SD_SEM$TechnicalName,each=7),
                                      "Efficacy"=rep(SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy,each=7),
                                      "Lower"=rep(SummaryTable_Efficacy_NeutRatio_SD_SEM$Lower,each=7),
                                      "Upper"=rep(SummaryTable_Efficacy_NeutRatio_SD_SEM$Upper,each=7),
                                      "MeanMethod"=rep(c("CensoredEstimate","CensoredEstimate","CensoredEstimate","Reported","Reported","Reported","PNT"),length(StudyList)),
                                      "SDMethod"=rep(c("byStudy_cens","MelbStudy_cens","Pooled_cens","byStudy_cens","MelbStudy_cens","Pooled_cens","PNT"),length(StudyList)),
                                      "EC50"=NA,
                                      "k"=NA,
                                      "SD"=NA,
                                      "MeanNeut"=NA,
                                      "PredictedEfficacy"=NA,
                                      "Lower_Pred"=NA,
                                      "Upper_Pred"=NA)
indextemp<-0
for (i in 1:length(StudyList)) {
  
    tempSummaryTable<-SummaryTable_Efficacy_NeutRatio_SD_SEM[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study!=StudyList[i],]

    ### Raw Efficacy Data Models - Fitting Logistic Model to Different combinations of Mean and SD
    FittedLogistic_RawEfficacy_MeanCens_SDCenstemp<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(tempSummaryTable$NumCont))],p[sum(!is.na(tempSummaryTable$NumCont))+1],p[sum(!is.na(tempSummaryTable$NumCont))+2],tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$NumVac[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfVac[!is.na(tempSummaryTable$NumCont)],
                                                                tempSummaryTable$NeutRatio_cens[!is.na(tempSummaryTable$NumCont)],    tempSummaryTable$SD[!is.na(tempSummaryTable$NumCont)])},c(log(tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)]/tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)]),LogisticEstimate),hessian=TRUE)
      
      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^tail(FittedLogistic_RawEfficacy_MeanCens_SDCenstemp$estimate,1)
      LeaveOneOut_LogisticTable$k[indextemp]<-exp(tail(FittedLogistic_RawEfficacy_MeanCens_SDCenstemp$estimate,2)[1])
      LeaveOneOut_LogisticTable$MeanNeut[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$SD[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$SD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],SummaryTable_Efficacy_NeutRatio_SD_SEM$SD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],log(LeaveOneOut_LogisticTable$k[indextemp]),log10(LeaveOneOut_LogisticTable$EC50[indextemp]))
      
      
    FittedLogistic_RawEfficacy_MeanCens_SDMelbtemp<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(tempSummaryTable$NumCont))],p[sum(!is.na(tempSummaryTable$NumCont))+1],p[sum(!is.na(tempSummaryTable$NumCont))+2],tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$NumVac[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfVac[!is.na(tempSummaryTable$NumCont)],
                                                                tempSummaryTable$NeutRatio_cens[!is.na(tempSummaryTable$NumCont)],    tempSummaryTable$MelbSD[!is.na(tempSummaryTable$NumCont)])},c(log(tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)]/tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)]),LogisticEstimate),hessian=TRUE)
    
      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^tail(FittedLogistic_RawEfficacy_MeanCens_SDMelbtemp$estimate,1)
      LeaveOneOut_LogisticTable$k[indextemp]<-exp(tail(FittedLogistic_RawEfficacy_MeanCens_SDMelbtemp$estimate,2)[1])
      LeaveOneOut_LogisticTable$MeanNeut[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$SD[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],log(LeaveOneOut_LogisticTable$k[indextemp]),log10(LeaveOneOut_LogisticTable$EC50[indextemp]))
      
    FittedLogistic_RawEfficacy_MeanCens_SDPooltemp<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(tempSummaryTable$NumCont))],p[sum(!is.na(tempSummaryTable$NumCont))+1],p[sum(!is.na(tempSummaryTable$NumCont))+2],tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$NumVac[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfVac[!is.na(tempSummaryTable$NumCont)],
                                                                tempSummaryTable$NeutRatio_cens[!is.na(tempSummaryTable$NumCont)],    tempSummaryTable$PooledSD[!is.na(tempSummaryTable$NumCont)])},c(log(tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)]/tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)]),LogisticEstimate),hessian=TRUE)
    
      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^tail(FittedLogistic_RawEfficacy_MeanCens_SDPooltemp$estimate,1)
      LeaveOneOut_LogisticTable$k[indextemp]<-exp(tail(FittedLogistic_RawEfficacy_MeanCens_SDPooltemp$estimate,2)[1])
      LeaveOneOut_LogisticTable$MeanNeut[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$SD[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_cens[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],log(LeaveOneOut_LogisticTable$k[indextemp]),log10(LeaveOneOut_LogisticTable$EC50[indextemp]))
      
    FittedLogistic_RawEfficacy_MeanRept_SDCenstemp<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(tempSummaryTable$NumCont))],p[sum(!is.na(tempSummaryTable$NumCont))+1],p[sum(!is.na(tempSummaryTable$NumCont))+2],tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$NumVac[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfVac[!is.na(tempSummaryTable$NumCont)],
                                                                tempSummaryTable$NeutRatio_Reported[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$SD[!is.na(tempSummaryTable$NumCont)])},c(log(tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)]/tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)]),LogisticEstimate),hessian=TRUE)

      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^tail(FittedLogistic_RawEfficacy_MeanRept_SDCenstemp$estimate,1)
      LeaveOneOut_LogisticTable$k[indextemp]<-exp(tail(FittedLogistic_RawEfficacy_MeanRept_SDCenstemp$estimate,2)[1])
      LeaveOneOut_LogisticTable$MeanNeut[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$SD[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$SD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],SummaryTable_Efficacy_NeutRatio_SD_SEM$SD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],log(LeaveOneOut_LogisticTable$k[indextemp]),log10(LeaveOneOut_LogisticTable$EC50[indextemp]))
      
    FittedLogistic_RawEfficacy_MeanRept_SDMelbtemp<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(tempSummaryTable$NumCont))],p[sum(!is.na(tempSummaryTable$NumCont))+1],p[sum(!is.na(tempSummaryTable$NumCont))+2],tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$NumVac[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfVac[!is.na(tempSummaryTable$NumCont)],
                                                                tempSummaryTable$NeutRatio_Reported[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$MelbSD[!is.na(tempSummaryTable$NumCont)])},c(log(tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)]/tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)]),LogisticEstimate),hessian=TRUE)
    
      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^tail(FittedLogistic_RawEfficacy_MeanRept_SDMelbtemp$estimate,1)
      LeaveOneOut_LogisticTable$k[indextemp]<-exp(tail(FittedLogistic_RawEfficacy_MeanRept_SDMelbtemp$estimate,2)[1])
      LeaveOneOut_LogisticTable$MeanNeut[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$SD[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],log(LeaveOneOut_LogisticTable$k[indextemp]),log10(LeaveOneOut_LogisticTable$EC50[indextemp]))
      
    FittedLogistic_RawEfficacy_MeanRept_SDPooltemp<-nlm(function(p){FittingLogistic_Raw(p[1:sum(!is.na(tempSummaryTable$NumCont))],p[sum(!is.na(tempSummaryTable$NumCont))+1],p[sum(!is.na(tempSummaryTable$NumCont))+2],tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$NumVac[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$InfVac[!is.na(tempSummaryTable$NumCont)],
                                                                tempSummaryTable$NeutRatio_Reported[!is.na(tempSummaryTable$NumCont)],tempSummaryTable$PooledSD[!is.na(tempSummaryTable$NumCont)])},c(log(tempSummaryTable$InfCont[!is.na(tempSummaryTable$NumCont)]/tempSummaryTable$NumCont[!is.na(tempSummaryTable$NumCont)]),LogisticEstimate),hessian=TRUE)
    
      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^tail(FittedLogistic_RawEfficacy_MeanRept_SDMelbtemp$estimate,1)
      LeaveOneOut_LogisticTable$k[indextemp]<-exp(tail(FittedLogistic_RawEfficacy_MeanRept_SDMelbtemp$estimate,2)[1])
      LeaveOneOut_LogisticTable$MeanNeut[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$SD[indextemp]<-SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]]
      LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]<-LogisticModel_PercentUninfected(SummaryTable_Efficacy_NeutRatio_SD_SEM$NeutRatio_Reported[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==StudyList[i]],log(LeaveOneOut_LogisticTable$k[indextemp]),log10(LeaveOneOut_LogisticTable$EC50[indextemp]))
      
      ###Predictive interval for last Model
      Cov<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPooltemp$hessian)[8:9,8:9]
      
      predictiveFunc <- function(p_temp) LogisticModel_PercentUninfected(LeaveOneOut_LogisticTable$MeanNeut[indextemp],LeaveOneOut_LogisticTable$SD[indextemp],p_temp[1],p_temp[2])
      grad1<-numericGradient(predictiveFunc, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPooltemp$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPooltemp$estimate,1)))[1]
      grad2<-numericGradient(predictiveFunc, c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPooltemp$estimate,2)[1],tail(FittedLogistic_RawEfficacy_MeanRept_SDPooltemp$estimate,1)))[2]
      G<-cbind(grad1,grad2)
      LeaveOneOut_LogisticTable$Lower_Pred[indextemp]=LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]-1.96*sqrt(G%*%Cov%*%t(G))
      LeaveOneOut_LogisticTable$Upper_Pred[indextemp]=LeaveOneOut_LogisticTable$PredictedEfficacy[indextemp]+1.96*sqrt(G%*%Cov%*%t(G))
      
  FittingTabletemp<-FittingTable[FittingTable$Study!=StudyList[i],]
  PossibleThresholds=FittingTabletemp$RatioTitres
  LL<-NULL
  for (j in 1:length(PossibleThresholds)) {
    LL[j]=PNC_LikelihoodObservingAbove(PossibleThresholds[j],tempSummaryTable$Study,tempSummaryTable$Efficacy,FittingTabletemp,tempSummaryTable$LOD_adj)
  }
  
  ##Find Best Fit
  MaxLikelihood=max(LL)
  PNC_Threshold=FittingTabletemp[which.max(LL),2]
 
      indextemp<-indextemp+1
      LeaveOneOut_LogisticTable$EC50[indextemp]<-10^PNC_Threshold
      
      
    }

LeaveOneOut_LogisticTable$Model<-paste(LeaveOneOut_LogisticTable$MeanMethod,LeaveOneOut_LogisticTable$SDMethod)

write.csv(LeaveOneOut_LogisticTable,file="LeaveOneOut_LogisticTable.csv")

Figure1C<-ggplot(data=LeaveOneOut_LogisticTable[LeaveOneOut_LogisticTable$MeanMethod=="Reported" & LeaveOneOut_LogisticTable$SDMethod=="Pooled_cens",], aes(x=100*Efficacy,y=100*PredictedEfficacy)) +
  geom_point(shape=1) +
  geom_errorbar(aes(ymin=100*Lower_Pred,ymax=100*Upper_Pred)) +
  geom_errorbarh(aes(xmin=Lower,xmax=Upper)) +
  theme_linedraw() +
  geom_abline(intercept=0,slope=1) +
  scale_x_continuous(limits=c(30,102)) +
  scale_y_continuous(limits=c(30,102)) +
  geom_text(aes(label=TechnicalName),vjust=0,hjust=0, nudge_x=0.05, nudge_y = 0.0) +
  xlab("Observed Efficacy") +
  ylab("Predicted Efficacy")

pdf("Figure1C.pdf",height=5,width=7)
Figure1C
dev.off()
Figure1C

###Supplement of all model estimates
TableOfEstimatedTitres<-read.csv("TableOfEstimatedTitres.csv") ##For PNC mdoel

TableOfEstimatedTitres$Group<-paste(TableOfEstimatedTitres$Method,TableOfEstimatedTitres$Mean,TableOfEstimatedTitres$SD)
LeaveOneOut_LogisticTable$Group<-paste("Logistic_Raw",LeaveOneOut_LogisticTable$Model)
LeaveOneOut_LogisticTable$Group[LeaveOneOut_LogisticTable$MeanMethod=="PNT"]<-"PNT IndividualTitres IndividualTitres"
LeaveOneOut_LogisticTable$Estimate_EC50<-LeaveOneOut_LogisticTable$EC50

FigureS1<-ggplot(data=TableOfEstimatedTitres[TableOfEstimatedTitres$Method %in% c("Logistic_Raw","PNT"),c("Group","Estimate_EC50","L_CI","U_CI")], aes(y=Estimate_EC50,x=Group,group=Group)) +
  geom_point(aes(color=Group),size=2,shape=1) +
  scale_y_log10(lim=c(0.03,1)) +
  geom_errorbar(aes(ymin=L_CI,ymax=U_CI,color=Group),width=0.2) +
  theme_linedraw() +
  geom_point(data=LeaveOneOut_LogisticTable[,c("Group","Estimate_EC50")],
             aes(color=Group),
             size=1.2,
             shape=4,
             position = position_dodge2(width = 0.1)) +
  xlab("Model") +
  ylab("Protective titre (/convalescent plasma)")

pdf("FigureS1.pdf",height=3.5,width=6.5)
FigureS1
dev.off()
FigureS1



