library('plyr')



#### Set directory to same directory as the r-script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### Import Tables Required for Estimating Mean and SD for each study
IndividualNeutData_NormalisedbyConv<-read.csv("IndividualNeutData_NormalisedbyConv.csv")
IndividualNeutData_RawVaccineAndConv<-read.csv("IndividualNeutData_RawVaccineAndConv.csv")
EfficacyTable<-read.csv("EfficacyMeanTable_20210226.csv")


####Test for difference in variance between vaccine groups
fligner.test(TitreLog ~ Study, data = IndividualNeutData_RawVaccineAndConv[IndividualNeutData_RawVaccineAndConv$Group=="Vaccine",])

###################################################################################################
#########Fitting Normal Distribution to Neut data with 

###Creating Table of all the estimated Means and SD after censoring for each study
#These will be used in model fitting and to estimate SEM

#negative log likelihood of normal distibrution model with censoring
Likelihood=function(p,log10_censT,data){-sum(log(dnorm(data[data>log10_censT],p[1],p[2])))-sum(log(pnorm(data[data<=log10_censT],p[1],p[2])))}

ListOfEstimatedMeanSDafterCensoring<-unique(IndividualNeutData_RawVaccineAndConv[,c("Study","Group")])

ListOfEstimatedMeanSDafterCensoring$SD<-NA
ListOfEstimatedMeanSDafterCensoring$EstimateMean<-NA
ListOfEstimatedMeanSDafterCensoring$NumberIndividuals<-NA
for (i in 1:nrow(ListOfEstimatedMeanSDafterCensoring)){
  
  tempdata=IndividualNeutData_RawVaccineAndConv$TitreLog[IndividualNeutData_RawVaccineAndConv$Study==ListOfEstimatedMeanSDafterCensoring$Study[i] & IndividualNeutData_RawVaccineAndConv$Group==ListOfEstimatedMeanSDafterCensoring$Group[i]]
  
  fitmdltemp<-nlm(function(p){Likelihood(p,log10(EfficacyTable$LOD[EfficacyTable$Study==ListOfEstimatedMeanSDafterCensoring$Study[i]]),tempdata)},c(mean(tempdata),sd(tempdata)))    


  ###Creating Table of all the estimated Means and SD after censoring for each study
  #These will be used in model fitting and to estimate SEM
  ListOfEstimatedMeanSDafterCensoring$SD[i]<-fitmdltemp$estimate[2]
  ListOfEstimatedMeanSDafterCensoring$EstimateMean[i]<-fitmdltemp$estimate[1]
  ListOfEstimatedMeanSDafterCensoring$NumberIndividuals[i]<-length(tempdata)
}



#### Next we wish to pool all the studies together (and still consider censoring) 
# and determine a pooled SD
ListLODtemp<-EfficacyTable$LOD[match(IndividualNeutData_NormalisedbyConv$Study,EfficacyTable$Study)]
ListLODtemp<-log10(ListLODtemp)-log10(IndividualNeutData_NormalisedbyConv$NeutConv)-IndividualNeutData_NormalisedbyConv$MeanRatio
names(ListLODtemp)<-NULL
PooledSDModelFit<-nlm(function(p){Likelihood(p,log10(ListLODtemp),IndividualNeutData_NormalisedbyConv$CentredRatio)},c(mean(IndividualNeutData_NormalisedbyConv$CentredRatio),sd(IndividualNeutData_NormalisedbyConv$CentredRatio)))    
PooledSD<-PooledSDModelFit$estimate[2]




### For the convalescence study in Melb - create duplicate row called 
#Conv to match Vaccine because for this group the ratio neut between 
# conv and vaccine is 1
temprow<-ListOfEstimatedMeanSDafterCensoring[ListOfEstimatedMeanSDafterCensoring$Study=="Convalescence",]
temprow$Group="Conv"
ListOfEstimatedMeanSDafterCensoring<-rbind(ListOfEstimatedMeanSDafterCensoring,temprow)


####Creating a Table where we will calculate ratio of neut vaccine vs conv
# for each study using these means from censored normal distribution fit
RowIndex<-ListOfEstimatedMeanSDafterCensoring[ListOfEstimatedMeanSDafterCensoring$Group!="Conv",c("Study","Group","SD")]
CalculatingRatioandSEMafterCensoring<-RowIndex
CalculatingRatioandSEMafterCensoring$NeutRatio_cens<-NA

for (i in 1:nrow(RowIndex)) {
  
  #Calculate Ratio Vaccine to Conv Neut
  CalculatingRatioandSEMafterCensoring$NeutVaccine_cens[i]<-ListOfEstimatedMeanSDafterCensoring$EstimateMean[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group==RowIndex$Group[i] ]
  CalculatingRatioandSEMafterCensoring$NeutConv_cens[i]<-ListOfEstimatedMeanSDafterCensoring$EstimateMean[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group=="Conv" ]
  CalculatingRatioandSEMafterCensoring$NeutRatio_cens[i]<-ListOfEstimatedMeanSDafterCensoring$EstimateMean[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group==RowIndex$Group[i] ]-ListOfEstimatedMeanSDafterCensoring$EstimateMean[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group=="Conv"]
  #Number of Individuals used in estimate of Conv SD
  CalculatingRatioandSEMafterCensoring$NumberIndividuals_Conv[i]=ListOfEstimatedMeanSDafterCensoring$NumberIndividuals[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group=="Conv"]
  #Number of Individuals used in estimate of Vaccine SD
  CalculatingRatioandSEMafterCensoring$NumberIndividuals_Vaccine[i]=ListOfEstimatedMeanSDafterCensoring$NumberIndividuals[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group==RowIndex$Group[i] ]
  #SEM of the ratio of neut from vaccine to Conv
  CalculatingRatioandSEMafterCensoring$SEM[i]<-sqrt(((ListOfEstimatedMeanSDafterCensoring$SD[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group==RowIndex$Group[i] ]^2)/CalculatingRatioandSEMafterCensoring$NumberIndividuals_Vaccine[i])+((ListOfEstimatedMeanSDafterCensoring$SD[ListOfEstimatedMeanSDafterCensoring$Study==RowIndex$Study[i] & ListOfEstimatedMeanSDafterCensoring$Group=="Conv"]^2)/CalculatingRatioandSEMafterCensoring$NumberIndividuals_Conv[i]))
  
  }



##### Building Summary Table with all SD and Means and Efficacy's
# This table will be used for most model fitting and figure creation
SummaryTable_Efficacy_NeutRatio_SD_SEM<-join(EfficacyTable,CalculatingRatioandSEMafterCensoring,by="Study")
SummaryTable_Efficacy_NeutRatio_SD_SEM$MelbSD<-CalculatingRatioandSEMafterCensoring$SD[CalculatingRatioandSEMafterCensoring$Study=="Convalescence"]
SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD<-PooledSD
  
write.csv(SummaryTable_Efficacy_NeutRatio_SD_SEM,"SummaryTable_Efficacy_NeutRatio_SD_SEM.csv")



####Export SD Estimates for each study to Supp. Table

SuppTable_TableOfSDperStudy<-join(SummaryTable_Efficacy_NeutRatio_SD_SEM[,c("Study","TechnicalName","LOD","NeutMean","NeutConv")],ListOfEstimatedMeanSDafterCensoring[ListOfEstimatedMeanSDafterCensoring$Group=="Vaccine",c("Study","NumberIndividuals","SD")],by="Study")
colnames(SuppTable_TableOfSDperStudy)[(ncol(SuppTable_TableOfSDperStudy)-1):ncol(SuppTable_TableOfSDperStudy)]<-c("NumberIndividualsVaccine","SDVaccine")
SuppTable_TableOfSDperStudy<-join(SuppTable_TableOfSDperStudy,ListOfEstimatedMeanSDafterCensoring[ListOfEstimatedMeanSDafterCensoring$Group=="Conv",c("Study","NumberIndividuals","SD")],by="Study")
colnames(SuppTable_TableOfSDperStudy)[(ncol(SuppTable_TableOfSDperStudy)-1):ncol(SuppTable_TableOfSDperStudy)]<-c("NumberIndividualsConv","SDConv")
PooledRow<-data.frame("Study"="Pooled*","TechnicalName"="Pooled*","NeutMean"=NA,"NeutConv"=NA,"LOD"=NA,"NumberIndividualsVaccine"=sum(ListOfEstimatedMeanSDafterCensoring$NumberIndividuals[ListOfEstimatedMeanSDafterCensoring$Group=="Vaccine"]),"SDVaccine"=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1],"NumberIndividualsConv"=NA,"SDConv"=NA)
SuppTable_TableOfSDperStudy<-rbind(SuppTable_TableOfSDperStudy,PooledRow)
SuppTable_TableOfSDperStudy[SuppTable_TableOfSDperStudy$Study=="Convalescence",c("NumberIndividualsVaccine","SDVaccine")]=NA
SuppTable_TableOfSDperStudy$SDVaccine<-round(SuppTable_TableOfSDperStudy$SDVaccine,2)
SuppTable_TableOfSDperStudy$SDConv<-round(SuppTable_TableOfSDperStudy$SDConv,2)
SuppTable_TableOfSDperStudy$LOD[SuppTable_TableOfSDperStudy$LOD==0]=NA

write.csv(SuppTable_TableOfSDperStudy,"SuppTable_TableOfSDperStudy.csv")



###Normality of combined distributions
#Test Normality of combined data
shapiro.test(CombinedTable_RawNeuts_LOD_centred$CentredTitre[CombinedTable_RawNeuts_LOD_centred$Group=="Vaccine" & CombinedTable_RawNeuts_LOD_centred$TitreLog>log10(CombinedTable_RawNeuts_LOD_centred$LOD)])


