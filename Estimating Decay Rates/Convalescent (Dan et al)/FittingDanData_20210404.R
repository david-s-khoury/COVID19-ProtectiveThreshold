library('ggplot2')
library('lmec')
library('plyr')

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

dataset <- read.csv("Dan2021_abf4063_Data_S1_ForFitting_20210402.csv")
dataset$Neut<-as.numeric(dataset$Neut)
dataset$cens<-dataset$Neut<=dataset$LOD
dataset$DonorID<-as.factor(dataset$DonorID)

ListOfTimePoints=aggregate(dataset[!is.na(dataset$Neut),c("DonorID")],by=list(dataset[!is.na(dataset$Neut),c("DonorID")]),FUN=length)

# subdataset<-dataset[dataset$DonorID %in% ListOfTimePoints$Group.1[ListOfTimePoints$x>1],]
# 
###Remove 4844 and 4843 because repeated measures on same day
dataset<-dataset[dataset$DonorID!="4844" & dataset$DonorID!="4843",]

###### Fitting LMEC

ylog = log(dataset$Neut[!is.na(dataset$Neut)]) #LogN Transformed data
cens =dataset$cens[!is.na(dataset$Neut)] #0 if real, 1 if it's censored (ie, 0 values)
cluster = match(dataset$DonorID[!is.na(dataset$Neut)],unique(dataset$DonorID[!is.na(dataset$Neut)])) #ID for each patient

id = unique(dataset$DonorID[!is.na(dataset$Neut)])
day = dataset$Days[!is.na(dataset$Neut)] #day 0 is day 26 in my fitting (the first sampling date)


#Initialize design matrix for the random effects
#A is for the random intercept
A <- matrix(1, nrow=NROW(cluster), ncol=1)

#B is for the random slope
M<-A*day #to tell lmec that we need random slopes
B <- matrix(NA, nrow=NROW(cluster), ncol=2)
B[, 1] <- A
B[, 2] <- M



#X is the fixed design matrix, it's the first column is for the intercept, second column for the slope, the third one is for the replicate fixed effect
intercept_matrix = matrix(rep(1, length(ylog)), ncol=1) #1 for the intercept

#Single Phase
XSingle = cbind(intercept_matrix,day)
#X = cbind(intercept_matrix,day,elements$Replicate)
fitSingle = lmec(yL=ylog,cens=cens, X=XSingle, Z=B, cluster=cluster, method='ML',maxstep =10000)

x=seq(25,240,by=0.1)
y=fitSingle$beta[1]+x*fitSingle$beta[2]
ylinear=exp(y)
fittedmodeldata=data.frame(x,ylinear,DonorID='none')

FigureS3<-ggplot(data=dataset,aes(x=Days,y=Neut,color=DonorID)) + geom_point() +geom_line() +
  geom_line(data=fittedmodeldata,aes(x=x,y=ylinear),color='black') +
  scale_y_log10() +
  theme_linedraw() +
  theme(legend.position='none', panel.background = element_blank()) +
  xlab("Days") +
  ylab("Neutralisation Titre")

pdf("FigureS3.pdf",height=3.5,width=5)
FigureS3
dev.off()
FigureS3


halflife=-log(2)/fitSingle$beta[2]
CIs=-log(2)/(fitSingle$beta[2]+c(-1,1)*1.96*sqrt(diag(fitSingle$varFix))[2])
