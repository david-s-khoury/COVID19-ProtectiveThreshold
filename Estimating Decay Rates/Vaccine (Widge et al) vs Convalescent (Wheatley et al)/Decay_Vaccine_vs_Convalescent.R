rm(list=ls(all=TRUE))
setwd("C:/Projects/Covid Decay with Vaccine")

dev.off()
library(nlme)
library(pracma)
library(lmec)

elements<-read.csv("Data_Convalescent_vs_Vaccine.csv", fileEncoding="UTF-8-BOM") ##For threshold and logistic

cens = elements$cens_a #censoring indicator, put 1 if the value is below threshod
y = log(elements$y) #dependent variable to fit (log of neutralization)
cluster = elements$cluster
id = unique(cluster)
day = elements$day #day26 (the first sampling time) is day 0 in the fitting
intercept_matrix = matrix(rep(1, length(y)), ncol=1)

#X is the design matrix for the fixed effect; group is a binary variable for convalescent vs vaccine data
X = cbind(intercept_matrix,day,elements$group,elements$group*day)


#B is the design matrix for the random effects, Have random effect for intercept and slope
B = cbind(intercept_matrix,day)

#call lmec function
fit_linear_different_rate = lmec(yL=(y),cens=(cens), X=X, Z=B, cluster=as.numeric(cluster), method='ML',varstruct = 'unstructured')


#X is the design matrix for the fixed effect; difference only in the intercept
X = cbind(intercept_matrix,day,elements$group)

#B is the design matrix for the random effects, Have random effect for intercept and slope
B = cbind(intercept_matrix,day)

#call lmec function
fit_linear_same_rate = lmec(yL=(y),cens=(cens), X=X, Z=B, cluster=as.numeric(cluster), method='ML',varstruct = 'unstructured')
