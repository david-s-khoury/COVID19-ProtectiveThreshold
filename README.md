# COVID19-ProtectiveThreshold
R codes for the article "Neutralizing antibody levels are highly predictive of immune protection from symptomatic SARS-CoV-2 infection"
There are four folders:
1. Comparing Severe vs Mild Infection
2. Data Extraction, Model Fitting, and Leave One Out Analysis
3. Estimating Decay Rates
4. Simulation Codes (Figure 2 and 3)
The scripts run independently from each folder.

Comparing Severe vs Mild Infection Folder
The R file "Fitting Mild and Severe cases.R" is the main script to implement the model comparison between mild vs severe cases. 
The basic idea is, whether both mild and severe cases can be fitted by a single parameter(s), and need separate parameter(s).
There are two data files here, "SummaryTable_Efficacy_NeutRatio_SD_SEM.csv" which stores the efficacy and titre data for mild infection; 
and "SummaryTable_Efficacy_NeutRatio_SD_SEM_Severe.csv" which contains efficacy and titre data for severe infection.

Data Extraction, Model Fitting, and Leave One Out Analysis
The R files need to be run in this order: "File1_FitExtractedDataForMeansAndStDev_20210402.R", "File2_FittingModels_20210404.R"; and "File3_LeaveOneOut_20210305.R".
The first script creates an estimate of the meand and SD for each study by fitting a censored normal distrbution based on the raw data. This will create an output file 
named "SummaryTable_Efficacy_NeutRatio_SD_SEM.csv", which is an input for the second script. Model fitting was done in the second script, to fit a logistic model and 
a binary classifier model. Finally, a leave one out validation analysis can be found the the third script.

Estimating Decay Rates
There are two sub-folders in this folder; "Convalescent (Dan et al)" and "Vaccine (Widge et al) vs Convalescent (Wheatley et al)". Inside each folder, there is an independent
dataset and R script to fit a censored linear mixed effect model to estimate the decay slope of the neutralisation.

Simulation Codes (Figure 2 and 3)
The main R script "20210511 generate figures_TimsCodeFig_2_3.R" is needed to generate Figure 2 and 3 from the main article.
