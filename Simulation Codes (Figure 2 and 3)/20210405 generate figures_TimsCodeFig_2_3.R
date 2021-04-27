# reshape2 package is required to be installed

# Clear workspace, load libraries, set working directory, create figure directory
rm(list=ls())
library(ggplot2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir.create("figures")

# Function to tag figure files with date, and set figure dimensions
myfun_date <- function() format(Sys.time(), "%Y %m %d %H%M%S")
fig_height<-6
fig_width<-12

# Set parameter values
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
std <- log(10^std10) # Standard deviation in natural log units
hl <- 108 # Half life of antibody decay
dr <- -log(2)/hl # Corresponding decay rate in days for half life above
k <- 3.0/log(10) # logistic k value in natural log units (divided by log10 to convert it to natural logarithm units)
logk <- log(k) # log(k)
days <- 0:250 # Number of days to model for figure 2


# Define logistic model ----------------------------------------------------------------
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

# Figure 2A, efficacy vs time ------------------------------------------------------------
# Set initial efficacy level
ef <- c(99,95,90,80,70,60,50)

# Calculate standard normal deviate (of natural log distribution) to achieve that 
# efficacy by fixed threshold model, to use as initial guess for logistic model.
# This equates to number of standard errors between mean antibody level and threshold of protection
threshold <- qnorm(1-ef/100, mean=0, sd = std)

# Calculate standard normal deviation (of natural log distribution) with logistic model
threshold_logistic <- NULL
for (i in 1:length(threshold)){
  threshold_logistic[i] <- -nlm(function(mu){abs(LogisticModel_PercentUninfected(mu,std,logk,0)-ef[i]/100)},-threshold[i])$estimate
}

# Decay's the mean antibody level with time (of natural log distribution)
mn <- 0 + dr*days

# Calculates the efficacy with time and saves in Results
Results <- NULL
for (i in 1:length(threshold)) {
  # Efficacy with time using  a simple protection threshold model
  efficacy <- (1-pnorm(threshold[i], mean=mn, sd=std))*100
  # Efficacy with time using the logistic model
  efficacy2 <- NULL
  for (mn_index in mn) {
    efficacy2 <- c(efficacy2,LogisticModel_PercentUninfected(mn_index, std, logk, threshold_logistic[i])*100)
  }
  
  #Store results
  Results <- rbind(Results, data.frame(Efficacy = efficacy, Efficacy_logistic = efficacy2, initialEfficacy = ef[i], d=days))
}

# Re-order efficacy levels for ggplot
Results$initialEfficacy <- factor(Results$initialEfficacy, 
                                  levels = sort(unique(Results$initialEfficacy), decreasing=TRUE))

#Plot and sve figure of efficacy versus time
ggplot(Results,
       aes(x=d, y = Efficacy_logistic, color=initialEfficacy, group=initialEfficacy)) + 
  geom_line() + 
  theme_classic() +
  geom_hline(yintercept = ef, linetype="dotted", size=0.1) +
  xlab("Days") +
  ylab("Efficacy (%)") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,1),breaks=seq(0,100,10), limits=c(0,100)) +
  labs(color="Initial Efficacy")
ggsave(paste0("figures/Figure 2A logistic ",myfun_date(),".jpg"), width = fig_width, height=fig_height, units="cm")
ggsave(paste0("figures/Figure 2A logistic ",myfun_date(),".pdf"), width = fig_width, height=fig_height, units="cm")

# Calculates vaccine starting with initial efficacy of 95 and 70 to maintain what efficacy by 250 days?
subset(Results, d %in% c(0,250) & initialEfficacy == 95)
subset(Results, d %in% c(0,250) & initialEfficacy == 70)

# Figure 2B, Starting efficacy vs time to %X% efficacy ------------------------------------------------------------
# Returns the times to achieve some final efficacy
myfun_calculate_points <- function(fef) {
  # All possible starting efficacies
  sef <- seq(fef,99,0.1)
  
  # Difference between starting and finishing efficacy z scores - threshold model to use as initial guess
  zef <- qnorm(1-sef/100, mean=0, std)
  
  # Difference between starting and finishing efficacy standard normal deviate of antibody level (z scores)
  zef_logistic <- NULL
  for (i in 1:length(sef)){
    zef_logistic[i] <- -nlm(function(mu){abs(LogisticModel_PercentUninfected(mu,std,logk,0)-sef[i]/100)},-zef[i])$estimate
  }
  
  zef <- zef-zef[1]
  zef_logistic <- zef_logistic-zef_logistic[1]
  
  # Calculates the time required for mean to move that z score distance
  times <- zef/dr
  times_logistic <- zef_logistic/dr

  #return(data.frame(fef=fef, sef = sef, times = times))
  return(data.frame(fef=fef, sef = sef, times = times_logistic))
}

Results <- NULL
fef <- c(50,70) # Final efficacy levels

for (i in fef){ # Calculate times required to decay to final efficacy
  Results <- rbind(Results, myfun_calculate_points(i))
}
# Re-levels the final efficacy to be in order for plotting legend
Results$fef <- factor(Results$fef, 
              levels = sort(unique(Results$fef), decreasing=TRUE))

# Dataframe for geom_ribbon shading in ggplot that follows
ribbon_data <- data.frame(times = 0:250)
ribbon_data$seventy <- approx(Results[Results$fef==70,]$times, Results[Results$fef==70,]$sef, ribbon_data$times)$y
ribbon_data$fifty <- approx(Results[Results$fef==50,]$times, Results[Results$fef==50,]$sef, ribbon_data$times)$y
ribbon_data$top <- 100

# Plots initial efficacy versus time until a booster is required (days)
ggplot(subset(Results, times<=250),
       aes(x=times, y= sef, color= fef, group=fef)) +
  geom_ribbon(data=ribbon_data, aes(x=times, ymin=fifty, ymax=seventy), inherit.aes = FALSE, fill="blue", alpha=0.4)+
  geom_ribbon(data=ribbon_data, aes(x=times, ymin=seventy, ymax=top), inherit.aes = FALSE, fill="red", alpha=0.4)+
  geom_line(size=1) + 
  theme_classic() +
  geom_hline(yintercept=ef, linetype = "dotted", size=0.1)+
  geom_vline(xintercept=seq(50,250,50), linetype = "dotted", size=0.1)+
  ylab("Initial Efficacy (%)") +
  xlab("Time until booster required (days)") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("red", "blue"),name = "Time to", labels = c("70% efficacy", "50% efficacy"),
                     guide = guide_legend(reverse = TRUE)) 
ggsave(paste("figures/Figure 2B Logistic model ",myfun_date(),".jpg"), width = fig_width, height=fig_height, units="cm")
ggsave(paste("figures/Figure 2B Logistic model ",myfun_date(),".pdf"), width = fig_width, height=fig_height, units="cm")

# How long does it take to reach 70% or 50% efficacy for a particular starting efficacy
subset(Results, sef == 90) # starting efficacy of 90%

# Figure 2C, Variant vs WT protection  ------------------------------------------------------------
myfun_calculate_mt_efficacies <- function(fold_change){
  # set wild type efficacy
  wt_ef <- seq(50,99,0.5)
  
  # Standard normal deviate of antibody level to achieve starting efficacy
  # Using a simple threshold model as an initial guess for logistic model
  wt_thresholds <- qnorm(1-wt_ef/100, mean=0, sd=std) 
  
  # Equivalent thresholds in the logistic model
  wt_thresholds_logistic <- NULL
  for (i in 1:length(wt_ef)){
    wt_thresholds_logistic[i] <- -nlm(function(mu){abs(LogisticModel_PercentUninfected(mu,std,logk,0)-wt_ef[i]/100)}
                                      ,-wt_thresholds[i])$estimate
  }
  
  # Calculates the corresponding standard normal deviate (threshold) for a variant
  # given some fold reduction in mean antibody level
  mt_thresholds <- wt_thresholds + log(fold_change) # Simple threshold model
  mt_thresholds_logistic <- wt_thresholds_logistic + log(fold_change) # Logistic model
  
  # Calculate efficacy of variant
    mt_ef_logistic <- NULL
  for (i in 1:length(mt_thresholds_logistic)) {
    mt_ef_logistic[i] <- LogisticModel_PercentUninfected(0,std,logk,mt_thresholds_logistic[i])*100
  }
  
  return(data.frame(fold_change=fold_change, wt_ef = wt_ef, mt_ef_logistic))
}

Results <- NULL
fold_change <- c(2,5,10)
for (i in fold_change){
  Results <- rbind(Results, myfun_calculate_mt_efficacies(i))
}

# Formatting changes for figure
Results$fold_change <- paste0(Results$fold_change,"-fold")
Results$fold_change <- factor(Results$fold_change, levels=c("2-fold","5-fold","10-fold"))

# Plots the efficacy of variant versus efficacy of wild-type
ggplot(Results,
       aes(x=wt_ef, y=mt_ef_logistic, color=fold_change)) +
  geom_line(size=1) +
  xlab("Efficacy against wild-type (%)") +
  ylab("Efficacy against variant (%)") +
  labs(color="Decrease of\nbinding to variant:") +
  geom_vline(xintercept=ef, linetype = "dotted", size=0.1)+
  geom_hline(yintercept=c(seq(0,90,10),99), linetype = "dotted", size=0.1)+
  theme_classic() + 
  scale_x_continuous(limits=c(50,100),expand = c(0, 0)) +
  scale_y_continuous(limits=c(0,100),expand = c(0, 0), breaks=seq(0,100,10)) +
  geom_abline(slope=1, linetype = "longdash")
ggsave(paste0("figures/Figure 2C Logistic model",myfun_date(),".jpg"), width = fig_width, height=fig_height, units="cm")
ggsave(paste0("figures/Figure 2C Logistic model",myfun_date(),".pdf"), width = fig_width, height=fig_height, units="cm")

# Reduction in efficacy with fold lower titre
subset(Results, wt_ef==95)
subset(Results, wt_ef==70)


# Figure 3 Other parameters nor previously defined --------------------------------------------------------
hl1 <- hl # Early half life in days
dr1 <- -log(2)/hl1 # Early corresponding decay rate (days)
hl2 <- 3650 # Late half life in days
dr2 <- -log(2)/hl2 # Late corresponding decay rate
days <- seq(0,1000,1) # Total number of days to model decay with time
threshold_mild <- 0.20 # Mild EC50 value
threshold_severe <- 0.030 # Severe EC50 value
threshold_severe_lower <- 0.0071 # Severe EC50 value lower confidence limit
threshold_severe_upper <- 0.13 # Severe EC50 value upper confidence limit
threshold_dif <- (log(threshold_mild) - log(threshold_severe))/std # EC50 difference in natural log standard errors

# Figure 3C --------------------------------------------------------------------------------------------------
mild_ef <- c(0.001,seq(0.01,0.99,0.01)) # All efficacy ranges for mild infection
severe_threshold_reducation <- threshold_mild/c(threshold_severe_upper,threshold_severe, threshold_severe_lower) # fold reduction in threshold including the confidence interval
Results <- NULL

# Calculates the corresponding treshold for severe infection for all mild infection efficacy
for (sthr in severe_threshold_reducation){
  # reduction in mild threshold to give severe threshold
  td <- (log(threshold_mild) - log(1/sthr * threshold_mild))/std
  
  mild_threshold <- qnorm(1-mild_ef,0,1) # mild efficacy threshold as initial guess for logistic model
  
  # Threshold for the mild efficacy using the logistic model
  mild_threshold_logistic <- NULL
  for (i in 1:length(mild_ef)){
    mild_threshold_logistic[i] <- -nlm(function(mu){abs(LogisticModel_PercentUninfected(mu,std,logk,0)-mild_ef[i])}
                                         ,-mild_threshold[i])$estimate
  }
  
  # Calculates the severe efficacy corresponding to the mild efficacy
  severe_ef_logistic <- NULL
  for (i in 1:length(mild_ef)){
    severe_ef_logistic[i] <- LogisticModel_PercentUninfected(0, std, logk, mild_threshold_logistic[i]-td)
  }
  
  Results <- rbind(Results, data.frame(mild = mild_ef*100
                                       , severe = severe_ef_logistic*100
                                       , severe_reduction = sthr))
}
Results2 <- reshape2::dcast(mild~severe_reduction, data=Results, value.var="severe")
names(Results2)[2:4] <- c("severe_lb","severe_estimate","severe_ub")

# plot efficacy against severe infection versus efficacy against any (mild) infection
ggplot(Results2,
       aes(x=mild, y=severe_estimate, ymin=severe_lb, ymax=severe_ub)) +
  geom_ribbon(fill="blue",alpha=0.5) +
  geom_line() +
  xlab("Efficacy against any infection (%)") +
  ylab("Efficacy against \nsevere infection (%)") +
  theme_classic() +
  geom_hline(yintercept = c(seq(50,90,10), 95), linetype="dotted", color="#AAAAAA") +
  geom_vline(xintercept = c(seq(50,90,10), 95), linetype="dotted", color="#AAAAAA") +
  scale_x_continuous(expand=c(0,0), limits=c(0,101)) + 
  scale_y_continuous(expand=c(0,0), limits=c(0,100))
ggsave(paste0("figures/Figure 3C Logistic model ",myfun_date(),".jpg"), width = fig_width, height=fig_height, units="cm")
ggsave(paste0("figures/Figure 3C Logistic model ",myfun_date(),".pdf"), width = fig_width, height=fig_height, units="cm")

# Figure 3D --------------------------------------------------------------------------------------------------
startingTitre <- c(1/4, 1/2, 1, 2, 4)
targetTimes<-365*c(1,1.5,2)

myfun_titre <- function(x,targettime, ct=250) {
  # x is the starting titre
  # ct is the cut of time after which the decay rate starts to reduce (default 250 days)
  
  decayRate <- rep(dr1, length(days))
  decayRate[days>targettime]=dr2
  
  slowing=(1/(targettime-ct))*(dr1-dr2)
  
  decayRate[days>ct & days<=targettime] <- dr1-slowing*(days[days>ct & days<=targettime]-ct)
  
  titre_with_time <- log(x)
  for (i in 2:length(days)){
    titre_with_time[i] <- titre_with_time[i-1]+decayRate[i]
  }
  
  return(exp(titre_with_time))
}
Results <- NULL
for (i in 1:length(startingTitre)){
  for (j in 1:length(targetTimes)) {
  st <- startingTitre[i]
  TargTime<-targetTimes[j]
  Results <- rbind(Results, data.frame(start_titre = st,
                                       target_time = TargTime,
                                       Days = days,
                                       titre_time = myfun_titre(st,TargTime)))
  }
}
Results$start_titre <- factor(Results$start_titre)
Results$target_time <- factor(Results$target_time)

Breaks_power <- seq(-6,3,1)
Breaks <- 2^Breaks_power
Breaks_label <- NULL
for (i in 1:length(Breaks_power)) {
  if (Breaks_power[i]>=0) {
    Breaks_label[i] <- as.character(2^Breaks_power[i])
  } else Breaks_label[i] <- paste0("1/",2^-(Breaks_power[i]))
  
}

Results$tt <- Results$Days<250

ggplot(Results[Results$start_titre %in% c(0.25,1,4),],
       aes(x=Days, y = titre_time, color=target_time, group=interaction(start_titre,target_time,tt), linetype=tt)) +
  geom_line() +
  geom_hline(yintercept = threshold_mild, color="green")+
  geom_hline(yintercept = threshold_severe, color="blue")+
  scale_y_log10(expand=c(0,0), breaks = Breaks
                #, labels=Breaks_label
                , labels=as.character(round(Breaks,3))
  ) +
  scale_x_continuous(expand=c(0,0), limits=c(min(days),max(days)+20))+
  annotate("rect",xmin=min(Results$Days), xmax=max(Results$Days),
           ymin=0.01, ymax=threshold_severe, alpha=0.1, fill="blue") +
  annotate("rect",xmin=min(Results$Days), xmax=max(Results$Days),
           ymin=threshold_severe, ymax=threshold_mild, alpha=0.1, fill="green") +
  geom_vline(xintercept=250, linetype="dotted", color="#AAAAAA", alpha=0.4) +
  geom_hline(yintercept=Breaks, linetype="dotted", color="#AAAAAA", alpha=0.4) +
  scale_linetype_manual(values=c("longdash","solid"))+
  ylab("Neutralisation level\n(fold of convalescent)") +
  xlab("Time (days)") +
  theme_classic() + theme(legend.position = "none")
ggsave(paste0("figures/Figure 3D ",myfun_date(),".jpg"), width = fig_width, height=fig_height, units="cm")
ggsave(paste0("figures/Figure 3D ",myfun_date(),".pdf"), width = fig_width, height=fig_height, units="cm")


# Figure 3E ------------------------------------------------------------------------------------
mild_ef <- c(70,80,90,95)
targetTimes<-365*c(1,1.5,2)

q_neut_titre_mild <- qnorm(1-mild_ef/100, 0, std) # q value for mild threshold in natural log units
q_neut_titre_mild_logistic <- NULL
for (i in 1:length(mild_ef)){
  q_neut_titre_mild_logistic[i] <- -nlm(function(mu){abs(LogisticModel_PercentUninfected(mu,std,logk,0)-mild_ef[i]/100)},-q_neut_titre_mild[i])$estimate
}

q_neut_titre_severe <- q_neut_titre_mild - threshold_dif
q_neut_titre_severe_logistic <- q_neut_titre_mild_logistic - threshold_dif

# in natural log units


Results <- NULL
for (i in 1:length(mild_ef)){
  for (j in 1:length(targetTimes)) {
  mef <- mild_ef[i]
  tragTime<-targetTimes[j]
  
  decay_of_mean_w_time <- log(myfun_titre(1,tragTime)) 
  
  # mildEfficacy <- 100*(1-pnorm(q_neut_titre_mild[i],decay_of_mean_w_time, std))
  # severeEfficacy <- 100*(1-pnorm(q_neut_titre_severe[i],decay_of_mean_w_time, std))
  # 
  mildEfficacyLogistic <- NULL
  for (j in 1:length(decay_of_mean_w_time)){
    mildEfficacyLogistic[j] <- LogisticModel_PercentUninfected(decay_of_mean_w_time[j], std, logk, q_neut_titre_mild_logistic[i])*100
  }
  
  severeEfficacyLogistic <- NULL
  for (j in 1:length(decay_of_mean_w_time)){
    severeEfficacyLogistic[j] <- LogisticModel_PercentUninfected(decay_of_mean_w_time[j], std, logk, q_neut_titre_severe_logistic[i])*100
  }
  
  Results <- rbind(Results, data.frame(mild_start_efficacy = mef,
                                       target_time = tragTime/365,
                                       Days = days,
                                       mild_efficacy = mildEfficacyLogistic,
                                       sever_efficacy = severeEfficacyLogistic
  ))
  }
}

Results2 <- reshape2::melt(Results, id.vars = c("mild_start_efficacy","target_time","Days"))
Results2$mild_start_efficacy <- factor(Results2$mild_start_efficacy)
Results2$target_time <- as.factor(Results2$target_time)

ggplot(Results2,
       aes(x=Days, y = value
           , color=target_time
           , group=interaction(variable,mild_start_efficacy,target_time)
           , linetype = variable)) + 
  geom_line() +
  ylab("Efficacy (%)") +
  xlab("Time (days)") +
  scale_linetype_manual(values=c("solid", "dashed"),name = "Severity", labels = c("Mild", "Severe")) +
  scale_color_discrete(name = "Time to\nplateau (y)") +
  scale_y_continuous(limits=c(0,100),breaks=c(0,20,40,60,80,100) , expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0))+
  theme_classic() +
  facet_wrap("mild_start_efficacy", scales='free') + 
  theme(strip.background = element_blank(), strip.text = element_blank())
ggsave(paste0("figures/Figure 3E Logistic",myfun_date(),".jpg"), width = fig_width, height=fig_height+1, units="cm")
ggsave(paste0("figures/Figure 3E Logistic",myfun_date(),".pdf"), width = fig_width, height=fig_height+1, units="cm")
