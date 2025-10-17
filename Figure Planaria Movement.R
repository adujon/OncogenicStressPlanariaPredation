
rm(list=ls()) #remove all objects into memmory
library(brms)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(job)
library(bayesplot)
library(ggplot2)
library(survival)
library(survminer)
library(ggpubr)
library(tidybayes)
library(ggsurvfit)

#place raw dataset in the same folder as the R script file
#place the cura_binom and girardia binom files in the same folder too

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Import the raw dataset ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read_excel("Raw datasets.xlsx", sheet = "Behavioural assay data")
data$Dose <- as.factor(data$Dose) #convert dose to factor
data$Batch <- as.factor(data$Batch) #convert dose to factor
#count the number of dead planaria
dead_planaria <- data[which(is.na(data$Activity) == T), ] #find planaria with missing values (missing values are due to dead planaria)
n_dead <- length(unique(dead_planaria$ID))
dead_planaria
n_dead

#subset the data to remove planaria that died
#P1B3-335-Girardia
#P2B1-335-Girardia
#P4A1-335-Girardia
#P4B2-480-Girardia
#P7B3-480-Girardia
data <- subset(data, Dead != 1)

#count the number of planaria activity measurements for both species
table(data$Species)

#Plot activity per day and species
ggplot(data, aes(x = Activity))+
  geom_histogram()+
  facet_wrap(~Day+Species)+
  xlab("Activity (mm)")

#Plot activity per dose and species
ggplot(data, aes(x = Activity))+
  geom_histogram()+
  facet_wrap(~Dose+Species)+
  xlab("Activity (mm)")

#compute the mean of the explanatory variables 
mean_size <- mean(na.omit(data$Size))
mean_day <- mean(data$Day)
mean_temperature <- mean(data$Temperature)
mean_activity <- mean(data$Activity)

#center the variables (this helps model convergence)
data$SizeC <- data$Size - mean_size
data$DayC <- data$Day - mean_day
data$TemperatureC <- data$Temperature - mean_temperature

#count the number of missing values
sum(is.na(data$SizeC)) #92 missing values


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Figure cura binomial  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the model predicting movement probability
cura_binom <- readRDS("cura_binom.rds")

#Define the conditions to calculate the marginal effects for the binomial model
conditions <- data.frame(SizeC = 5- mean_size)

#compute the marginal effects
cura_binom_marginal <- conditional_effects(cura_binom, effects = 'DayC:Dose', conditions = conditions)

#create a plot
cura_binom_marginal <- cura_binom_marginal[[1]]
cura_binom_marginal$DayC <- cura_binom_marginal$DayC+mean_day
  
#create a color ramp for the plots
my_colors <- RColorBrewer::brewer.pal(10, "Spectral")[c(1,3,4,7,8,10)]
my_colors <- rev(my_colors)
my_colors[2] <- "skyblue"

p1 <- ggplot(cura_binom_marginal, aes(x = DayC, y = estimate__, group = Dose))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Dose), alpha = 0.15)+
  geom_line(aes(colour = Dose), linewidth = 2)+
  theme_classic()+
  xlab("Day")+
  ylab("Probability of moving (marginal effect)")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.1),)+
  ggtitle(expression(paste("(A) ", italic("Cura pinguis"), " - Movement probability")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors)+
  scale_fill_manual(values = my_colors)
p1



library(RColorBrewer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----Figure Girardia  binomial ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load the model
girardia_binom <- readRDS("girardia_binom.rds")

# Define the condition where AnimalSize is set to 0.5
conditions <- data.frame(SizeC = 8.9- mean_size,
                         Split = 0)

#compute the marginal effects
girardia_binom_marginal <- conditional_effects(girardia_binom , effects = 'DayC:Dose', conditions = conditions)

#create a plot
girardia_binom_marginal <- girardia_binom_marginal[[1]]
girardia_binom_marginal$DayC <- girardia_binom_marginal$DayC+mean_day

p2 <- ggplot(girardia_binom_marginal, aes(x = DayC, y = estimate__, group = Dose))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Dose), alpha = 0.15)+
  geom_line(aes(color = Dose), linewidth = 2)+
  theme_classic()+
  xlab("Day")+
  ylab("Probability of moving (marginal effect)")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.1),)+
  ggtitle(expression(paste("(B) ", italic("Girardia tigrina"), " - Movement probability")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors)+
  scale_fill_manual(values = my_colors)
p2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Figure Cura truncated Gaussian ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cura_trunc_gaussian <- readRDS("cura_gaussian.rds")

#define the conditions to fit the marginal effects
conditions <- data.frame(SizeC = 5.0- mean_size)

#compute the marginal effects
cura_trunc_gaussian_marginal <- conditional_effects(cura_trunc_gaussian, effects = 'DayC:Dose', conditions = conditions)

#create a plot
cura_trunc_gaussian_marginal <- cura_trunc_gaussian_marginal[[1]]
cura_trunc_gaussian_marginal$DayC <- cura_trunc_gaussian_marginal$DayC+mean_day

p3 <- ggplot(cura_trunc_gaussian_marginal, aes(x = DayC, y = estimate__, group = Dose))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Dose), alpha = 0.15)+
  geom_line(aes(color = Dose), linewidth = 2)+
  theme_classic()+
  xlab("Day")+
  ylab("Mean distance (mm, marginal effect)")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 600))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.1),)+
  ggtitle(expression(paste("(C) ", italic("Cura pinguis"), " - Distance travelled")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors)+
  scale_fill_manual(values = my_colors)
p3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----Figure Girardia truncated gaussian ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

girardia_trunc_gaussian <- readRDS("girardia_gaussian.rds")

#define the conditions to fit the marginal effects
conditions <- data.frame(SizeC = 8.9- mean_size,
                         Split = 0)

#compute the marginal effects
girardia_trunc_gaussian_marginal <- conditional_effects(girardia_trunc_gaussian, effects = 'DayC:Dose', conditions = conditions)

#create a plot
girardia_trunc_gaussian_marginal <- girardia_trunc_gaussian_marginal[[1]]
girardia_trunc_gaussian_marginal$DayC <- girardia_trunc_gaussian_marginal$DayC+mean_day

p4 <- ggplot(girardia_trunc_gaussian_marginal, aes(x = DayC, y = estimate__, group = Dose))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Dose), alpha = 0.5)+
  geom_line(aes(color = Dose), linewidth=2)+
  theme_classic()+
  xlab("Day")+
  ylab("Mean distance (mm, marginal effect)")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 600))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.1),)+
  ggtitle(expression(paste("(D) ", italic("Girardia tigrina"), " - Distance travelled")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors)+
  scale_fill_manual(values = my_colors)
p4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Figure Cura models combined  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Define the conditions to calculate the marginal effects for the binomial model
conditions <- data.frame(SizeC = 5- mean_size, TemperatureC = -0.337)

#compute the marginal effects for bimomial model
cura_binom_me <-  conditional_effects(cura_binom, effects = 'DayC:Dose', conditions = conditions, ndraws = 2000, spaghetti = T)
cura_binom_me <- attr(cura_binom_me[[1]], "spaghetti")

#compute the marginal effects for truncated gaussian model
cura_trunc_me  <- conditional_effects(cura_trunc_gaussian, effects = 'DayC:Dose', conditions = conditions, ndraws = 2000, spaghetti = T)
cura_trunc_me <- attr(cura_trunc_me[[1]], "spaghetti")

#check the two models marginal effects are setup the same
all.equal(cura_binom_me$sample_, cura_binom_me$sample__)

#calculate the activity by multiplying the marginal effect estimates of both models
cura_final <- cura_binom_me #duplicating one of the two dataset for convenience
cura_final$estimate__ <- cura_binom_me$estimate__ * cura_trunc_me$estimate__ 

#calculate the mean and 95%CI for each combinations
pts <- unique(cura_final$DayC)
ds <- unique(cura_final$Dose)
df_cura_final <- data.frame(Dose = NULL, DayC = NULL, estimate__ = NULL, lower__ = NULL, upper__  = NULL)

for(i in pts){
  
  #subset the data for the current day
  tmp_day <- subset(cura_final, DayC == i)
  
  for(j in ds){
    
    #subset the data for current dose
    tmp_df <- subset(tmp_day, Dose == j)
    tmp_estimate <- mean(tmp_df$estimate__)
    tmp_lower <- quantile(tmp_df$estimate_, 0.025)
    tmp_upper <- quantile(tmp_df$estimate_, 0.975)
    
    #create temporary vector to append data frame
    tmp_vec <- data.frame(Dose = tmp_df$Dose[1], 
                          DayC = tmp_df$DayC[1], 
                          estimate__ = tmp_estimate,
                          lower__  = tmp_lower,
                          upper__ = tmp_upper )
    
    #append the dataframe
    df_cura_final <- rbind(df_cura_final, tmp_vec)
    
  }
  
}

#plot the final data for cura
df_cura_final$DayC <- df_cura_final$DayC + mean_day
p5 <- ggplot(df_cura_final, aes(x = DayC, y = estimate__, group = Dose))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Dose), alpha = 0.5)+
  geom_line(aes(color = Dose), linewidth=2)+
  theme_classic()+
  xlab("Day")+
  ylab("Mean distance (mm, marginal effect)")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 600))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.95),)+
  ggtitle(expression(paste("(E) ", italic("Cura pinguis"), " - Estimated activity")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors)+
  scale_fill_manual(values = my_colors)
p5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Figure Girardia models combined  ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Define the conditions to calculate the marginal effects for the binomial model
conditions <- data.frame(SizeC = 8.9 - mean_size, Split = 0, TemperatureC = 0.3410083)

#compute the marginal effects for bimomial model
girardia_binom_me <-  conditional_effects(girardia_binom, effects = 'DayC:Dose', conditions = conditions, ndraws = 2000, spaghetti = T)
girardia_binom_me <- attr(girardia_binom_me[[1]], "spaghetti")

#compute the marginal effects for truncated gaussian model
girardia_trunc_me  <- conditional_effects(girardia_trunc_gaussian, effects = 'DayC:Dose', conditions = conditions, ndraws = 2000, spaghetti = T)
girardia_trunc_me <- attr(girardia_trunc_me[[1]], "spaghetti")

#check the two models marginal effects are setup the same
all.equal(girardia_binom_me$sample_, girardia_binom_me$sample__)

#calculate the activity by multiplying the marginal effect estimates of both models
girardia_final <- girardia_binom_me #duplicating one of the two dataset for convenience
girardia_final$estimate__ <- girardia_binom_me$estimate__ * girardia_trunc_me$estimate__ 

#calculate the mean and 95%CI for each combinations
pts <- unique(girardia_final$DayC)
ds <- unique(girardia_final$Dose)
df_girardia_final <- data.frame(Dose = NULL, DayC = NULL, estimate__ = NULL, lower__ = NULL, upper__  = NULL)

for(i in pts){
  
  #subset the data for the current day
  tmp_day <- subset(girardia_final, DayC == i)
  
  for(j in ds){
    
    #subset the data for current dose
    tmp_df <- subset(tmp_day, Dose == j)
    tmp_estimate <- mean(tmp_df$estimate__)
    tmp_lower <- quantile(tmp_df$estimate_, 0.025)
    tmp_upper <- quantile(tmp_df$estimate_, 0.975)
    
    #create temporary vector to append data frame
    tmp_vec <- data.frame(Dose = tmp_df$Dose[1], 
                          DayC = tmp_df$DayC[1], 
                          estimate__ = tmp_estimate,
                          lower__  = tmp_lower,
                          upper__ = tmp_upper )
    
    #append the dataframe
    df_girardia_final <- rbind(df_girardia_final, tmp_vec)
    
  }
  
}

#plot the final data for girardia
df_girardia_final$DayC <- df_girardia_final$DayC + mean_day
p6 <- ggplot(df_girardia_final, aes(x = DayC, y = estimate__, group = Dose))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Dose), alpha = 0.5)+
  geom_line(aes(color = Dose), linewidth=2)+
  theme_classic()+
  xlab("Day")+
  ylab("Mean distance (mm, marginal effect)")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 600))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.95),)+
  ggtitle(expression(paste("(F) ", italic("Girardia tigrina"), " - Estimated activity")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors)+
  scale_fill_manual(values = my_colors)
p6


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Probability of splitting for Girardia ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#susbet the data to the Girardia species
data_girardia <- subset(data, Species == "Girardia")

#reformat the data for the survival analysis
data_cox <- data.frame(ID = NULL, Time = NULL, Size = NULL, Dose = NULL, Censored = NULL)

for(i in unique(data_girardia$ID)){
  
  #subset the current planaria
  tmp <- subset(data_girardia, ID == i)
  
  #check if the planaria split
  split <- if(sum(tmp$Split) > 0) {
    split <- 1} else {
      split <- 0   
    }
  
  #code the time
  if(split == 1){
    
    time <- tmp$Day[which(tmp$Split > 0)[1]]
    censored <- 1 #event occurred
    
  } else {
    
    time <- 5
    censored <- 0 #censored individual
    
  }
  
  #append the data frame
  tmp_vec <- data.frame(ID = i, Time = time,  Size= tmp$Size[1], Dose = tmp$Dose[1], Censored = censored)  
  data_cox <- rbind(data_cox, tmp_vec)  
  
}

#create the survival plot
fit <- survfit(Surv(Time, Censored)~Dose, data = data_cox)

#survival plot 
p7 <- ggsurvfit(fit ,type = "risk")+
    add_censor_mark() +
    add_confidence_interval() +
    scale_ggsurvfit()+
    xlab("Time (Days)")+
    ylab("Porportion of fissioned planaria")+
    theme(legend.position = "right")+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.direction = "horizontal",
        legend.position = c(0.45, 0.91),)+
  ggtitle(expression(paste("(G) ", italic("Girardia tigrina"), " - Fission")))+
  theme(plot.title = element_text(size=16))+
  scale_color_manual(values = my_colors, labels = c("0", "335", "480", "670", "960"))+
  scale_fill_manual(values = my_colors, labels = c("0", "335", "480", "670", "960"))+
  guides(color=guide_legend("Dose"), fill = guide_legend("Dose"))
p7

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Combine all the plots ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pf <- ggarrange(p1, p2,p3, p4,p5, p6, p7, ncol = 2, nrow = 4)
pf

#change working directory to figure folder and export the figure
setwd("./Figures")
png("Figure 4.png", width = 25, height = 40, unit = "cm", res = 450)
pf
dev.off()
