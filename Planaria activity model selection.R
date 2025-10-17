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
library(posterior)

#place raw dataset file in the same folder as this script

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Import the raw dataset ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- read_excel("Raw datasets.xlsx", sheet = "Behavioural assay data")
data$Dose <- as.factor(data$Dose) #convert dose to factor
data$Batch <- as.factor(data$Batch) #convert dose to factor
data$Moved <- ifelse(data$Activity == 0,0,1) #create a binary variable determining if planaria moved or now

#count the number of dead planaria
dead_planaria <- data[which(is.na(data$Activity) == T), ] #find planaria with missing values (missing values are due to dead planaria)
n_dead <- length(unique(dead_planaria$ID))
dead_planaria
n_dead

#subset the data to remove planarias that died
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----  Fit binomial models for Cura pinguis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#subset the activity data for the native species
data_cura <- subset(data, Species == "Cura")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with interaction between Dose and Day - No plate as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bf_binom_cura <- bf(Moved  ~    0 +  Dose *  DayC + SizeC + TemperatureC + Batch+  (DayC|ID))

model_binom_cura <- brm(bf_binom_cura , 
                   data = data_cura, 
                   family = bernoulli(),
                   cores = 4,
                   chains = 4,
                    seed = 1234)      

model_binom_cura

pp_check(model_binom_cura, ndraws = 100)

plot(conditional_effects(model_binom_cura, effects = 'DayC:Dose'), ask = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with interaction between Dose and Day - Plate included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bf_binom_cura_pre <- bf(Moved  ~    0+  Dose *  DayC + SizeC + TemperatureC + Batch+  (DayC|ID) + (1|PlateID))

model_binom_cura_pre <- brm(bf_binom_cura_pre , 
                        data = data_cura, 
                        family = bernoulli(),
                        cores = 4,
                        chains = 4,
                        seed = 1234)      

model_binom_cura_pre

pp_check(model_binom_cura_pre, ndraws = 100)

plot(conditional_effects(model_binom_cura_pre, effects = 'DayC:Dose'), ask = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with no interaction between Dose and Day - Plate ot included as random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bf_binom_cura_ni <- bf(Moved  ~   0 +  Dose +  DayC + SizeC + TemperatureC + Batch+  (DayC|ID))

model_binom_cura_ni <- brm(bf_binom_cura_ni, 
                        data = data_cura, 
                        family = bernoulli(),
                        cores = 4,
                        chains = 4,
                        seed = 1234)      

model_binom_cura_ni

pp_check(model_binom_cura_ni, ndraws = 100)

plot(conditional_effects(model_binom_cura_ni, effects = 'DayC:Dose'), ask = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with no interaction between Dose and Day - Plate included as random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bf_binom_cura_ni_pre <- bf(Moved  ~   0 +  Dose +  DayC + SizeC + TemperatureC + Batch + (DayC|ID) + (1|PlateID))

model_binom_cura_ni_pre <- brm(bf_binom_cura_ni_pre, 
                           data = data_cura, 
                           family = bernoulli(),
                           cores = 4,
                           chains = 4,
                           seed = 1234)       

model_binom_cura_ni_pre

pp_check(model_binom_cura_ni_pre, ndraws = 100)

plot(conditional_effects(model_binom_cura_ni_pre, effects = 'DayC:Dose'), ask = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Compare the fit of the four models models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model_comp <- loo(model_binom_cura,model_binom_cura_ni,
    model_binom_cura_pre,model_binom_cura_ni_pre)

#check if the different in elpd is greater than the standard error
abs(model_comp[[2]][,1])  / model_comp[[2]][,2] 

#all models have a relatively similar performance, therefore we select the most
#parsimonious one

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Plot the marginal effects for the best model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

best_cura_model_binom <- conditional_effects(model_binom_cura_ni,effects = c("DayC:Dose"))
plot(best_cura_model)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----  Fit the binomial models for Girardia tigrina ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_tigrina <- subset(data, Species == "Girardia") #subset the dataset for the invasive species
data_tigrina <- data_tigrina[-which(is.na(data_tigrina$Size) == T),] #remove points with no size measurements

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##----- Fit a model with interaction between Dose and Day - Plate not included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bf_binom_girardia <- bf(Moved  ~   0 + Dose *  DayC + SizeC + TemperatureC + Batch+ Split + (DayC|ID))

model_binom_girardia <- brm(bf_binom_girardia , 
                  data = data_tigrina, 
                  family = bernoulli(),
                  cores = 4,
                  chains = 4,
                  seed = 1234)      

model_binom_girardia

pp_check(model_binom_girardia, ndraws = 100)

plot(conditional_effects(model_binom_girardia, effects = 'DayC:Dose'), ask = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##----- Fit a model with interaction between Dose and Day - Plate included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bf_binom_girardia_pre <- bf(Moved  ~   0 + Dose *  DayC + SizeC + TemperatureC + Batch+ Split + (DayC|ID) + (1|PlateID))

model_binom_girardia_pre <- brm(bf_binom_girardia_pre , 
                            data = data_tigrina, 
                            family = bernoulli(),
                            cores = 4,
                            chains = 4,
                            seed = 1234)      

model_binom_girardia_pre

pp_check(model_binom_girardia_pre, ndraws = 100)

plot(conditional_effects(model_binom_girardia_pre, effects = 'DayC:Dose'), ask = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with no interaction between Dose and Day - Plate not included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
               
bf_binom_girardia_ni <- bf(Moved  ~   0 + Dose +  DayC + SizeC + TemperatureC + Batch+ Split + (DayC|ID))

model_binom_girardia_ni <- brm(bf_binom_girardia_ni , 
                            data = data_tigrina, 
                            family = bernoulli(),
                            cores = 4,
                            chains = 4,
                            seed = 1234)       

model_binom_girardia_ni

pp_check(model_binom_girardia_ni, ndraws = 100)

plot(conditional_effects(model_binom_girardia_ni, effects = 'DayC:Dose'), ask = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with no interaction between Dose and Day - Plate included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             

bf_binom_girardia_ni_pre <- bf(Moved  ~   0 + Dose +  DayC + SizeC + TemperatureC + Batch+ Split + (DayC|ID) +  (1|PlateID))

model_binom_girardia_ni_pre <- brm(bf_binom_girardia_ni_pre , 
                               data = data_tigrina, 
                               family = bernoulli(),
                               cores = 4,
                               chains = 4,
                               seed = 1234)       

model_binom_girardia_ni_pre

pp_check(model_binom_girardia_ni_pre, ndraws = 100)

plot(conditional_effects(model_binom_girardia_ni_pre, effects = 'DayC:Dose'), ask = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Compare the fit of the four models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model_comp <- loo(model_binom_girardia,model_binom_girardia_ni,
    model_binom_girardia_pre,model_binom_girardia_ni_pre)

#check if the model performance
abs(model_comp[[2]][,1]) /model_comp[[2]][,2]
  
#Model performances are very close from each other therefore we use the most parsimonious one
#this is the model without interaction between day and dose and no plate as a random effect

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Plot the marginal effects for the best model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

best_girardia_model_binom <- conditional_effects(model_binom_girardia_ni,effects = c("DayC:Dose"))
plot(best_girardia_model_binom)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----  Fit the activity model for Cura pinguis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#select Cura data
data_cura <- subset(data, Species == "Cura")

#select days for which planaria that moved
data_cura <- subset(data_cura, Activity > 0)

#create an histogram of Cura activities for each day and dose
p_cura <- ggplot(data_cura, aes(x = Activity))+
  geom_histogram()+
  facet_wrap(~Day+Dose)
p_cura #truncated Gaussian distribution is needed here

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with interaction between Dose and Day - Plate not included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bf_gaussian_cura <- bf(Activity | trunc(lb = 0) ~  0 + Dose *  DayC +  SizeC + TemperatureC + Batch +  (DayC|ID))

model_gaussian_cura <- brm(bf_gaussian_cura , 
                   data = data_cura , 
                   family = gaussian(),
                   cores = 4,
                   chains = 4,
                   seed = 1234)      

model_gaussian_cura

pp_check(model_gaussian_cura, ndraws = 100)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model with interaction between Dose and Day - Plate included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bf_gaussian_cura_pre <- bf(Activity | trunc(lb = 0) ~  0 + Dose *  DayC +  SizeC + TemperatureC + Batch +  (DayC|ID) + (1|PlateID))

model_gaussian_cura_pre <- brm(bf_gaussian_cura_pre , 
                           data = data_cura , 
                           family = gaussian(),
                           cores = 4,
                           chains = 4,
                           seed = 1234)      

model_gaussian_cura_pre

pp_check(model_gaussian_cura_pre, ndraws = 100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model without interaction between Dose and Day - Plate not included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bf_gaussian_cura_ni <- bf(Activity | trunc(lb = 0) ~  0 + Dose +  DayC + SizeC + TemperatureC + Batch + (DayC|ID))

model_gaussian_cura_ni <- brm(bf_gaussian_cura_ni, 
                           data = data_cura, 
                           family = gaussian(),
                           cores = 4,
                           chains = 4)      

model_gaussian_cura_ni

pp_check(model_gaussian_cura_ni, ndraws = 100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Fit a model without interaction between Dose and Day - Plate included as a random effect ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bf_gaussian_cura_ni_pre <- bf(Activity | trunc(lb = 0) ~  0 + Dose +  DayC + SizeC + TemperatureC + Batch + (DayC|ID) + (1|PlateID))

model_gaussian_cura_ni_pre <- brm(bf_gaussian_cura_ni_pre, 
                              data = data_cura, 
                              family = gaussian(),
                              cores = 4,
                              chains = 4)      

model_gaussian_cura_ni_pre

pp_check(model_gaussian_cura_ni_pre, ndraws = 100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Compare the fit of the four models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model_comp_gaussian <- loo(model_gaussian_cura,model_gaussian_cura_ni,
    model_gaussian_cura_pre,model_gaussian_cura_ni_pre)

#check if the difference in elpd is large enought
abs(model_comp_gaussian[[2]][,1]) / abs(model_comp_gaussian[[2]][,2])

#including plate ID as a random intercept doesn't improve the model much
#including an interaction between day and dose doesn't improve the model much

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Plot the marginal effects for the best model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

best_cura_model_gaussian <- conditional_effects(model_gaussian_cura_ni,effects = c("DayC:Dose"))
plot(best_cura_model_gaussian)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----  Fit the activity models for Girardia tigrina ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##---- Plot the marginal effects for the best model ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

best_girardia_model_gaussian <- conditional_effects(model_gaussian_girardia,effects = c("DayC:Dose"))
plot(best_girardia_model_gaussian)
