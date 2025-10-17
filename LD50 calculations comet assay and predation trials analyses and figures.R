rm(list=ls())
library(brms)
library(parallel)
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(ggpubr)
library(posterior)
library(readxl) #package to open excel files into R 

#place raw dataset in the same folder as the R script file

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#----  LD50 for Girardia tigrina planarias ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the raw data
data_planaria <- read_excel("Raw datasets.xlsx", sheet = "LD50 data")

#subset the raw data for Girardia mesurements
dataG <- subset(data_planaria, Species == 'Girardia')

#detect the number of cores available on the computer
my.cores <- detectCores()

#fit the model
m1_brm <- brm(Dead~Dose, family = "bernoulli",
              data = dataG,
              warmup = 500,iter = 8000, thin=5,
              chains = 4, init = "random",
              seed = 12345,
              cores = my.cores,
              verbose = 1,
              control = list(adapt_delta = 0.99, max_treedepth = 15))

#spread the draws and compute 95% CI for the LD50
spreadG <- m1_brm %>%
  spread_draws(b_Intercept, b_Dose)

#compute the LD 50 and 95%CI
LD50G <- -spreadG[,4]/spreadG[,5] 
mean(LD50G[,1]) 
quantile(LD50G[,1], prob = 0.025)
quantile(LD50G[,1], prob = 0.975)

#calculate the mortality curve for Girardia tigrina
girardia_curve <- m1_brm %>%
  spread_draws(b_Intercept,  b_Dose) %>%
  mutate(MSESC = list(seq(0, 6400, 1))) %>% #the observed value range of MSESC
  unnest(MSESC) %>%
  mutate(pred = exp(b_Intercept + b_Dose*MSESC)/(1+exp(b_Intercept + b_Dose*MSESC))) %>%
  group_by(MSESC) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low = quantile(pred, prob = 0.025),
            pred_high = quantile(pred, prob = 0.975))
girardia_curve$Species <- "Girardia"  

#compute odd ratios for an increase of 100 mJ
slopes_100 <- exp(data.frame(spreadG[,5]*100)[,1])
quantile(slopes_100, c(0.025,0.5,0.975))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- LD50 for Cura pinguis planaria ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#subset the raw data for Cura mesurements
dataC <- subset(data_planaria, Species == 'Cura')

#detect the number of cores available on the computer
my.cores <- detectCores()

#fit the model
m2_brm <- brm(Dead~Dose, family = "bernoulli",
              data = dataC,
              warmup = 500,iter = 8000, thin=5,
              chains = 4, init = "random",
              seed = 12345,
              cores = my.cores,
              verbose = 1,
              control = list(adapt_delta = 0.99, max_treedepth = 15))

#spread the draws and compute 95% CI for the LD50
spreadC  <- m2_brm %>%
  spread_draws(b_Intercept, b_Dose)

#compute the LD 50 and 95%CI
LD50C <- -spreadC[,4]/spreadC[,5] 
mean(LD50C[,1]) 
quantile(LD50C[,1], prob = 0.025)
quantile(LD50C[,1], prob = 0.975)

cura_curve <-  m2_brm %>%
  spread_draws(b_Intercept,  b_Dose) %>%
  mutate(MSESC = list(seq(0, 6400, 1))) %>% #the observed value range of MSESC
  unnest(MSESC) %>%
  mutate(pred = exp(b_Intercept + b_Dose*MSESC)/(1+exp(b_Intercept + b_Dose*MSESC))) %>%
  group_by(MSESC) %>%
  summarise(pred_m = mean(pred, na.rm = TRUE),
            pred_low = quantile(pred, prob = 0.025),
            pred_high = quantile(pred, prob = 0.975))
  
cura_curve$Species <- "Cura"  

#compute odd ratios for an increase of 100 mJ
slopes_100 <- exp(data.frame(spreadC[,5]*100)[,1])
quantile(slopes_100, c(0.025,0.5,0.975))

#combine the results from the two species
data_model <- rbind(girardia_curve, cura_curve)

#plot the dose response curves for both species
p1 <- ggplot(data = data_model, aes(x = MSESC, y = pred_m, group = Species))+
  geom_ribbon(data = data_model,aes(ymin = pred_low, ymax = pred_high, fill = Species), alpha=0.2) +
  geom_line(data = data_model, aes(x = MSESC, y = pred_m, colour = Species),size = 1.5) +
  ylab("Probability of death after 96 hours") +
  xlab(bquote(bold("UVB dose ("~mJ/cm^2~")")))+
  labs(title = "(A) Median lethal dose")+
  geom_hline(yintercept = 0.5,linetype="dashed")+
  geom_vline(xintercept = c(1340,1930),linetype="dashed")+
  theme_classic()+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = c(0.75, 0.15),
        plot.margin = margin(r= 1, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))+
  scale_color_manual(values = c("#ff8030","#56B4E9"), name = "Species", labels = c(expression(italic("Cura pinguis")),expression(italic("Girardia tigrina"))))+
  scale_fill_manual(values = c( "#ff8030", "#56B4E9"), name = "Species", labels = c(expression(italic("Cura pinguis")),expression(italic("Girardia tigrina"))))
p1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Comet assay analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the raw dataset form the master file
data<- read_excel("Raw datasets.xlsx", sheet = "Comet assay data")
data$Dose <- factor(data$Dose, levels = c("0", "335", "480", "670", "960"))

#create a boxplot of the tDNA damage for each dose and species
p2<- ggplot(data, aes(y = tDNA, x = Dose, fill = Species))+
  geom_boxplot()+
  xlab(bquote(bold("UVB dose ("~mJ/cm^2~")")))+
  ylab("Proportion of DNA in the tail of the comet (tDNA)")+
  labs(title = "(B) Comet assay")+
  scale_fill_manual(values=c("#ff8030","#56B4E9"),labels = c(expression(italic("Cura pinguis")),expression(italic("Giradia trigrina"))))+
  theme_classic()+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        legend.key=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = c(0.2, 0.9),
        plot.margin = margin(r= 1, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))
p2

#fit a Gaussian model to investigate for differences in damages between species
model_damages <- brm(tDNA ~ Dose * Species,
                      family = "Gaussian",
                      data = data,
                      warmup = 500,iter = 8000, thin=5,
                      chains = 4, init = "random",
                      seed = 12345,
                      cores = my.cores,
                      verbose = 1,
                      control = list(adapt_delta = 0.99, max_treedepth = 15))
model_damages

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Predation trials analysis ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the dataset
data_predation <- read_excel("Raw datasets.xlsx", sheet = "Predation experiment data")
data_predation$Batch <- as.factor(data_predation$Batch)

#center the dose variable to help with model convergence
data_predation$DoseC <- data_predation$Dose - mean(data_predation$Dose)

#fit a model with batch as a random intercept
model_1 <- brm(Predated ~ 0 + Day + DoseC  + Batch,
               family = bernoulli(),
               data = data_predation,
               control = list(adapt_delta = 0.9),
               cores = parallel::detectCores())

#fit a model with batch as a random intercept
model_2 <- brm(Predated ~ 0 + Day * DoseC  + Batch,
               family = bernoulli(),
               data = data_predation,
               control = list(adapt_delta = 0.9),
               cores = parallel::detectCores())

#fit a model with  plates as a random intercept
model_3 <- brm(Predated ~ 0 + Day + DoseC + Batch  + (1|PlateID),
               family = bernoulli(),
               data = data_predation,
               control = list(adapt_delta = 0.9),
               cores = parallel::detectCores())

#fit a model with  plates as a random intercept
model_4 <- brm(Predated ~ 0 + Day * DoseC + Batch  + (1|PlateID),
               family = bernoulli(),
               data = data_predation,
               control = list(adapt_delta = 0.9),
               cores = parallel::detectCores())

#compare the performance of the four models
model_comp <- loo(model_1, model_2, model_3, model_4)

#check if the difference in elpd is large enough
abs(model_comp[[2]][,1]) / abs(model_comp[[2]][,2])

#plate ID should be retained in the model
#interaction doesn't improve the model

#plot the batch effect for model 4
conditional_effects(model_3, "Batch")
#no particular temporal trends in the batch effects

#select best model
best_model <- model_3

##### calculate the odds ratio for each categories
#extract the dras
mu_24H <- extract_variable(best_model, variable = "b_Day24H")
mu_48H <- extract_variable(best_model, variable = "b_Day48H")
mu_72H <- extract_variable(best_model, variable = "b_Day72H")
mu_Dose <- extract_variable(best_model, variable = "b_DoseC")

#get the probabilities fore each day
p_mu_24H <- exp(mu_24H)/ (1+exp(mu_24H))
p_mu_48H <- exp(mu_48H)/ (1+exp(mu_48H))
p_mu_72H <- exp(mu_72H)/ (1+exp(mu_72H))

#compute the odds ratios
ratio_48H <- p_mu_48H / p_mu_24H 
ratio_72H <- p_mu_72H / p_mu_24H 

#compute the odds ratio 95%CI
round(quantile(ratio_48H, c(0.025,0.5,0.975)), digits = 1)
round(quantile(ratio_72H, c(0.025,0.5,0.975)), digits = 1)
round(quantile(exp(mu_Dose*100), c(0.025,0.5,0.975)), digits = 2) #odds ratio for an increase of 100 mJ


dosage <- conditional_effects(best_model,"DoseC",conditions = data.frame(Day = c("24H", "48H", "72H")))[[1]]
dosage$DoseC <- dosage$DoseC + mean(data_predation$Dose)
  
p3 <- ggplot(data = dosage, aes(x = DoseC, y = estimate__, group = Day, color = Day))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = Day), alpha = 0.15, colour = NA)+
  geom_vline(xintercept = c(0, 335,480,690,960),linetype = "dashed")+
  geom_line(size = 2)+
  xlab(bquote(bold("UVB dose ("~mJ/cm^2~")")))+
  ylab("Average predation probability")+
  scale_color_manual(values = c("purple", "seagreen", "darkorange"))+
  scale_fill_manual(values = c("purple", "seagreen", "darkorange"))+
  theme_classic()+
  labs(colour = "Time (hours)", fill = "Time (hours)")+
  labs(title = "(C) UVB irradiated")+
  xlim(c(0,960))+
  ylim(c(0,1))+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        legend.key=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = c(0.84, 0.9),
        plot.margin = margin(r= 15, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))
p3             

#combine plots to export them for Figure 2
pf <- ggarrange(p1,p2,p3, align = c("hv"), nrow = 2, ncol = 2)
pf

#Export Figure 2 to be used in the manuscript
png("Figure 2.png", width = 30, height = 30, unit = "cm", res = 300)
pf
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#---- Supplementary metarial ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#fit a model with  plates as a random intercept
model_sup_mat <- brm(Predated ~ 0 + Day + DoseC + Irradiation + (1|PlateID),
               family = bernoulli(),
               data = data_predation,
               control = list(adapt_delta = 0.9, max_treedepth = 12),
               iter = 2000,
               cores = parallel::detectCores(),
               )

#compute the marginal effect for non irradiated planaria (UVB Dose = 0 mJ/cm2)
control <- conditional_effects(model_sup_mat,"Day:Irradiation", conditions = data.frame(DoseC = 0 - mean(data_predation$Dose)))[[1]]
colnames(control)[11] <- "lower"
colnames(control)[12] <- "upper"

p1 <- ggplot(data = control, aes(x = Day, y = estimate__, Group = Irradiation, color = Irradiation))+
  geom_point(data = control, size = 5,position = position_dodge(0.3))+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, position = position_dodge(0.3))+
  xlab(bquote("Time (hours)"))+
  ylim(c(0,1))+
  ylab("Average predation probability")+
  theme_classic()+
  scale_color_manual(values = c("purple", "seagreen"))+
  labs(title = "(A) No UVB irradiation")+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        axis.title=element_text(size=14,face="bold"),
        legend.position = c(0.25, 0.93),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.margin = margin(r= 1, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))
p1             


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot for 24 hour post irradiation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


dosage_24 <- conditional_effects(model_sup_mat,"DoseC:Irradiation",conditions = data.frame(Day = c("24H")))[[1]]
dosage_24$DoseC <- dosage_24$DoseC + mean(data_predation$Dose)

p2 <- ggplot(data = dosage_24, aes(x = DoseC, y = estimate__, group = Irradiation, color = Irradiation))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, group = Irradiation, fill = Irradiation), alpha = 0.15, colour = NA)+
  geom_vline(xintercept = c(0, 335,480,690,960),linetype = "dashed")+
  geom_line(data = dosage_24, size = 2)+
  xlab(bquote(bold("UVB dose ("~mJ/cm^2~")")))+
  ylab("Average predation probability")+
  scale_color_manual(values = c("purple", "seagreen"))+
  scale_fill_manual(values = c("purple", "seagreen"))+
  theme_classic()+
  labs(colour = "Irradiation", fill = "Irradiation")+
  labs(title = "(B) Predation at 24 hours")+
  xlim(c(0,960))+
  ylim(c(0,1))+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        legend.key=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = c(0.74, 0.9),
        plot.margin = margin(r= 15, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))
p2             

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot for 48 hours post irradiation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


dosage_48 <- conditional_effects(model_sup_mat,"DoseC:Irradiation",conditions = data.frame(Day = c("48H")))[[1]]
dosage_48$DoseC <-  dosage_48$DoseC + mean(data_predation$Dose)

p3 <- ggplot(data = dosage_48, aes(x = DoseC, y = estimate__, group = Irradiation, color = Irradiation))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, group = Irradiation, fill = Irradiation), alpha = 0.15, colour = NA)+
  geom_vline(xintercept = c(0, 335,480,690,960),linetype = "dashed")+
  geom_line(data = dosage_48, size = 2)+
  xlab(bquote(bold("UVB dose ("~mJ/cm^2~")")))+
  ylab("Average predation probability")+
  scale_color_manual(values = c("purple", "seagreen"))+
  scale_fill_manual(values = c("purple", "seagreen"))+
  theme_classic()+
  labs(colour = "Irradiation", fill = "Irradiation")+
  labs(title = "(C) Predation at 48 hours")+
  xlim(c(0,960))+
  ylim(c(0,1))+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        legend.key=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = c(0.76, 0.98),
        plot.margin = margin(r= 15, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))
p3        

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot for 72 hours post irradiation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dosage_72 <- conditional_effects(model_sup_mat,"DoseC:Irradiation",conditions = data.frame(Day = c("72H")))[[1]]
dosage_72$DoseC <-  dosage_72$DoseC  + mean(data_predation$Dose)

p4 <- ggplot(data = dosage_72, aes(x = DoseC, y = estimate__, group = Irradiation, color = Irradiation))+
  geom_ribbon(aes(ymin = lower__, ymax = upper__, group = Irradiation, fill = Irradiation), alpha = 0.15, colour = NA)+
  geom_vline(xintercept = c(0, 335,480,690,960),linetype = "dashed")+
  geom_line(data = dosage_72, size = 2)+
  xlab(bquote(bold("UVB dose ("~mJ/cm^2~")")))+
  ylab("Average predation probability")+
  scale_color_manual(values = c("purple", "seagreen"))+
  scale_fill_manual(values = c("purple", "seagreen"))+
  theme_classic()+
  labs(colour = "Irradiation", fill = "Irradiation")+
  labs(title = "(D) Predation at 72 hours")+
  xlim(c(0,960))+
  ylim(c(0,1))+
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold", size = 14),
        axis.text.y=element_text(face="bold", size = 14),
        legend.key=element_blank(),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.position = c(0.76, 0.98),
        plot.margin = margin(r= 15, l = 15, t = 5, unit = "pt"),
        strip.text.y = element_text(hjust=0.5,vjust = 0.5,angle=180,face="bold", size = 12))
p4     


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combine and export the plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pf <- ggarrange(p1,p2,p3,p4, align = c("hv"), nrow = 2, ncol = 2)
pf

png("Supplementary Figure 1.png", width = 30, height = 30, unit = "cm", res = 300)
pf
dev.off()