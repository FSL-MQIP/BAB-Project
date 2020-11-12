getwd()
setwd("~/../Desktop/Martin")
#Load pacakge
library(tidyr)
library(mice)
library(mitml)
library(fitdistrplus)
library(lme4)
library(dplyr)

library(lmtest)
library(MuMIn)
library(ggplot2)
library(car)
##Craigs Station##
#Load file
Spore = read.csv("Initial spore concentration.csv", header = T)
Spore = gather(Spore,"Month","Count",-"ID")
colnames(Spore)[1] = "ID"
Spore$Count[Spore$Count == "" ] = NA
Spore$Count[Spore$Month == "Oct" & Spore$ID == "A"] = round(mean(as.numeric(Spore$Count[Spore$ID == "A"]), na.rm = T),0)
Spore$Count[Spore$Count == "<18"] = "4"

#farm_survey
farm_info_comb = read.csv("Farm Practices Survey_Edited_11_11_20.csv")
str(farm_info)
#merge
Spore_comb = left_join(Spore, farm_info_comb, by = "ID")


##clean merged data
Spore_comb = Spore_comb[-which(is.na(Spore_comb$Count)), ]
Spore_comb$Month = factor(Spore_comb$Month, 
                      levels = c("Jan","Feb","March","April","May","June","July","Aug","Sept","Oct","Nov","Dec"))

#removing farm D from analysis
Spore_noD_comb = Spore_comb[Spore_comb$ID != "D",]
Spore_noD_comb$Count = as.numeric(Spore_noD_comb$Count)
Spore_noD_comb$Bedding = as.factor(Spore_noD_comb$Bedding)
Spore_noD_comb$Hold_Variable= as.factor(Spore_noD_comb$Hold_Variable)
Spore_noD_comb$Clean_agent = as.factor(Spore_noD_comb$Clean_agent)
Spore_noD_comb$Score_variable = as.factor(Spore_noD_comb$Score_variable)
Spore_noD_comb$Udder_clip = as.factor(Spore_noD_comb$Udder_clip)
colnames(Spore_noD_comb)[4] = "Farm Name"
#View(Spore_noD_comb)
str(Spore_noD_comb)


## Multi-model inference
#The conditional R2 is the proportion of total variance explained through both fixed and random effects.
#Model 1: bedding only (4 levels) - combined bed_type and bed_top_freq
fit_bed = lmer(Count ~ Bedding + (1|ID), data= Spore_noD_comb)
r.squaredGLMM(fit_bed)
summary(fit_bed)

#Model 2:clean_agent only (4 levels) - combined detergent and bleach
fit_clean = lmer(Count ~ Clean_agent + (1|ID), data= Spore_noD_comb)
r.squaredGLMM(fit_clean)
summary(fit_clean)

#Model 3: Score_variable only (5 levels) -combined teat_end_cond, teat_end_clean, udder_clean
fit_score = lmer(Count ~ Score_variable + (1|ID), data= Spore_noD_comb)
r.squaredGLMM(fit_score)
summary(fit_score)

#Model 4: hold_variable (3 levels) - combined hold_scra and hold_flush
fit_hold = lmer(Count ~ Hold_Variable + (1|ID), data= Spore_noD_comb)
r.squaredGLMM(fit_hold)
summary(fit_hold)

#Model 5:udder_clip (2 levels) - no change 
fit_udder_clip = lmer(Count ~ Udder_clip + (1|ID), data= Spore_noD_comb)
r.squaredGLMM(fit_udder_clip)
summary(fit_udder_clip)


