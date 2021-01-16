getwd()
setwd("~/../Desktop/BAB project")
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
library(emmeans)
library(afex)
library(extrafont)
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
#View(farm_info_comb)
#str(farm_info)

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

View(Spore_noD_comb)
#str(Spore_noD_comb)

#draw histogram of the spore count

font_import()
loadfonts(device="win")       #Register fonts for Windows bitmap output
fonts()                       #vector of font family names
##  [1] "Andale Mono"                  "AppleMyungjo"                
##  [3] "Arial Black"                  "Arial"                       
##  [5] "Arial Narrow"                 "Arial Rounded MT Bold"  
        
h= hist(log10(Spore_noD_comb$Count), breaks = 10,xlim = c(0.5,3), xlab = "log10(Spore Count) MPN/L",
        ylab = "Frequency", main = NULL,  cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
xfit = seq(min(log10(Spore_noD_comb$Count)), max(log10(Spore_noD_comb$Count)), length = 40) 
yfit = dnorm(xfit, mean = mean(log10(Spore_noD_comb$Count)), sd = sd(log10(Spore_noD_comb$Count))) 
yfit = yfit * diff(h$mids[1:2]) * length(log10(Spore_noD_comb$Count)) 
lines(xfit, yfit, col = "red", lwd = 1.1)

#Spore_noD_comb$Count = log10(Spore_noD_comb$Count)
ggplot(Spore_noD_comb, aes(x=ID, y=Count)) + geom_boxplot()+ theme_replace()+
  labs(title=NULL,y = "Spore Count (MPN/L)", x = "farm ID")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(text=element_text(family="Times", size=12))
 
## Multi-model inference
#The conditional R2 is the proportion of total variance explained through both fixed and random effects.
#Model 1: bedding only (4 levels) - combined bed_type and bed_top_freq
Spore_noD_comb$Bedding = relevel(Spore_noD_comb$Bedding, ref = "Manure Solids and high"  )
levels(Spore_noD_comb$Bedding)

fit_bed = lmer(log10(Count) ~ Bedding + (1|ID), data= Spore_noD_comb) #base level is Manure solids and high 
r.squaredGLMM(fit_bed)
summary(fit_bed) 
anova(fit_bed)
emmeans(fit_bed,pairwise~ Bedding, type = "response")
multcomp::cld(emmeans(fit_bed, ~ Bedding))

#Model 2:clean_agent only (4 levels) - combined detergent and bleach
Spore_noD_comb$Clean_agent = relevel(Spore_noD_comb$Clean_agent, ref = "det_Y_ble_Y" )
levels(Spore_noD_comb$Clean_agent)
fit_clean = lmer(log10(Count) ~ Clean_agent + (1|ID), data= Spore_noD_comb) #Base level is det_Y_ble_Y
r.squaredGLMM(fit_clean)
summary(fit_clean)
anova(fit_clean)
emmeans(fit_clean,pairwise~ Clean_agent, type = "response")
multcomp::cld(emmeans(fit_clean, ~ Clean_agent))
#Model 3: Score_variable only (5 levels) -combined teat_end_cond, teat_end_clean, udder_clean
levels(Spore_noD_comb$Score_variable)
fit_score = lmer(log10(Count) ~ Score_variable + (1|ID), data= Spore_noD_comb)#base level is udder_clean_N,teat_clean_N, teat_cond_N
r.squaredGLMM(fit_score)
summary(fit_score)
anova(fit_score)
emmeans(fit_score,pairwise~ Score_variable, type = "response")
multcomp::cld(emmeans(fit_score, ~ Score_variable))

#Model 4: hold_variable (3 levels) - combined hold_scra and hold_flush
fit_hold = lmer(log10(Count) ~ Hold_Variable + (1|ID), data= Spore_noD_comb) #Base level is scra_N and flush_Y
r.squaredGLMM(fit_hold)
summary(fit_hold)
anova(fit_hold)
emmeans(fit_hold,pairwise~ Hold_Variable, type = "response")
multcomp::cld(emmeans(fit_hold, ~ Hold_Variable))

#Model 5:udder_clip (2 levels) - no change 
fit_udder_clip = lmer(log10(Count) ~ Udder_clip + (1|ID), data= Spore_noD_comb) #Base level is udder_clip_high
r.squaredGLMM(fit_udder_clip)
summary(fit_udder_clip)
anova(fit_udder_clip)
emmeans(fit_udder_clip,pairwise~ Udder_clip, type = "response")
multcomp::cld(emmeans(fit_udder_clip, ~ Udder_clip))

###mean and standard deviation calculation
t.test(Spore_noD_comb$Count)
sd(Spore_noD_comb$Count)
max(Spore_noD_comb$Count)
aggregate(Spore_noD_comb$Count, list(Spore_noD_comb$ID), FUN=mean)
aggregate(log10(Spore_noD_comb$Count), list(Spore_noD_comb$ID), FUN=sd)
aggregate(Spore_noD_comb$Count, list(Spore_noD_comb$ID), FUN= IQR)
#comparison between farm F and C
farm_c = subset(Spore_noD_comb, ID == "C")
farm_f = subset(Spore_noD_comb, ID == "F")
summary(farm_f$Count)
