getwd()
setwd("~/../Desktop/BAB project")

############Load pacakge############
library(tidyr)
library(mice)
library(mitml)
library(fitdistrplus)
library(lme4)
library(dplyr)
library(lmtest)
library(MuMIn)
library(ggplot2)
library(ggsignif)
library(car)
library(emmeans)
library(afex)
library(extrafont)

##########Import data###############

#Import initial spore concentration
Spore = read.csv("Initial spore concentration.csv", header = T)
Spore = gather(Spore,"Month","Count",-"ID")
colnames(Spore)[1] = "ID"
Spore$Count[Spore$Count == "" ] = NA
Spore$Count[Spore$Month == "Oct" & Spore$ID == "A"] = round(mean(as.numeric(Spore$Count[Spore$ID == "A"]), na.rm = T),0)
Spore$Count[Spore$Count == "<18"] = "4.5"

#Import farm_survey
farm_info_comb = read.csv("Farm Practices Survey_Edited_11_11_20.csv")

#merge two data sets
Spore_comb = left_join(Spore, farm_info_comb, by = "ID")

##clean merged data
Spore_comb = Spore_comb[-which(is.na(Spore_comb$Count)), ]
Spore_comb$Month = factor(Spore_comb$Month, 
                      levels = c("Jan","Feb","March","April","May","June","July","Aug","Sept","Oct","Nov","Dec"))

#removing farm D from analysis
Spore_noD_comb = Spore_comb[Spore_comb$ID != "D",]

#Reorder the farm ID
Farm_ID_new = c(A="A", B="B", C="C", 
                E="D", F="E", G="F", H="G")

Spore_noD_comb$ID = as.character(Farm_ID_new[Spore_noD_comb$ID])

#correct the property of each variable
Spore_noD_comb$Count = as.numeric(Spore_noD_comb$Count)
Spore_noD_comb$Bedding = as.factor(Spore_noD_comb$Bedding)
Spore_noD_comb$Hold_Variable= as.factor(Spore_noD_comb$Hold_Variable)
Spore_noD_comb$Clean_agent = as.factor(Spore_noD_comb$Clean_agent)
Spore_noD_comb$Score_variable = as.factor(Spore_noD_comb$Score_variable)
Spore_noD_comb$Udder_clip = as.factor(Spore_noD_comb$Udder_clip)


##########draw histogram of the spore count############
#change the front of the label 
windowsFonts(times = windowsFont("Times New Roman")) 

par(family = "times", font = 1.5, font.lab = 1.5, font.axis = 1.5)

#1) draw the overall distribution
tiff(filename = "figure1.tif",width = 5.512*1200, height = 6*1200,res = 1200)
h= hist(log10(Spore_noD_comb$Count), breaks = 10,xlim = c(0.5,3), ylim = c(0,20), xlab = expression("BAB spore count (log"[10]*" MPN/L)"),
        ylab = "Frequency", main = NULL,  cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
xfit = seq(min(log10(Spore_noD_comb$Count)), max(log10(Spore_noD_comb$Count)), length = 40) 
yfit = dnorm(xfit, mean = mean(log10(Spore_noD_comb$Count)), sd = sd(log10(Spore_noD_comb$Count))) 
yfit = yfit * diff(h$mids[1:2]) * length(log10(Spore_noD_comb$Count)) 
lines(xfit, yfit, col = "red", lwd = 1.5)
dev.off()

#2) draw the spore count by farm ID
Spore_noD_comb$log_count = log10(Spore_noD_comb$Count)
tiff(filename = "figure2.tif",width = 5.512*1200, height = 6*1200,res = 1200)

ggplot(Spore_noD_comb, aes(x=ID, y=log_count)) + geom_boxplot()+ theme_replace()+
  labs(title=NULL,y = expression("BAB spore count (log"[10]*" MPN/L)"), x = "farm ID")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+ geom_text(x=5, y=720, label="a",  size = 4)+ 
  geom_text(x=1, y=180, label="b",  size = 4)+
  geom_text(x=2, y=280, label="b",  size = 4)+
  geom_text(x=3, y=110, label="b",  size = 4)+
  geom_text(x=4, y=360, label="b",  size = 4)+
  geom_text(x=6, y=460, label="b",  size = 4)+
  geom_text(x=7, y=140, label="b",  size = 4)

dev.off()


############ Test whether each farm varies significantly different ########
anova_farm_mean =  lm(Count ~ ID, data = Spore_noD_comb)
summary(anova_farm_mean)
anova(anova_farm_mean)
emmeans(anova_farm_mean,pairwise~ID) #Tukey test

############# Multi-model inference #####################
#The conditional R2 is the proportion of total variance explained through both fixed and random effects.

#1)Model 1: bedding only (4 levels) - combined bed_type and bed_top_freq
Spore_noD_comb$Bedding = relevel(Spore_noD_comb$Bedding, ref = "Manure Solids and high"  )
levels(Spore_noD_comb$Bedding)

fit_bed = lmer(log10(Count) ~ Bedding + (1|ID), data= Spore_noD_comb) #base level is Manure solids and high 
r.squaredGLMM(fit_bed)
summary(fit_bed) 
anova(fit_bed)
emmeans(fit_bed,pairwise~ Bedding, type = "response")
multcomp::cld(emmeans(fit_bed, ~ Bedding))
plot(Spore_noD_comb$Bedding, Spore_noD_comb$Count, main="Distribution of spore count in bedding variable",
     ylab="Spore Count (MPN/L)", xlab = "level in bedding variable")

Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Bedding) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE), sdSportCt = sd(logCount, na.rm=TRUE))%>%
  group_by(Bedding)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())

#2) Model 2:clean_agent only (4 levels) - combined detergent and bleach
Spore_noD_comb$Clean_agent = relevel(Spore_noD_comb$Clean_agent, ref = "det_Y_ble_Y" )
levels(Spore_noD_comb$Clean_agent)
fit_clean = lmer(log10(Count) ~ Clean_agent + (1|ID), data= Spore_noD_comb) #Base level is det_Y_ble_Y
r.squaredGLMM(fit_clean)
summary(fit_clean)
anova(fit_clean)
emmeans(fit_clean,pairwise~ Clean_agent, type = "response")
multcomp::cld(emmeans(fit_clean, ~ Clean_agent))
plot(Spore_noD_comb$Clean_agent,Spore_noD_comb$Count, main="Distribution of spore count in clean_agent",
     ylab="Spore Count (MPN/L)", xlab = "level in clean_agent")

Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Clean_agent) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE), sdSportCt = sd(logCount, na.rm=TRUE))%>%
  group_by(Clean_agent)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())

#3) Model 3: Score_variable only (5 levels) -combined teat_end_cond, teat_end_clean, udder_clean
levels(Spore_noD_comb$Score_variable)
fit_score = lmer(log10(Count) ~ Score_variable + (1|ID), data= Spore_noD_comb)#base level is udder_clean_N,teat_clean_N, teat_cond_N
r.squaredGLMM(fit_score)
summary(fit_score)
anova(fit_score)
emmeans(fit_score,pairwise~ Score_variable, type = "response")
multcomp::cld(emmeans(fit_score, ~ Score_variable))
plot(Spore_noD_comb$Score_variable, Spore_noD_comb$Count, main="Distribution of spore count in score variable",
     ylab="Spore Count (MPN/L)", xlab = "level in score variable")

#4) Model 4: hold_variable (3 levels) - combined hold_scra and hold_flush
fit_hold = lmer(log10(Count) ~ Hold_Variable + (1|ID), data= Spore_noD_comb) #Base level is scra_N and flush_Y
r.squaredGLMM(fit_hold)
summary(fit_hold)
anova(fit_hold)
emmeans(fit_hold,pairwise~ Hold_Variable, type = "response")
multcomp::cld(emmeans(fit_hold, ~ Hold_Variable))
plot(Spore_noD_comb$Hold_Variable, Spore_noD_comb$Count, main="Distribution of spore count in hold_variable",
     ylab="Spore Count (MPN/L)", xlab = "level in hold_variable")


#5) Model 5:udder_clip (2 levels) - no change 

fit_udder_clip = lmer(log10(Count) ~ Udder_clip + (1|ID), data= Spore_noD_comb) #Base level is udder_clip_high
r.squaredGLMM(fit_udder_clip)
summary(fit_udder_clip)
anova(fit_udder_clip)
emmeans(fit_udder_clip,pairwise~ Udder_clip, type = "response")
multcomp::cld(emmeans(fit_udder_clip, ~ Udder_clip))
plot(Spore_noD_comb$Udder_clip, Spore_noD_comb$Count, main="Distribution of spore count in udder_clip variable",
     ylab="Spore Count (MPN/L)", xlab = "level in udder_clip")

######Calculation of the mean, standard deviation and SEM of spore counts from all farms####
t.test(Spore_noD_comb$Count)
#mean
overall_mean = mean(Spore_noD_comb$Count)
overall_log_mean = log10(mean(Spore_noD_comb$Count))
overall_mean
overall_log_mean
#standard deviation
overall_sd = sd(Spore_noD_comb$Count)
overall_log_mean = log10(sd(Spore_noD_comb$Count))
overall_sd
overall_log_mean
#standard error of the mean = SEM
overall_sem = sd(Spore_noD_comb$Count)/sqrt(84)
overall_sem
log_overall_sem = log10(sd(Spore_noD_comb$Count))/sqrt(84)
log_overall_sem
#mean_by_farm
mean_per_farm = aggregate(Spore_noD_comb$Count, list(Spore_noD_comb$ID), FUN=mean)
mean_per_farm
#sd_by_farm
sd_per_farm = aggregate(Spore_noD_comb$Count, list(Spore_noD_comb$ID), FUN= sd)
log_sd_per_farm = aggregate(log10(Spore_noD_comb$Count), list(Spore_noD_comb$ID), FUN=sd)
#IQR_by_farm
log_IQR_per_farm = aggregate(Spore_noD_comb$Count, list(Spore_noD_comb$ID), FUN= IQR)
#comparison between farm F and C
farm_c = subset(Spore_noD_comb, ID == "C")
farm_f = subset(Spore_noD_comb, ID == "F")
summary(farm_f$Count)

