#get the path
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

#Load file
Spore = read.csv("Initial spore concentration.csv", header = T)
Spore = gather(Spore,"Month","Count",-"ID")
colnames(Spore)[1] = "ID"
Spore$Count[Spore$Count == "" ] = NA
Spore$Count[Spore$Month == "Oct" & Spore$ID == "A"] = round(mean(as.numeric(Spore$Count[Spore$ID == "A"]), na.rm = T),0)
Spore$Count[Spore$Count == "<18"] = "4"

#farm_survey
farm_info_comb = read.csv("Farm Practices Survey.csv")
#View(farm_info_comb)
#str(farm_info)

#merge
Spore_comb = left_join(Spore, farm_info_comb, by = c("ID"="Farm"))


##clean merged data
Spore_comb = Spore_comb[-which(is.na(Spore_comb$Count)), ]
Spore_comb$Month = factor(Spore_comb$Month, 
                          levels = c("Jan","Feb","March","April","May","June","July","Aug","Sept","Oct","Nov","Dec"))

Spore_noD_comb = Spore_comb[Spore_comb$ID != "D",]
Spore_noD_comb$Count = as.numeric(Spore_noD_comb$Count)
Spore_noD_comb$Number.of.Milking.Cows = as.numeric(Spore_noD_comb$Number.of.Milking.Cows)
Spore_noD_comb$What.type.of.bedding.is.used.for.lactating.cows. = as.factor(Spore_noD_comb$What.type.of.bedding.is.used.for.lactating.cows.)
Spore_noD_comb$How.Often.is.Bedding.Topped.up..per.week.= as.numeric(Spore_noD_comb$How.Often.is.Bedding.Topped.up..per.week.)
Spore_noD_comb$How.often.are.alleyways.scraped.per.day. = as.numeric(Spore_noD_comb$How.often.are.alleyways.scraped.per.day.)
Spore_noD_comb$Are.paper.towels.or.laundered.towels.used.during.milking.preparation. = as.factor(Spore_noD_comb$Are.paper.towels.or.laundered.towels.used.during.milking.preparation.)
Spore_noD_comb$Is.detergent.used.on.Towels. = as.factor(Spore_noD_comb$Is.detergent.used.on.Towels.)
Spore_noD_comb$Is.bleach.used.on.towels. = as.factor(Spore_noD_comb$Is.bleach.used.on.towels.)
Spore_noD_comb$Are.towels.dried. = as.factor(Spore_noD_comb$Are.towels.dried.)
Spore_noD_comb$Is.the.holding.area.cleaned.during.milking..If.so..how.is.it.cleaned. = as.factor(Spore_noD_comb$Is.the.holding.area.cleaned.during.milking..If.so..how.is.it.cleaned.)
Spore_noD_comb$Are.udders.clipped.or.flamed..If.so..how.often..Per.year = as.numeric(Spore_noD_comb$Are.udders.clipped.or.flamed..If.so..how.often..Per.year)
Spore_noD_comb$Is.udder.cleanliness.scored.routinely..If.so..how.often..Per.year = as.numeric(Spore_noD_comb$Is.udder.cleanliness.scored.routinely..If.so..how.often..Per.year)
Spore_noD_comb$Is.teat.end.cleanliness.scored.routinely..If.so..how.often..Per.year = as.numeric(Spore_noD_comb$Is.teat.end.cleanliness.scored.routinely..If.so..how.often..Per.year)
Spore_noD_comb$Is.teat.end.condition.scored.routinely..If.so..how.often..Per.year = as.numeric(Spore_noD_comb$Is.teat.end.condition.scored.routinely..If.so..how.often..Per.year)

#View(Spore_noD_comb)
#str(Spore_noD_comb)

####power analysis(Means and SD)
#power = 0.8
#alpha = 0.05, the range of n is [22,4252]
#alpha = 0.1, the range of n is [18,3348]
#1) bedding type
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, What.type.of.bedding.is.used.for.lactating.cows.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(What.type.of.bedding.is.used.for.lactating.cows.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())

#manure solids and sand: if alpha = 0.05: n = 74; if alpha = 0.1: n = 58
#manure solids and sawdust: if alpha = 0.05: n = 100; if alpha = 0.1: n = 78
#sand and sawdust: if alpha = 0.05: n = 1234; if alpha = 0.1: n = 972


#2) bedding topped up frequency
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, How.Often.is.Bedding.Topped.up..per.week.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(How.Often.is.Bedding.Topped.up..per.week.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())
#freq = 2 v.s freq = 7 per week: if alpha = 0.05: n = 78;if alpha = 0.1: n = 62

#3) alleyways scraping frequency
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, How.often.are.alleyways.scraped.per.day.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(How.often.are.alleyways.scraped.per.day.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())

#4) paper towles
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Are.paper.towels.or.laundered.towels.used.during.milking.preparation.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Are.paper.towels.or.laundered.towels.used.during.milking.preparation.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt),num=n())


#5) bleach
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Is.bleach.used.on.towels.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Is.bleach.used.on.towels.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())
#No v.s Yes: if alpha = 0.05: n = 582;if alpha = 0.1: n = 458

#6) detergent
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Is.detergent.used.on.Towels.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Is.detergent.used.on.Towels.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())
#Yes v.s NA:if alpha = 0.05: n = 26;if alpha = 0.1: n = 20

#7) dried towels
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Are.towels.dried.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Are.towels.dried.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())
#Yes v.s NA: if alpha = 0.05: n = 26;if alpha = 0.1: n = 20

#8) Is holding area cleaned
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Is.the.holding.area.cleaned.during.milking..If.so..how.is.it.cleaned.) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Is.the.holding.area.cleaned.during.milking..If.so..how.is.it.cleaned.)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())
#Scraped v.s scraped/water flush: if alpha = 0.05: n = 656;if alpha = 0.1: n = 516
#scraped v.s water flush:if alpha = 0.05: n = 144;if alpha = 0.1: n = 112
#water flush v.s scraped/water flush:if alpha = 0.05: n = 40 ;if alpha = 0.1: n = 32


#9) udders clipped
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Are.udders.clipped.or.flamed..If.so..how.often..Per.year) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Are.udders.clipped.or.flamed..If.so..how.often..Per.year)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt), num=n())
#freq = 2 per year v.s freq = 4 per year: if alpha = 0.05: n = 4252;if alpha = 0.1: n = 3348

#10) udder cleanliness score
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Is.udder.cleanliness.scored.routinely..If.so..how.often..Per.year) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Is.udder.cleanliness.scored.routinely..If.so..how.often..Per.year)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt),num=n())
#freq = 0 per year v.s freq = 6 per year: if alpha = 0.05: n = 22;if alpha = 0.1: n = 18

#11) teat end cleanliness score
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Is.teat.end.cleanliness.scored.routinely..If.so..how.often..Per.year) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Is.teat.end.cleanliness.scored.routinely..If.so..how.often..Per.year)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt),num=n())
#freq = 0 per year v.s freq = 12 per year: if alpha = 0.05: n = 284;if alpha = 0.1: n = 224

#12) teat end condition scored 
Spore_noD_comb %>% 
  mutate(logCount=log10(Count))%>%
  group_by(ID, Is.teat.end.condition.scored.routinely..If.so..how.often..Per.year) %>%
  summarize(avgSportCt = mean(logCount, na.rm=TRUE))%>%
  group_by(Is.teat.end.condition.scored.routinely..If.so..how.often..Per.year)%>%
  summarize(farmSportCtAve=mean(avgSportCt), farmSportCtSD=sd(avgSportCt),num=n())


