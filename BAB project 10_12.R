##BAB calculated##
getwd()
setwd("~/../Desktop/Martin")
##Craigs Station##
#Load file
Spore = read.csv("Initial spore concentration.csv", header = T)
#View(Spore)

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

#######################9/27/2020#############################################
################## Build a model to predict sprore concentration ############


##1) merge the data with "Farm Practices Survey.csv"
Spore = read.csv("Initial spore concentration.csv", header = T)
colnames(Spore)[1] = "ID"
Spore = gather(Spore,"Month","Count",-"ID")
Spore$Count[Spore$Count == "" ] = NA
Spore$Count[Spore$Month == "Oct" & Spore$ID == "A"] = round(mean(as.numeric(Spore$Count[Spore$ID == "A"]), na.rm = T),0)
Spore$Count[Spore$Count == "<18"] = "4"

#clean survey data first
farm_info = read.csv("Farm Practices Survey_Edited 09222020.csv")
farm_info = farm_info[1:7,]
colnames(farm_info) = c("name","ID","bed_type","bed_top_freq","alle_scra_freq","detergent","bleach","hold_scra","hold_flush","udder_clip_freq","udder_clean_score","teat_end_clean","teat_end_cond")
str(farm_info)

farm_info$bed_type = as.character(farm_info$bed_type)
farm_info$bed_type[farm_info$bed_type == "Manure solids"] = "Manure Solids"

farm_info = subset(farm_info, 
                   select = -alle_scra_freq) # Remove alleway from variable list

farm_info$udder_clip_freq[c(2,3,4,7)] = "High"
farm_info$udder_clip_freq[-c(2,3,4,7)] = "Low" # Change udder clip freq to binary

farm_info$udder_clean_score[c(1,3,6,7)] = "Yes" # Change udder clean score to binary

farm_info$teat_end_clean[c(1,4,6,7)] = "Yes" # Change teat end clean to binary

farm_info$teat_end_cond[c(1,4,6)] = "Yes" # Change teat end cond to binary

farm_info$bed_type[6] = "Sand" #Lawnel's farm predominantly uses sand

farm_info$bed_top_freq[c(4,5,7)] ="High"
farm_info$bed_top_freq[-c(4,5,7)] = "Low" # Change bedding top frequency to binary

write.table(farm_info, file = "farm_info.xls", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)


#merge
Spore2 = left_join(Spore, farm_info, by = "ID")

#clean merged data
Spore3 = Spore2[-which(is.na(Spore2$Count)), ]
Spore3$Month = factor(Spore3$Month, 
                      levels = c("Jan","Feb","March","April","May","June","July","Aug","Sept","Oct","Nov","Dec"))

##removing farm D from analysis
Spore_noD = Spore3[Spore3$ID != "D",]
Spore_noD$Count = as.numeric(Spore_noD$Count)

Spore_noD$detergent = as.factor(Spore_noD$detergent)
Spore_noD$detergent = droplevels(Spore_noD$detergent)
Spore_noD$udder_clip_freq = as.factor(Spore_noD$udder_clip_freq)
Spore_noD$hold_scra = as.factor(Spore_noD$hold_scra)
Spore_noD$bed_type = as.factor(Spore_noD$bed_type)
Spore_noD$hold_flush= as.factor(Spore_noD$hold_flush)
Spore_noD$udder_clean_score = as.factor(Spore_noD$udder_clean_score)
Spore_noD$teat_end_cond = as.factor(Spore_noD$teat_end_cond)
Spore_noD$teat_end_clean = as.factor(Spore_noD$teat_end_clean)




##2)Exploratory analysis
#the vif value
fit = glmer(Count ~ bed_type + bed_top_freq + detergent  
       + bleach + hold_scra + hold_flush + udder_clip_freq + udder_clean_score 
       +teat_end_clean + teat_end_cond + (1|ID), data= Spore_noD, family = poisson)
vif(fit)


#2.1.1chart for 2 categorical variables:
# some mosaic plots = visualization of contigency table
#the p-value in the mosaic table is not reliable, need more data
library(vcd)
mosaic(~ hold_scra + alle_scra_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra+ hold_flush, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + bed_type, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + detergent, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + bleach, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + udder_clip_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_scra + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ alle_scra_freq + hold_flush, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + bed_type, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + detergent, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + bleach, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + udder_clip_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ alle_scra_freq + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ hold_flush + bed_type, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_flush + detergent, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_flush + udder_clip_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_flush + udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_flush + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ hold_flush + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ bed_type + detergent, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bed_type + bleach, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bed_type + udder_clip_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bed_type + udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bed_type + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bed_type + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ detergent + bleach, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ detergent + udder_clip_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ detergent + udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ detergent + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ detergent + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ bleach + udder_clip_freq, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bleach + udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bleach + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ bleach + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ udder_clip_freq +udder_clean_score, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ udder_clip_freq +teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ udder_clip_freq +teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ udder_clean_score + teat_end_clean, data = Spore_noD, legend = TRUE, shade = TRUE)
mosaic(~ udder_clean_score + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

mosaic(~ teat_end_clean + teat_end_cond, data = Spore_noD, legend = TRUE, shade = TRUE)

#mosaic(~ Sex + Survived, Titanic, legend = TRUE, shade = TRUE)

# contingency tables
with(Spore_noD, table(detergent, bed_type)) 

# fisher exact test
table1 = with(Spore_noD, table(alle_scra_freq, hold_scra))
fisher.test(table1)
#have association, they are dependent. 
#They are different since odds ratio = 0

# 2.1.2 chart for single vategorical variable:

# barplot
library(ggplot2)
ggplot(Spore_noD, aes(x = hold_scra)) + geom_bar() + ggtitle("Barplot for the variable 'hold_score'")

# pie chart
cat_table = table(Spore_noD$hold_scra)
pie_df = data.frame(group = names(cat_table), value = as.numeric(cat_table))
ggplot(pie_df, aes(x = "", y = value, fill = group)) + 
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggtitle("Pie chart for the variable 'hold_score'")

#2.1.3 Single categorical variables ~ Y count
# For two groups: boxplot for y vs x, and two sample t-test
boxplot(Count ~ hold_scra, data = Spore_noD)
boxplot(Count ~ alle_scra_freq, data = Spore_noD)
boxplot(Count ~ hold_flush, data = Spore_noD)


# parametric: t-test
t.test(Count ~ hold_scra, data = Spore_noD)
# non-parametric: wilcoxon rank sum test
wilcox.test(Count ~ hold_scra, data = Spore_noD)

# For more than two groups: boxplot and ANOVA
Spore_noD$udder_clean_score = factor(Spore_noD$udder_clean_score)
boxplot(Count ~ udder_clean_score, data = Spore_noD)
# parametric: ANOVA
summary(aov(Count ~ udder_clean_score, data = Spore_noD)) 
# non-parametric: kruskal-wallis test
kruskal.test(Count ~ udder_clean_score, data = Spore_noD) 

boxplot(Count ~ bed_type, data = Spore_noD)
boxplot(Count ~ detergent, data = Spore_noD)
boxplot(Count ~ bleach, data = Spore_noD)
boxplot(Count ~ udder_clip_freq, data = Spore_noD)
boxplot(Count ~ udder_clean_score, data = Spore_noD)
boxplot(Count ~ teat_end_clean, data = Spore_noD)
boxplot(Count ~ teat_end_cond, data = Spore_noD)




#2.2 For numerical variables:
# scatterplot
ggplot(Spore_noD, aes(x = bed_top_freq, y = Count)) + geom_point() + ggtitle("Scatterplot of Count against bed_top_freq")
ggplot(Spore_noD, aes(x = bed_top_freq, y = Count, color = hold_scra)) + geom_point() + ggtitle("Scatterplot of Count against bed_top_freq by hold_scra")
hist(Spore_noD$bed_top_freq)

##3) try lmer or glmer model
#If we put all the variables in the function and use summary(), there will be cross-impact between variables
#ex: one variable is important, but aftering putting another one, it becomes not important

# this is a function to check whether each variable is significant or not
check_significant = function(x){
  myformula = as.formula(paste("Count ~ ", x, "+ (1|ID)"))
  model = glmer(myformula, data = Spore_noD, family = poisson)
  null_model = glmer(Count~1 + (1|ID), data = Spore_noD, family = poisson)
  pval <- anova(model,null_model)[[8]][2]
  if ( pval < 0.1) {
    is.significant = "Yes"
  } else {
    is.significant = "No"
  }
  return(c(x, is.significant, pval ))
}


A = glmer(Count~detergent + (1|ID), data = Spore_noD, family = poisson)
B = glmer(Count~1 + (1|ID), data = Spore_noD, family = poisson)
anova(A,B)

# categorical variables
cat_variables = c("bed_type",
                  "detergent", 
                  "bleach",
                  "hold_scra",
                  "hold_flush",
                  "udder_clean_score",
                  "teat_end_clean",
                  "teat_end_cond",
                  "alle_scra_freq",
                  "udder_clip_freq") 

# numerical variables
num_variables =  "bed_top_freq"

# get the significant check table
significant_check_table = as.data.frame(t(sapply(c(cat_variables, num_variables), check_significant)))
names(significant_check_table) = c("varname", "is.significant", "p-value")
significant_check_table

#the result shows that udder_clean_score, teat_end_clean, and alle_scra_freq are significant variables

#######################################AIC/Backward Selection#######################
#1) Forward selection based on AIC
m.detergent = glmer(Count ~ detergent + (detergent|ID), data = Spore_noD, family = poisson)
m.alle_scra_freq = glmer(Count ~ alle_scra_freq + (alle_scra_freq|ID), data = Spore_noD, family = poisson)
m.udder_clip_freq = glmer(Count ~ udder_clip_freq + (udder_clip_freq|ID), data = Spore_noD, family = poisson)

summary(m.detergent);summary(m.alle_scra_freq);summary(m.udder_clip_freq)

#another way to include ID as a random variable gives different result:
m.detergent_1 = glmer(Count ~ detergent + (1|ID), data = Spore_noD, family = poisson)
m.alle_scra_freq_1 = glmer(Count ~ alle_scra_freq + (1|ID), data = Spore_noD, family = poisson)
m.udder_clip_freq_1 = glmer(Count ~ udder_clip_freq + (1|ID), data = Spore_noD, family = poisson)

summary(m.detergent_1);summary(m.alle_scra_freq_1);summary(m.udder_clip_freq_1)

# AIC of individual models
AICc(m.detergent);AICc(m.alle_scra_freq);AICc(m.udder_clip_freq) #only include udder_clip_freq
AICc(m.detergent_1);AICc(m.alle_scra_freq_1);AICc(m.udder_clip_freq_1) #only include alle_scra_freq


# Forward stepwise selection round 2
m2.base = glmer(Count ~ udder_clip_freq + (udder_clip_freq|ID), data = Spore_noD, family = poisson)
m2.detergent = glmer(Count ~ detergent + udder_clip_freq + (udder_clip_freq|ID) + (detergent|ID), data = Spore_noD, family = poisson)
m2.alle_scra_freq = glmer(Count ~ alle_scra_freq + udder_clip_freq+ (alle_scra_freq|ID) + (udder_clip_freq|ID), data = Spore_noD, family = poisson)

summary(m2.detergent);summary(m2.alle_scra_freq)

AICc(m2.base);AICc(m2.detergent);AICc(m2.alle_scra_freq)

# check improvement
lrtest(m2.detergent, m2.base)
lrtest(m2.alle_scra_freq, m2.base)


# Model with 3 sig variables
m =glmer(Count ~ udder_clean_score+
                 teat_end_clean +
                 (udder_clean_score+  teat_end_clean|ID), 
           data = Spore_noD, family = poisson)
summary(m)
#allow regression coefficient to be different in each farm

m =glmer(Count ~ udder_clean_score +
           teat_end_clean +
           (1|ID) , data = Spore_noD, family = poisson)
summary(m)

ggplot(Spore_noD, aes(Count)) +  geom_histogram(fill='white',color='black') + 
  facet_wrap(~udder_clean_score)

hist(log(Spore_noD$Count))

table(Spore_noD$ID)
str(Spore_noD)

#vif model inflation in the variation whether variables are confounded 

#since model selection has different selection criteria 
#we want to minimize sum of square error but do not want to increase the coefficient of determination ( R^2)  


##2) Backward/Stepwise selection
#we want to check whether forward and backward selection gives te same result
#Compared with forward selection , backward selection includes the situation when mutilple variables are put together

# Step 1: remove one predictor from the full model

m.full = glmer(Count ~ detergent +
                 towel_dried + 
                 alle_scra_freq+
                 udder_clip_freq +
                 (detergent|ID) +
                 (towel_dried|ID) +
                 (alle_scra_freq|ID) +
                 (udder_clip_freq|ID) , data = Spore_noD, family = poisson)

m.full.subtract.detergent = glmer(Count ~
                                    towel_dried + 
                                    alle_scra_freq+
                                    udder_clip_freq +
                                    (towel_dried|ID) +
                                    (alle_scra_freq|ID) +
                                    (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.towel_dried = glmer(Count ~ detergent +
                                      alle_scra_freq+
                                      udder_clip_freq +
                                      (detergent|ID) +
                                      (alle_scra_freq|ID) +
                                      (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.alle_scra_freq = glmer(Count ~ detergent +
                                         towel_dried + 
                                         udder_clip_freq +
                                         (detergent|ID) +
                                         (towel_dried|ID) +
                                         (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.udder_clip_freq = glmer(Count ~ detergent +
                                          towel_dried + 
                                          alle_scra_freq+
                                          (detergent|ID) +
                                          (towel_dried|ID) +
                                          (alle_scra_freq|ID) , data = Spore_noD, family = poisson)
AICc(m.full);AICc(m.full.subtract.detergent);AICc(m.full.subtract.towel_dried);AICc(m.full.subtract.alle_scra_freq);AICc(m.full.subtract.udder_clip_freq)

#The result shows that detergent and is_the_towel_dried are redundant variables - only need one

#To compare which to drop?deterent/towel_dried? - see step 2:

# Step 2: using m.full.subtract.detergent () 
m.full.subtract.detergent.towel_dried = glmer(Count ~
                                                alle_scra_freq+
                                                udder_clip_freq +
                                                (alle_scra_freq|ID) +
                                                (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.detergent.alle_scra_freq = glmer(Count ~
                                                   towel_dried + 
                                                   udder_clip_freq +
                                                   (towel_dried|ID) +
                                                   (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.detergent.udder_clip_freq = glmer(Count ~
                                                    towel_dried + 
                                                    alle_scra_freq+
                                                    (towel_dried|ID) +
                                                    (alle_scra_freq|ID) , data = Spore_noD, family = poisson)
AICc(m.full.subtract.detergent.towel_dried);AICc(m.full.subtract.detergent.alle_scra_freq);AICc(m.full.subtract.detergent.udder_clip_freq)

# Step 2: using m.full.subtract.towel_dried
m.full.subtract.towel_dried.detergent = glmer(Count ~ 
                                                alle_scra_freq+
                                                udder_clip_freq +
                                                (alle_scra_freq|ID) +
                                                (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.towel_dried.alle_scra_freq = glmer(Count ~ detergent +
                                                     udder_clip_freq +
                                                     (detergent|ID) +
                                                     (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.towel_dried.udder_clip_freq = glmer(Count ~ detergent +
                                                      alle_scra_freq+
                                                      (detergent|ID) +
                                                      (alle_scra_freq|ID) , data = Spore_noD, family = poisson)
AICc(m.full.subtract.towel_dried.detergent);AICc(m.full.subtract.towel_dried.alle_scra_freq);AICc(m.full.subtract.towel_dried.udder_clip_freq)

#The result shows that removing either one has the same effect
#so detergent and is_the_towel_dried are redundant variables!

# Step 3: After knowing the existence of redundant variables, remove one predictor(towel_dried ) from the full mode

  ## Run Backward Stepwise selection again to check the model
  
  # Step 3.1: remove towel_dried from the previous full model
 
#new full model
m.full = glmer(Count ~ detergent +
                   alle_scra_freq+
                   udder_clip_freq +
                   (detergent|ID) +
                   (alle_scra_freq|ID) +
                   (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.detergent = glmer(Count ~
                                    alle_scra_freq+
                                    udder_clip_freq +
                                    (alle_scra_freq|ID) +
                                    (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.alle_scra_freq = glmer(Count ~ detergent +
                                         udder_clip_freq +
                                         (detergent|ID) +
                                         (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.udder_clip_freq = glmer(Count ~ detergent +
                                          alle_scra_freq+
                                          (detergent|ID) +
                                          (alle_scra_freq|ID) , data = Spore_noD, family = poisson)
AICc(m.full);AICc(m.full.subtract.detergent);AICc(m.full.subtract.alle_scra_freq);AICc(m.full.subtract.udder_clip_freq)

#The AIC result shows that m.full.subtract.detergent has the smallest AIC value

# Step 3.2: using m.full.subtract.detergent
#the new full model then becomes m.full.subtract.detergent
#The two variables left are alle_scra_freq and udder_clip_freq
#compare 
m.full.subtract.detergent.alle_scra_freq = glmer(Count ~
                                                   udder_clip_freq +
                                                   (udder_clip_freq|ID) , data = Spore_noD, family = poisson)
m.full.subtract.detergent.udder_clip_freq = glmer(Count ~
                                                    alle_scra_freq+
                                                    (alle_scra_freq|ID) , data = Spore_noD, family = poisson)
AICc(m.full.subtract.detergent);AICc(m.full.subtract.detergent.alle_scra_freq);AICc(m.full.subtract.detergent.udder_clip_freq)
#The final result shows that m.full.subtract.detergent.alle_scra_freq has the smallest AIC c

########
#compare the two models 
m.full.subtract.detergent.alle_scra_freq #This is the final model selected from the backward selection after removing the redundant variable
m.full.subtract.detergent#This is base model selected after removing the redundant variable

#Compare these models and found m.full.subtract.detergent.alle_scra_freq is better 
#It shows we just need to include udder_clip_freq, so the forward and backward selection agrees!



