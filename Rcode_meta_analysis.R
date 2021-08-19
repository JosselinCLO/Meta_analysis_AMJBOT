##### Creator: XXXXX
##### Data: 29.06.2021
##### Article title: Short- and long-term consequences of genome doubling: a meta-analysis

# This script will allow you to reproduce the analyses from the paper
# For questions, please contact XXXXXX

# load the R package used for the analyses

library(ape)
library(matrixcalc)
library(metafor)
library(Matrix)
library(MASS)
library(pwr)
library(multcomp)
library(psych)
library(data.table)
library(lme4)

# load the dataset and the phylogenetic tree

data=read.table("Evolution_traits_polyploids2.csv",header = T, sep=";", dec=",")
data$Ratio = log(data$Ratio) # We transformed the data with a log transformation

tree<- read.tree("Phylo_correction") # The phylogenetic tree used for phylogenetic non-independence

# Division of the main dataset in subsets 

Vegetative = subset(data, data$Categorie == "Vegetative")
Phenological = subset(data, data$Categorie == "Phenological")
Floral = subset(data, data$Categorie == "Floral")
Macroscopic = rbind(Vegetative,Phenological, Floral)
Fitness = subset(data, data$Categorie == "Fitness")
Cellular_s = subset(data, data$Categorie == "Cellular - size")
Cellular_n = subset(data, data$Categorie == "Cellular - number")

####### 1. Phylogenetic analyses ##########

# I presented here how to test if including a phylogenetic correlation matrix
# improve the likelihood of the data, for the macroscopic traits

MetaData_f_All_Species_Data <- unique(Macroscopic$Species_Phylo) # Catch the species names in the Macroscopic dataset
Tree_f_ALL<-drop.tip(tree, tree$tip.label[-na.omit(match(MetaData_f_All_Species_Data, tree$tip.label))]) # Modifying the phylogenetic tree for the subset of species
plot(Tree_f_ALL) # Plot the tree
forcedC_Moderators <- as.matrix(forceSymmetric(vcv(Tree_f_ALL, corr=TRUE))) # For using rma.mv functions, the correlation matrix have to be Symmetric

vi = rnorm(nrow(Macroscopic), 0.001, 0.000000001)
# In meta-analyses, it is important to consider the sampling error of the raw data.
# As we are using a ratio, we have no sampling error, so we simulated it to be null on average

# Model with the phylogenetic correction

Model_REML_f_ALL1 = rma.mv(Ratio~Tetra_type,vi, data = Macroscopic, random = c( ~ 1 | Species, ~ 1 | Species_Phylo), R = list(Species_Phylo = forcedC_Moderators), method = "ML")

# Model without the phylogenetic correction

Model_REML_f_ALL2 = rma.mv(Ratio~Tetra_type,vi, data = Macroscopic, random = c( ~ 1 | Species), method = "ML")

## We tested the significance with a Likelihood ratio test between the full and restricted models

anova(Model_REML_f_ALL1,Model_REML_f_ALL2) 

## The LRT is not significant, so we removed the phylogenetic matrix from analyses
## and we only keep the random effect "Species"
## For testing the phylogenetic effect on the other subset, one just has to replace "Macroscopic"
## By the subset of its choice in the previous lines of code

####### 2. Analyses among cytotypes ##########

# I present here the analyses found in the main manuscript

### 2.1 Comparisons among diploid and established tetraploid values

### Macroscopic traits

Macroscopic_est = subset(Macroscopic, Macroscopic$Tetra_type=="Est")

summary(lm(Macroscopic_est$Value_diplo~Macroscopic_est$Value_tetra))

### Fitness traits

Fitness_est = subset(Fitness, Fitness$Tetra_type=="Est")

summary(lm(Fitness_est$Value_diplo~Fitness_est$Value_tetra))

### Cellular - size traits

Cellular_sest = subset(Cellular_s, Cellular_s$Tetra_type=="Est")

summary(lm(Cellular_sest$Value_diplo~Cellular_sest$Value_tetra))

### Cellular - number traits

Cellular_nest = subset(Cellular_n, Cellular_n$Tetra_type=="Est")

summary(lm(Cellular_nest$Value_diplo~Cellular_nest$Value_tetra))

## In all cases, the diploid value of a trait is an excellent predictor
## of the established tetraploid value

### 2.2 Are the lnRR from the different categories of traits significantly different ?

# For the macroscopic subcategories

modele1 = lmer(Macroscopic$Ratio~Macroscopic$Categorie+(1|Macroscopic$Species), REML=F)
modele2 = lmer(Macroscopic$Ratio~1+(1|Macroscopic$Species), REML=F)

anova(modele1, modele2)
# No, the subcategories do not vary significantly

# Fof the main categories:

# We first name the subcategories of Macroscopic traits "Macroscopic"
for(i in 1:(nrow(Macroscopic))){
  Macroscopic$Categorie[i] = "Macroscopic"
}

# We merged the new subsets and divided the new dataset for neo- and established tetraploids

data_comp = rbind(Macroscopic, Fitness, Cellular_n, Cellular_s)

data_comp_neo = subset(data_comp, data_comp$Tetra_type == "Neo")
data_comp_est = subset(data_comp, data_comp$Tetra_type == "Est")

# Do categories differ in neo lineages

modele1 = lmer(data_comp_neo$Ratio~data_comp_neo$Categorie+(1|data_comp_neo$Species), REML=F)
modele2 = lmer(data_comp_neo$Ratio~1+(1|data_comp_neo$Species), REML=F)

anova(modele1, modele2) # Yes they differ

# Do categories differ in established lineages

modele1 = lmer(data_comp_est$Ratio~data_comp_est$Categorie+(1|data_comp_est$Species), REML=F)
modele2 = lmer(data_comp_est$Ratio~1+(1|data_comp_est$Species), REML=F)

anova(modele1, modele2) # Yes they differ

### 2.3 Are the lnRR within categories of traits significantly different among new and established lineages

# Macroscopic traits

modele1 = lmer(Macroscopic$Ratio~Macroscopic$Tetra_type+(1|Macroscopic$Species), REML=F)
modele2 = lmer(Macroscopic$Ratio~1+(1|Macroscopic$Species), REML=F)

anova(modele1, modele2) # No significant differences

## We can extract the values for Figure 1.

Est_mor = summary(modele1)$coefficients[1]
Neo_mor = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_mor = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_mor = 1.96*(summary(modele1)$coefficients[4])

# Fitness traits

modele1 = lmer(Fitness$Ratio~Fitness$Tetra_type+(1|Fitness$Species), REML=FALSE)
modele2 = lmer(Fitness$Ratio~1+(1|Fitness$Species), REML=F)

anova(modele1, modele2) # Yes they differ significantly

## We can extract the values for Figure 1.

Est_fit = summary(modele1)$coefficients[1]
Neo_fit = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_fit = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_fit = 1.96*(summary(modele1)$coefficients[4])

# Cellular - size traits

modele1 = lmer(Cellular_s$Ratio~Cellular_s$Tetra_type+(1|Cellular_s$Species), REML=FALSE)
modele2 = lmer(Cellular_s$Ratio~1+(1|Cellular_s$Species), REML=F)

anova(modele1, modele2) # Yes they differ significantly

## We can extract the values for Figure 1.

Est_ces = summary(modele1)$coefficients[1]
Neo_ces = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_ces = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_ces = 1.96*(summary(modele1)$coefficients[4])

##### Cellular - number #####

modele1 = lmer(Cellular_n$Ratio~Cellular_n$Tetra_type+(1|Cellular_n$Species), REML=FALSE)
modele2 = lmer(Cellular_n$Ratio~1+(1|Cellular_n$Species), REML=F)

anova(modele1, modele2) # No they do not differ

## We can extract the values for Figure 1.

Est_cen = summary(modele1)$coefficients[1]
Neo_cen = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_cen = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_cen = 1.96*(summary(modele1)$coefficients[4])

# Make the figure 1

valeur1 = c(1,3,5,7)
valeur2 = c(2,4,6,8)

val1 = c(Neo_mor,Neo_ces,Neo_cen,Neo_fit)
IC1 = c(IC_Neo_mor,IC_Neo_ces,IC_Neo_cen,IC_Neo_fit)

val2 = c(Est_mor,Est_ces,Est_cen,Est_fit)
IC2 = c(IC_Est_mor,IC_Est_ces,IC_Est_cen,IC_Est_fit)

plot(val1~valeur1, pch = 16, xlim=c(1,8), col="cornflowerblue", ylim=c(-0.6,0.6), cex=1.5, ylab = "lnRR", xlab="", font.lab=2, cex.lab= 1.3)
points(val2~valeur2, pch = 16,col="tomato", cex=1.5)
abline(0,0, lty=2)
abline(v=2.5, lty=2)
abline(v=4.5, lty=2)
abline(v=6.5, lty=2)
arrows(valeur1,val1+IC1, valeur1,val1-IC1,lwd=1.5, angle=0,length=0.1,code=3, lty=1,col="cornflowerblue")
arrows(valeur2,val2+IC2, valeur2,val2-IC2,lwd=1.5, angle=0,length=0.1,code=3, lty=1,col="tomato")
legend(6.6,0.6,c("Neo - 4x","Est - 4x"), fill=c("cornflowerblue","tomato"), horiz=F, cex=1, box.lty=0, text.font = 2)

### 2.3 Are the lnRR within subcategories of traits significantly different among new and established lineages

# Floral traits

modele1 = lmer(Floral$Ratio~Floral$Tetra_type+(1|Floral$Species), REML=F)
modele2 = lmer(Floral$Ratio~1+(1|Floral$Species), REML=F)

anova(modele1, modele2) # No they do not differ significantly

## We can extract the values for Figure 2.

Est_flo = summary(modele1)$coefficients[1]
Neo_flo = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_flo = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_flo = 1.96*(summary(modele1)$coefficients[4])

# Vegetative traits 

modele1 = lmer(Vegetative$Ratio~Vegetative$Tetra_type+(1|Vegetative$Species), REML=F)
modele2 = lmer(Vegetative$Ratio~1+(1|Vegetative$Species), REML=F)

anova(modele1, modele2) # No they do not differ significantly

## We can extract the values for Figure 2.

Est_veg = summary(modele1)$coefficients[1]
Neo_veg = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_veg = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_veg = 1.96*(summary(modele1)$coefficients[4])

# Phenological traits

modele1 = lmer(Phenological$Ratio~Phenological$Tetra_type+(1|Phenological$Species), REML=F)
modele2 = lmer(Phenological$Ratio~1+(1|Phenological$Species), REML=F)

anova(modele1, modele2) # No they do not differ significantly

## We can extract the values for Figure 2.

Est_phe = summary(modele1)$coefficients[1]
Neo_phe = summary(modele1)$coefficients[1] + summary(modele1)$coefficients[2]

IC_Est_phe = 1.96*(summary(modele1)$coefficients[3])
IC_Neo_phe = 1.96*(summary(modele1)$coefficients[4])

anova(modele1, modele2)

## Make the figure 2

valeur = c(1:6)
effectif = c("(182)","(79)","(51)","(38)","(10)","(6)")
valeur1 = c(1,3,5)
valeur2 = c(2,4,6)

val1 = c(Neo_veg,Neo_flo,Neo_phe)
IC1 = c(IC_Neo_veg,IC_Neo_flo,IC_Neo_phe)

val2 = c(Est_veg,Est_flo,Est_phe)
IC2 = c(IC_Est_veg,IC_Est_flo,IC_Est_phe)

plot(val1~valeur1, pch = 16, xlim=c(1,6), col="cornflowerblue", ylim=c(-0.3,0.6), cex=1.5, ylab = "lnRR", xlab="", font.lab=2, cex.lab= 1.3)
points(val2~valeur2, pch = 16,col="tomato", cex=1.5)
abline(0,0, lty=2)
abline(v=2.5, lty=2)
abline(v=4.5, lty=2)
arrows(valeur1,val1+IC1, valeur1,val1-IC1,lwd=1.5, angle=0,length=0.1,code=3, lty=1,col="cornflowerblue")
arrows(valeur2,val2+IC2, valeur2,val2-IC2,lwd=1.5, angle=0,length=0.1,code=3, lty=1,col="tomato")
legend(5.2,0.6,c("Neo - 4x","Est - 4x"), fill=c("cornflowerblue","tomato"), horiz=F, cex=1, box.lty=0, text.font = 2)
text(valeur, -0.25, effectif)


#### Analyses on the Coefficient of variation ####

# We first transofmred back the data

Macroscopic$Ratio = exp(Macroscopic$Ratio)
Fitness$Ratio = exp(Fitness$Ratio)
Cellular_n$Ratio = exp(Cellular_n$Ratio)
Cellular_s$Ratio = exp(Cellular_s$Ratio)

# Macroscopic traits 

Macroscopic_1 = subset(Macroscopic, Macroscopic$Tetra_type=="Est")
Macroscopic_2 = subset(Macroscopic, Macroscopic$Tetra_type=="Neo")

# We list the articles used in each subcategories

list_neo = as.numeric(row.names(table(Macroscopic_2$Study_ID)))
list_est = as.numeric(row.names(table(Macroscopic_1$Study_ID)))

CV_morpho_neo = c(NULL)
CV_morpho_est = c(NULL)

# Neo lineages

k = 0

for(i in 1:length(list_neo)){
  
  # Subset for each article
  
  sub.data = subset(Macroscopic_2, Macroscopic_2$Study_ID == list_neo[i])
  
  if(nrow(sub.data) > 1){ # We compute the C.V if there are more than 2 values in the article
    
    k = k + 1
    
    CV_morpho_neo[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))
  }
}

# Established lineages

k = 0

for(i in 1:length(list_est)){
  
  sub.data = subset(Macroscopic_1, Macroscopic_1$Study_ID == list_est[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_morpho_est[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))}
}

CV_morphol_neo = mean(CV_morpho_neo)
CV_morphol_est = mean(CV_morpho_est)

sd_morphol_neo = sd(CV_morpho_neo)
sd_morphol_est = sd(CV_morpho_est)

wilcox.test(CV_morpho_est, CV_morpho_neo) # No significant differences

# Fitness traits

Fitness_1 = subset(Fitness, Fitness$Tetra_type=="Est")
Fitness_2 = subset(Fitness, Fitness$Tetra_type=="Neo")

list_neo = as.numeric(row.names(table(Fitness_2$Study_ID)))
list_est = as.numeric(row.names(table(Fitness_1$Study_ID)))

CV_fitness_neo = c(NULL)
CV_fitness_est = c(NULL)

# Neo lineages

k = 0

for(i in 1:length(list_neo)){
  
  sub.data = subset(Fitness_2, Fitness_2$Study_ID == list_neo[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_fitness_neo[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))
  }
}

# Established lineages

k = 0

for(i in 1:length(list_est)){
  
  sub.data = subset(Fitness_1, Fitness_1$Study_ID == list_est[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_fitness_est[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))}
}

CV_fit_neo = mean(CV_fitness_neo)
CV_fit_est = mean(CV_fitness_est)

sd_fit_neo = sd(CV_fitness_neo)
sd_fit_est = sd(CV_fitness_est)

wilcox.test(CV_fitness_est, CV_fitness_neo) # No significant differences

# Cellular-size traits

size_1 = subset(Cellular_s, Cellular_s$Tetra_type=="Est")
size_2 = subset(Cellular_s, Cellular_s$Tetra_type=="Neo")

list_neo = as.numeric(row.names(table(size_2$Study_ID)))
list_est = as.numeric(row.names(table(size_1$Study_ID)))

CV_size_neo = c(NULL)
CV_size_est = c(NULL)

# Neo lineages

k = 0

for(i in 1:length(list_neo)){
  
  sub.data = subset(size_2, size_2$Study_ID == list_neo[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_size_neo[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))
  }
}

# Established lineages

k = 0

for(i in 1:length(list_est)){
  
  sub.data = subset(size_1, size_1$Study_ID == list_est[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_size_est[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))}
}

CV_sizel_neo = mean(CV_size_neo)
CV_sizel_est = mean(CV_size_est)

sd_sizel_neo = sd(CV_size_neo)
sd_sizel_est = sd(CV_size_est)

wilcox.test(CV_size_est, CV_size_neo) # No significant differences

# Cellular-number traits 

number_1 = subset(Cellular_n, Cellular_n$Tetra_type=="Est")
number_2 = subset(Cellular_n, Cellular_n$Tetra_type=="Neo")

list_neo = as.numeric(row.names(table(number_2$Study_ID)))
list_est = as.numeric(row.names(table(number_1$Study_ID)))

CV_number_neo = c(NULL)
CV_number_est = c(NULL)

# Neo lineages

k = 0

for(i in 1:length(list_neo)){
  
  sub.data = subset(number_2, number_2$Study_ID == list_neo[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_number_neo[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))
  }
}

# Established lineages

k = 0

for(i in 1:length(list_est)){
  
  sub.data = subset(number_1, number_1$Study_ID == list_est[i])
  
  if(nrow(sub.data) > 1){
    
    k = k + 1
    
    CV_number_est[k] = abs(sd(sub.data$Ratio)/mean(sub.data$Ratio))}
}

CV_numberl_neo = mean(CV_number_neo)
CV_numberl_est = mean(CV_number_est)

sd_numberl_neo = sd(CV_number_neo)
sd_numberl_est = sd(CV_number_est)

wilcox.test(CV_number_est, CV_number_neo) # No significant differences

# We can make the Figure 3

valeur = c(1:8)
effectif = c("(36)","(17)","(18)","(4)","(6)","(2)","(9)","(6)")
valeur1 = c(1,3,5,7)
valeur2 = c(2,4,6,8)

val1 = c(CV_morphol_neo,CV_sizel_neo,CV_numberl_neo,CV_fit_neo)
IC1 = c(sd_morphol_neo,sd_sizel_neo,sd_numberl_neo,sd_fit_neo)

val2 = c(CV_morphol_est,CV_sizel_est,CV_numberl_est,CV_fit_est)
IC2 = c(sd_morphol_est,sd_sizel_est,sd_numberl_est,sd_fit_est)

plot(val1~valeur1, pch = 16, xlim=c(1,8), col="cornflowerblue", ylim=c(-0.1,0.6), cex=1.5, ylab = "CV(RR)", xlab="", font.lab=2, cex.lab= 1.3)
points(val2~valeur2, pch = 16,col="tomato", cex=1.5)
abline(0,0, lty=2)
abline(v=2.5, lty=2)
abline(v=4.5, lty=2)
abline(v=6.5, lty=2)
arrows(valeur1,val1+IC1, valeur1,val1-IC1,lwd=1.5, angle=0,length=0.1,code=3, lty=1,col="cornflowerblue")
arrows(valeur2,val2+IC2, valeur2,val2-IC2,lwd=1.5, angle=0,length=0.1,code=3, lty=1,col="tomato")
legend(6.6,0.6,c("Neo - 4x","Est - 4x"), fill=c("cornflowerblue","tomato"), horiz=F, cex=1, box.lty=0, text.font = 2)
text(valeur, -0.075, effectif)
