
## Rcode from Weber et al, "Speciation dynamics and (non-)parallel evolution along an ecological gradient in African cichlid fish' ##

## Table of contents ##

# 1) Set path and load libraries (l.18-27)
# 2) Vector calculations 
#   2.a) Morphology  (l.30-147)  
#   2.b) Genetics (outliers) (l.151-278)
#   2.c) Genetics (non-outliers) (l.282-415)
# 3) Vector correlations 
#   3.a) ThetaP and ThetaG outliers / ThetaG (l.419-440)
#   3.b) DeltaLP and DeltaLG outliers / DeltaLG && LP and LG outliers / LG (l.443-464)
#   3.c) DM and ThetaP && FST and ThetaG outliers / ThetaG (l.467-513)
#   3.d) SGV and ThetaP / ThetaG outliers / ThetaG (l.516-562)

###########################################
## 1) Setting path and loading libraries ##
###########################################

setwd("path/to/input/files")
library(plyr)
library(NISTunits)
library(ecodist)
library(ape)
library(ade4)
library(xlsx)

########################################
## 2.a) vector analyses:  Morphology  ##
########################################
 
## Load morphological data
traitsP<-read.table("Morphological_data_vector_analyses.txt", header=TRUE, sep="\t")

# size correction: Body depth / Standard length
traitsP$BDSL <-traitsP$BD/traitsP$SL

## subset data per system ##
Kalambo1 <-subset (traitsP, pop == "KaL" | pop == "Ka1")
Kalambo2 <-subset (traitsP, pop == "KaL" | pop == "Ka2")
Chitili <-subset (traitsP, pop == "ChL" | pop == "Ch1")
Lunzua <-subset (traitsP, pop == "LzL" | pop == "Lz1")
Lufubu_Abur <-subset (traitsP, species=="Abur" & pop == "Lf2" | pop == "LfL")
Rusizi_Abur <-subset (traitsP, species=="Abur" & (pop == "RuL" | pop == "Ru1"))
Rusizi_Hsta <-subset (traitsP, species=="Hsta" & (pop == "RuL" | pop == "Ru1"))
Lufubu_Cten <-subset (traitsP, species=="Cten" & pop == "Lf2" | pop == "NdL")
Mbulu <-subset (traitsP, pop == "ch_" | pop == "mb_")

## calculate t-value for each trait per system ##

cols_to_testP <- c(8,10:45)

P_vectors_KaLamabo1 <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Kalambo1[[colname]] ~ Kalambo1$pop)$statistic
    return(data.frame(Kalambo1=t_val))
  })

P_vectors_KaLamabo2 <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Kalambo2[[colname]] ~ Kalambo2$pop)$statistic
    return(data.frame(Kalambo2=t_val))
  })

P_vectors_Chitili <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Chitili[[colname]] ~ Chitili$pop)$statistic
    return(data.frame(Chitili=t_val))
  })

P_vectors_Lunzua <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Lunzua[[colname]] ~ Lunzua$pop)$statistic
    return(data.frame(Lunzua=t_val))
  })

P_vectors_Lufubu_Abur <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Lufubu_Abur[[colname]] ~ Lufubu_Abur$pop)$statistic
    return(data.frame(Lufubu_Abur=t_val))
  })

P_vectors_Rusizi_Abur <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Rusizi_Abur[[colname]] ~ Rusizi_Abur$pop)$statistic
    return(data.frame(Rusizi_Abur=t_val))
  })

P_vectors_Rusizi_Hsta <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Rusizi_Hsta[[colname]] ~ Rusizi_Hsta$pop)$statistic
    return(data.frame(Rusizi_Hsta=t_val))
  })

P_vectors_Lufubu_Cten <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Lufubu_Cten[[colname]] ~ Lufubu_Cten$pop)$statistic
    return(data.frame(Lufubu_Cten=t_val))
  })

P_vectors_Mbulu <- ldply(
  cols_to_testP,
  function(colname) {
    t_val = t.test(Mbulu[[colname]] ~ Mbulu$pop)$statistic
    return(data.frame(Mbulu=t_val))
  })

## concatenate all dataframes in one to have phenotypic vectors ##
vectors_P <- cbind(P_vectors_KaLamabo1,P_vectors_KaLamabo2,P_vectors_Chitili,P_vectors_Lunzua,P_vectors_Lufubu_Abur,P_vectors_Rusizi_Abur,P_vectors_Rusizi_Hsta,P_vectors_Lufubu_Cten,P_vectors_Mbulu)

################################################################################################

## calculating angles between vectors (ThetaP) and convert in degrees ##
## Pearson correlations ##
P_pearson_pairwise <- cor(vectors_P[,1:9], vectors_P[,1:9], use = 'pairwise.complete.obs')
P_angles_in_radians <- acos(P_pearson_pairwise)
P_angles_in_degrees <- NISTradianTOdeg(P_angles_in_radians)

################################################################################################

## calculating multivariate Euclidean distance - vector lengths ##
norm_vec <- function(x) sqrt(sum(x^2))

vector_lengths_P <- cbind(norm_vec(vectors_P$Kalambo1),norm_vec(vectors_P$Kalambo2), norm_vec(vectors_P$Chitili), 
                        norm_vec(vectors_P$Lunzua), norm_vec(vectors_P$Lufubu_Abur), norm_vec(vectors_P$Rusizi_Abur),
                        norm_vec(vectors_P$Rusizi_Hsta), norm_vec(vectors_P$Lufubu_Cten), norm_vec(vectors_P$Mbulu))

colnames(vector_lengths_P)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")

## calculate delta length P ##
delta_Lp <- combn(vector_lengths_P,2,FUN=diff)

## converting vector to a distance matrix ##
delta_lp_pairwise <- full(delta_Lp) 
colnames(delta_lp_pairwise)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")
rownames(delta_lp_pairwise)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")

################################################################################################


###############################################
## 2.b) vector analyses: Genetics (outliers) ##
###############################################

## load data ##
traitsG<-read.table("Genetic_data_outliers_vector_analyses.txt")
names(traitsG) <- c(0:191)
colnames(traitsG)[1] <-"IND"
colnames(traitsG)[192] <-"SP_POP"
colnames(traitsG)[193] <-"SPECIES"
colnames(traitsG)[194] <-"POP"

#######################################

## subset data per system ##

G_Kalambo1 <-subset (traitsG, POP == "KaL" | POP == "Ka1")
G_Kalambo2 <-subset (traitsG, POP == "KaL" | POP == "Ka2")
G_Chitili <-subset (traitsG, POP == "ChL" | POP == "Ch1")
G_Lunzua <-subset (traitsG, POP == "LzL" | POP == "Lz1")
G_Lufubu_Abur <-subset (traitsG, SPECIES=="Abur" & POP == "Lf2" | POP == "LfL")
G_Rusizi_Abur <-subset (traitsG, SPECIES=="Abur" & (POP == "RuL" | POP == "Ru1"))
G_Rusizi_Hsta <-subset (traitsG, SPECIES=="Hsta" & (POP == "RuL" | POP == "Ru1"))
G_Lufubu_Cten <-subset (traitsG, SPECIES=="Chorei" & POP == "Lf2" | POP == "NdL")
G_Mbulu <-subset (traitsG, POP == "ChilaL" | POP == "MbuluR")

#######################################

## calculate t-value for each trait per system ##

# define which columns to test: first 78 principal components
cols_to_testG <- c(2:79)

G_vectors_KaLamabo1 <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Kalambo1[[colname]] ~ G_Kalambo1$POP)$statistic
    return(data.frame(Kalambo1=t_val))
  })

G_vectors_KaLamabo2 <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Kalambo2[[colname]] ~ G_Kalambo2$POP)$statistic
    return(data.frame(Kalambo2=t_val))
  })

G_vectors_Chitili <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Chitili[[colname]] ~ G_Chitili$POP)$statistic
    return(data.frame(Chitili=t_val))
  })

G_vectors_Lunzua <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Lunzua[[colname]] ~ G_Lunzua$POP)$statistic
    return(data.frame(Lunzua=t_val))
  })

G_vectors_Lufubu_Abur <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Lufubu_Abur[[colname]] ~ G_Lufubu_Abur$POP)$statistic
    return(data.frame(Lufubu_Abur=t_val))
  })

G_vectors_Rusizi_Abur <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Rusizi_Abur[[colname]] ~ G_Rusizi_Abur$POP)$statistic
    return(data.frame(Rusizi_Abur=t_val))
  })

G_vectors_Rusizi_Hsta <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Rusizi_Hsta[[colname]] ~ G_Rusizi_Hsta$POP)$statistic
    return(data.frame(Rusizi_Hsta=t_val))
  })

G_vectors_Lufubu_Cten <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Lufubu_Cten[[colname]] ~ G_Lufubu_Cten$POP)$statistic
    return(data.frame(Lufubu_Cten=t_val))
  })

G_vectors_Mbulu <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_Mbulu[[colname]] ~ G_Mbulu$POP)$statistic
    return(data.frame(Mbulu=t_val))
  })

## concatenate all dataframes in one to have genomic vectors (outliers) ##
vectors_G <- cbind(G_vectors_KaLamabo1,G_vectors_KaLamabo2,G_vectors_Chitili,G_vectors_Lunzua,G_vectors_Lufubu_Abur,G_vectors_Rusizi_Abur,G_vectors_Rusizi_Hsta,G_vectors_Lufubu_Cten,G_vectors_Mbulu)

#######################################

## calculating angles between vectors (ThetaG) and convert in degrees ##
## Pearson correlations ##
G_pearson_pairwise <- cor(vectors_G[,1:9], vectors_G[,1:9], use = 'pairwise.complete.obs')
G_angles_in_radians <- acos(G_pearson_pairwise)
G_angles_in_degrees <- NISTradianTOdeg(G_angles_in_radians)

#######################################

## calculating multivariate Euclidean distance - vector lengths ##
norm_vec <- function(x) sqrt(sum(x^2)) 

vector_lengths_G <- cbind(norm_vec(vectors_G$Kalambo1),norm_vec(vectors_G$Kalambo2), norm_vec(vectors_G$Chitili),
                          norm_vec(vectors_G$Lunzua), norm_vec(vectors_G$Lufubu_Abur), norm_vec(vectors_G$Rusizi_Abur),
                          norm_vec(vectors_G$Rusizi_Hsta), norm_vec(vectors_G$Lufubu_Cten), norm_vec(vectors_G$Mbulu))

colnames(vector_lengths_G)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")

#######################################

## calculate delta length G ##
delta_LG <- combn(vector_lengths_G,2,FUN=diff)

## converting vector to a distance matrix ##
delta_lG_pairwise <- full(delta_LG) 
colnames(delta_lG_pairwise)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")
rownames(delta_lG_pairwise)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")

################################################################################################


###################################################
## 2.c) vector analyses: Genetics (non-outliers) ##
###################################################

## load data ##
traitsG_non_out<-read.table("Genetic_data_non_outliers_vector_analyses.txt")
names(traitsG_non_out) <- c(0:191)
colnames(traitsG_non_out)[1] <-"IND"
colnames(traitsG_non_out)[192] <-"SP_POP"
colnames(traitsG_non_out)[193] <-"SPECIES"
colnames(traitsG_non_out)[194] <-"POP"

#######################################

## subset data per system ##
G_non_out_Kalambo1 <-subset (traitsG_non_out, POP == "KaL" | POP == "Ka1")
G_non_out_Kalambo2 <-subset (traitsG_non_out, POP == "KaL" | POP == "Ka2")
G_non_out_Chitili <-subset (traitsG_non_out, POP == "ChL" | POP == "Ch1")
G_non_out_Lunzua <-subset (traitsG_non_out, POP == "LzL" | POP == "Lz1")
G_non_out_Lufubu_Abur <-subset (traitsG_non_out, SPECIES=="Abur" & POP == "Lf2" | POP == "LfL")
G_non_out_Rusizi_Abur <-subset (traitsG_non_out, SPECIES=="Abur" & (POP == "RuL" | POP == "Ru1"))
G_non_out_Rusizi_Hsta <-subset (traitsG_non_out, SPECIES=="Hsta" & (POP == "RuL" | POP == "Ru1"))
G_non_out_Lufubu_Cten <-subset (traitsG_non_out, SPECIES=="Chorei" & POP == "Lf2" | POP == "NdL")
G_non_out_Mbulu <-subset (traitsG_non_out, POP == "ChilaL" | POP == "MbuluR")

#######################################

## calculate t-value for each trait per system ##

G_non_out_vectors_KaLamabo1 <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Kalambo1[[colname]] ~ G_non_out_Kalambo1$POP)$statistic
    return(data.frame(Kalambo1=t_val))
  })

G_non_out_vectors_KaLamabo2 <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Kalambo2[[colname]] ~ G_non_out_Kalambo2$POP)$statistic
    return(data.frame(Kalambo2=t_val))
  })

## PC5 invariant, replace with a very close value
G_non_out_Chitili[1,6]<-(-0.0010999999999999)

G_non_out_vectors_Chitili <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Chitili[[colname]] ~ G_non_out_Chitili$POP)$statistic
    return(data.frame(Chitili=t_val))
  })

## PC3 invariant, replace with a very close value
G_non_out_Lunzua[1,4]<-(-0.0020999999999999)

G_non_out_vectors_Lunzua <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Lunzua[[colname]] ~ G_non_out_Lunzua$POP)$statistic
    return(data.frame(Lunzua=t_val))
  })

## PC3 invariant in one pop, replace with a very close value
G_non_out_Lufubu_Abur[13,4]<-(-0.0008)
G_non_out_Lufubu_Abur[1,4]<-(-0.0017)

G_non_out_vectors_Lufubu_Abur <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Lufubu_Abur[[colname]] ~ G_non_out_Lufubu_Abur$POP)$statistic
    return(data.frame(Lufubu_Abur=t_val))
  })

G_non_out_vectors_Rusizi_Abur <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Rusizi_Abur[[colname]] ~ G_non_out_Rusizi_Abur$POP)$statistic
    return(data.frame(Rusizi_Abur=t_val))
  })

G_non_out_vectors_Rusizi_Hsta <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Rusizi_Hsta[[colname]] ~ G_non_out_Rusizi_Hsta$POP)$statistic
    return(data.frame(Rusizi_Hsta=t_val))
  })

G_non_out_vectors_Lufubu_Cten <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Lufubu_Cten[[colname]] ~ G_non_out_Lufubu_Cten$POP)$statistic
    return(data.frame(Lufubu_Cten=t_val))
  })

G_non_out_vectors_Mbulu <- ldply(
  cols_to_testG,
  function(colname) {
    t_val = t.test(G_non_out_Mbulu[[colname]] ~ G_non_out_Mbulu$POP)$statistic
    return(data.frame(Mbulu=t_val))
  })

## concatenate all dataframes in one ##
vectors_G_non_out <- cbind(G_non_out_vectors_KaLamabo1,G_non_out_vectors_KaLamabo2,G_non_out_vectors_Chitili,G_non_out_vectors_Lunzua,G_non_out_vectors_Lufubu_Abur,G_non_out_vectors_Rusizi_Abur,G_non_out_vectors_Rusizi_Hsta,G_non_out_vectors_Lufubu_Cten,G_non_out_vectors_Mbulu)

#######################################

## calculating angles between vectors (ThetaG non outliers) and convert in degrees ##
## Pearson correlations ##
G_non_out_pearson_pairwise <- cor(vectors_G_non_out[,1:9], vectors_G_non_out[,1:9], use = 'pairwise.complete.obs')
G_non_out_angles_in_radians <- acos(G_non_out_pearson_pairwise)
G_non_out_angles_in_degrees <- NISTradianTOdeg(G_non_out_angles_in_radians)

#######################################

## calculating multivariate Euclidean distance - vector lengths ##
norm_vec <- function(x) sqrt(sum(x^2)) 

vector_lengths_G_non_out <- cbind(norm_vec(vectors_G_non_out$Kalambo1),norm_vec(vectors_G_non_out$Kalambo2), norm_vec(vectors_G_non_out$Chitili),
                        norm_vec(vectors_G_non_out$Lunzua), norm_vec(vectors_G_non_out$Lufubu_Abur), norm_vec(vectors_G_non_out$Rusizi_Abur),
                        norm_vec(vectors_G_non_out$Rusizi_Hsta), norm_vec(vectors_G_non_out$Lufubu_Cten), norm_vec(vectors_G_non_out$Mbulu))

colnames(vector_lengths_G_non_out)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")

#######################################

## calculate delta length G non-outliers ##
delta_LG_non_out <- combn(vector_lengths_G_non_out,2,FUN=diff)

## converting vector to a distance matrix ##
delta_lG_non_out_pairwise <- full(delta_LG_non_out) 
colnames(delta_lG_non_out_pairwise)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")
rownames(delta_lG_non_out_pairwise)[1:9]<-c("Kalambo1","Kalambo2","Chitili","Lunzua","Lufubu_Abur","Rusizi_Abur","Rusizi_Hsta","Lufubu_Cten","Mbulu")

################################################################################################

####################################################################
###  3.a) correlations of Theta P and Theta G / Theta G outliers ###
####################################################################

# Mantel test to test the correlation between thetaP and thetaG (outliers) [symmatric matrices are required]
mantel.test(P_angles_in_degrees,G_angles_in_degrees, nperm=9999, alternative="two.sided", graph=F)

# Mantel test to test the correlation between thetaP and thetaG (non-outliers)
mantel.test(P_angles_in_degrees,G_non_out_angles_in_degrees, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
# load data #
all_data <- read.xlsx("summary_statistics_all_populations.xlsx",sheetIndex = 1,header=T)

# linear model thetaP and thetaG outliers #
LM_thetaP_ThetaG <- lm(all_data$Theta_P~all_data$Theta_G_outliers_78)
summary(LM_thetaP_ThetaG) 

# linear model thetaP and thetaG non-outliers #
LM_thetaP_ThetaG_non_out <- lm(all_data$Theta_P~all_data$Theta_G_78)
summary(LM_thetaP_ThetaG_non_out)


################################################################################################

##################################################################################################
###  3.b) correlations of Delta LP and Delta LG / Delta LG outliers && LP and LG / LG outliers ###
##################################################################################################


# Mantel test to test the correlation between Delta LP and Delta LG (outliers)
mantel.test(delta_lp_pairwise,delta_lG_pairwise, nperm=9999, alternative="two.sided", graph=F)

# Mantel test to test the correlation between Delta LP and Delta LG (non-outliers)
mantel.test(delta_lp_pairwise,delta_lG_non_out_pairwise, nperm=9999, alternative="two.sided", graph=F)

########################

# Linear model LP and LG (outliers)
lm_LP_LG <-lm(vector_lengths_P[1,]~vector_lengths_G[1,])
summary(lm_LP_LG)

# Linear model LP and LG (non-outliers)
lm_LP_LG_non_out <-lm(vector_lengths_P[1,]~vector_lengths_G_non_out[1,])
summary(lm_LP_LG_non_out)


################################################################################################

############################################################################################################################
###  3.c) correlations of Mahalanobis distance (DM) and Theta P // Genetic distance (FST) and Theta G / Theta G outliers ###
############################################################################################################################

# load data #
# Mahalanobis distances (DM) between lake (i.e. ancestral) populations
DM<-read.xlsx("summary_statistics_all_populations.xlsx", sheetIndex=2, header=T)
row.names(DM)<-DM$NA.
DM<-DM[,-1]
DM<-as.matrix(DM)

# Genetic distances (FST) between lake (i.e. ancestral) populations
FST<-read.xlsx("summary_statistics_all_populations.xlsx", sheetIndex=3, header=T)
row.names(FST)<-FST$NA.
FST<-FST[,-1]
FST<-as.matrix(FST)

##############

## Mantel test to test the correlation between Theta P and DM
mantel.test(P_angles_in_degrees,DM, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
## linear model Theta P and DM ##
LM_DM_thetaP <- lm(all_data$Theta_P~all_data$Mahalanobis_distance_all_sp)
summary(LM_DM_thetaP) 

##############

## Mantel test to test the correlation between Theta G (outliers) and FST
mantel.test(G_angles_in_degrees,FST, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
## linear model Theta G (outliers) and FST ##
LM_FST_thetaG <- lm(all_data$Theta_G_outliers_78~all_data$FST)
summary(LM_FST_thetaG)

##############

## Mantel test to test the correlation between Theta G (non-outliers) and FST
mantel.test(G_non_out_angles_in_degrees,FST, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
## linear model Theta G (non-outliers) and FST ##
LM_FST_thetaG_non_out <- lm(all_data$Theta_G_78~all_data$FST)
summary(LM_FST_thetaG_non_out)

################################################################################################

#######################################################################################################
###  3.d) correlations of standing genetic variation (SGV) and Theta P / Theta G / Theta G outliers ###
#######################################################################################################

# load data #
# Standing genetic variation (SGV) between lake (i.e. ancestral) populations
SGV<-read.xlsx("summary_statistics_all_populations.xlsx", sheetIndex=4, header=T)
row.names(SGV)<-SGV$NA.
SGV<-SGV[,-1]
SGV<-as.matrix(SGV)

# Remove Kalambo 1 / Kalambo 2 comparison as outlier (100% SGV), keep Kalambo 1 system
SGV_Ka1 <-SGV[-2,-2]
P_angles_in_degrees_Ka1 <- P_angles_in_degrees[-2,-2]
G_angles_in_degrees_Ka1 <- G_angles_in_degrees[-2,-2]
G_non_out_angles_in_degrees_Ka1 <- G_non_out_angles_in_degrees[-2,-2]

##############

# Mantel test to test the correlation between Theta P and SGV
mantel.test(P_angles_in_degrees_Ka1,SGV_Ka1, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
## linear model Theta P and SGV ##
LM_SGV_thetaP <- lm(all_data$Theta_P~all_data$X.SGV)
summary(LM_SGV_thetaP) 

##############

# Mantel test to test the correlation between Theta G (outliers) and SGV
mantel.test(G_angles_in_degrees_Ka1,SGV_Ka1, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
## linear model Theta G (outliers) and SGV ##
LM_SGV_thetaG <- lm(all_data$Theta_G_outliers_78~all_data$X.SGV)
summary(LM_SGV_thetaG)

##############

# Mantel test to test the correlation between Theta G (non-outliers) and SGV
mantel.test(G_non_out_angles_in_degrees_Ka1,SGV_Ka1, nperm=9999, alternative="two.sided", graph=F)

## correlation coefficients ##
## linear model Theta G (non-outliers) and SGV ##
LM_SGV_thetaG_non_out <- lm(all_data$Theta_G_78~all_data$X.SGV)
summary(LM_SGV_thetaG_non_out)

