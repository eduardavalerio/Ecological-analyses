### --------------------------------------------------------------------------------------------------
#   BETA DIVERSITY
#
#   Author : Rodrigo Rodrigues Domingues, minor modif by Eduarda Valério de Jesus 
#   Last updated 2024/02/10 
#
#   Script to perform beta diversity between locations and bathymetry based on species
#   presence/absence data 
### --------------------------------------------------------------------------------------------------


#Beta diversity metrics  

#beta diversity: the extent of change in community composition, or degree of community differentiation,
#in relation to an environmental gradient.
#variation in community composition from one location to another.

#beta diversidade -> ratio between gama diversity and alfa diversity
#ALFA = number of species in a community, GAMA = number of species in general (all the communities)
#if alfa = 5, gama = 5, beta = 1.
#if alfa = 5, gama = 15, beta = 3. 

#SORENSEN SIMILARITY (BETAsor)
#Standardizes the results between 0 (completely equal composition) and 1 (completely diferent composition).
#function: BETAsor = (BETA - 1)/(N - 1), N = number of total comunities 
#turnover -> species exchange
#nestedness -> loss of species 
#turnover + nestedness -> exchange and loss of species 

#SIMPSON SIMILARITY (BETAsim)
#measures only the variation due to species turnover

#to calculate the nestedness -> BETAnes = BETAsor - BETAsim

#install and load packages  
install.packages("tidyverse")
install.packages("betapart")
install.packages("wesanderson")
install.packages("reshape2")
install.packages('readxl')
library(tidyverse)
library(betapart)
library(wesanderson)
library(reshape2)
library(readxl)
library(viridisLite)
library(viridis)

#save the output on a txt file
sink('output-betadiv-porifera-novo.txt')

######################################## bathymetry ################################################
#read the data frame
bat <- read.csv("porifera_bat_beta_novo.csv", header = TRUE, row.names = 1)
#turn ito a matrix
bat <- as.matrix(batimetria)

#total betadiversity

batimetria <- ifelse(batimetria > 0, 1, 0) #use this code if the data is in abundance and not presence/absence 

beta_total <- beta.multi(batimetria, index.family = "sorensen")
beta_total

#Beta diversity pair by pair
beta_par <- beta.pair(batimetria, index.family = "sorensen")
beta_par

#each type of beta diversity separately
beta_par$beta.sim
beta_par$beta.sne
beta_par$beta.sor

#plot 
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
sim <- get_upper_tri(as.matrix(beta_par$beta.sim))

sor <- get_upper_tri(as.matrix(beta_par$beta.sor))

sne <- get_upper_tri(as.matrix(beta_par$beta.sne))

melted_cormat <- reshape2::melt(sim)
melted_cormat2 <- reshape2::melt(sor)
melted_cormat3 <- reshape2::melt(sne)

melted_cormat <- rbind(melted_cormat2, melted_cormat, melted_cormat3)

melted_cormat$metric <- factor(rep(c("Total", "Turnover", "Nestedness"), 
                                   each = nrow(melted_cormat2)), 
                               c("Total", "Turnover", "Nestedness"))

pal <- wes_palette("Zissou1", 100, type = "continuous")

tiff(file = "betadiv_por_bat.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
g <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = viridis(11), 
                       name="Beta Diversity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust = 0),
        text = element_text(size = 16)) +
  coord_fixed() +
  xlab("") +
  ylab("") +
  facet_wrap(metric ~.) +
  ggtitle("Components of beta diversity between depths")
g
dev.off()

#################################################### locations ######################################################
#read the data frame 
localidades<-read.csv("porifera_reg_beta_novo.csv", header = TRUE, row.names = 1)
localidades<-as.matrix(localidades)

#total beta diversity
 
localidades <- ifelse(localidades > 0, 1, 0) #use this code if the data is in abundance and not presence/absence 


beta_total <- beta.multi(localidades, index.family = "sorensen")
beta_total

#Beta diversity pair by pair 
beta_par <- beta.pair(localidades, index.family = "sorensen")
beta_par

#each type of beta diversity separately
beta_par$beta.sim
beta_par$beta.sne
beta_par$beta.sor

#plot 
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
sim <- get_upper_tri(as.matrix(beta_par$beta.sim))

sor <- get_upper_tri(as.matrix(beta_par$beta.sor))

sne <- get_upper_tri(as.matrix(beta_par$beta.sne))

melted_cormat <- reshape2::melt(sim)
melted_cormat2 <- reshape2::melt(sor)
melted_cormat3 <- reshape2::melt(sne)

melted_cormat <- rbind(melted_cormat2, melted_cormat, melted_cormat3)

melted_cormat$metric <- factor(rep(c("Total", "Turnover", "Nestedness"), 
                                   each = nrow(melted_cormat2)), 
                               c("Total", "Turnover", "Nestedness"))

pal <- wes_palette("Zissou1", 100, type = "continuous")

tiff(file = "betadiv_por_reg.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
g <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = viridis(3), 
                       name="Beta Diversity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust = 0),
        text = element_text(size = 16)) +
  coord_fixed() +
  xlab("") +
  ylab("") +
  facet_wrap(metric ~.) +
  ggtitle("Components of beta diversity across locations")
g
dev.off()

sink()
