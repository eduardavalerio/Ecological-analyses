### --------------------------------------------------------------------------------------------------
#   Beta diversity based on bathymetry and locations
#
#   Author : Bruno Vilela, minor modif by Eduarda Valério de Jesus 
#   Last updated 2024/01/18
#
#   Script to perform beta diversity based on differents dephts and locations on the coast
#   of the state of São Paulo
### --------------------------------------------------------------------------------------------------


#beta diversity metrics - Bruno Vilela (https://youtu.be/nNXuMbMC8kc?si=8TNzgv4C8z7F87B7)

#beta diversity: the extent of change in community composition, or degree of community differentiation,
#in relation to an environmental gradient.
#variation in community composition from one location to another.

#beta diversity -> ratio between gamma diversity and alpha diversity
#ALFA = number of distinct species in a community, GAMA = number of distinct species overall (all communities analyzed)
#if alpha = 5, gamma = 5, beta will be 1.
#if alpha = 5, gamma = 15, beta = 3.

#SORENSEN SIMILARITY (BETAsor)
#Standardizes the results between 0 (totally different composition) and 1 (totally the same composition).
#FUNCTION: BETAsor = (BETA - 1)/(N - 1), where N = number of total communities
#Turnover -> species exchange
#Nestedness -> loss of species
#Turnover + nestedness -> exchange and loss of species
#same beta diversity can represent different situations

#SIMPSON SIMILARITY (BETAsim)
#measures only the variation due to the exchange of species (turnover)

#to calculate nestedness -> BETAnes = BETAsor - BETAsim

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

#------------------------------------bathymetry----------------------------------------
#data
batimetria<-read.csv("betadiv_por_bat.csv", header = TRUE, row.names = 1)
batimetria<-as.matrix(batimetria)

#total beta diversity

#use this code if the data is about abundance 
batimetria <- ifelse(batimetria > 0, 1, 0)

beta_total <- beta.multi(batimetria, index.family = "sorensen")
beta_total

#beta diversity pairwise  
beta_par <- beta.pair(batimetria, index.family = "sorensen")
beta_par
#see the parameters of beta diversity separately
beta_par$beta.sim
beta_par$beta.sne
beta_par$beta.sor

#Plot 
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
  scale_fill_gradientn(colours = pal, 
                       name="Beta Diversity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        text = element_text(size = 16)) +
  coord_fixed() +
  xlab("") +
  ylab("") +
  facet_wrap(metric ~.) +
  ggtitle("Components of beta diversity between depths")
g
dev.off()

#------------------------------------locations----------------------------------------
#data 
localidades<-read.csv("betadiv_por_reg.csv", header = TRUE, row.names = 1)
localidades<-as.matrix(localidades)

#total beta diversity

#use this code if the data is about abundance 
localidades <- ifelse(localidades > 0, 1, 0)

beta_total <- beta.multi(localidades, index.family = "sorensen")
beta_total

#beta diversity pairwise
beta_par <- beta.pair(localidades, index.family = "sorensen")
beta_par
#see the parameters of beta diversity separately
beta_par$beta.sim
beta_par$beta.sne
beta_par$beta.sor

#Plot 
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
  scale_fill_gradientn(colours = pal, 
                       name="Beta Diversity") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        text = element_text(size = 16)) +
  coord_fixed() +
  xlab("") +
  ylab("") +
  facet_wrap(metric ~.) +
  ggtitle("Components of beta diversity across locations")
g
dev.off()

