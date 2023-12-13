#métricas de beta diversidade - Bruno Vilela 

#beta diversidade: a extensão da mudança na composição da comunidade, ou grau de diferenciação da comunidade,
#em relação a um gradiente ambiental.
#variação da composição da comunidade de um local para outro.

#beta diversidade -> razão entre a diversidade gama e a diversidade alfa 
#ALFA = número de espécies distintas em uma comunidade, GAMA = número de espécies distintas no geral(todas comunidades analisadas)
#se alfa = 5, gama = 5, beta será 1.
#se alfa = 5, gama = 15, beta = 3. 

#DISSIMILARIDADE DE SORENSEN (BETAsor)
#Padroniza os resultados entre 0 (composição totalmente distintas) e 1 (composição totalmente iguais).
#FUNÇÃO: BETAsor = (BETA - 1)/(N - 1), sendo N = número de comunidades totais 
#Turnover -> troca de espécies
#Aninhamento -> perda de espécies
#Turnover + aninhamento -> troca e perda de espécies
#mesma beta diversidade pode representar diferentes situações

#DISSIMILARIDADE DE SIMPSON (BETAsim)
#mede apenas a variação devido a troca de espécies (turnover)

#para calcular o aninhamento -> BETAnes = BETAsor - BETAsim

#CALCULAR A BETA DIVERSIDADE 
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

######################################## batimetria ################################################
#para ler os próprios dados 
batimetria<-read.csv("betadiv_por_bat.csv", header = TRUE, row.names = 1)
batimetria<-as.matrix(batimetria)

#Beta diversidade total

#usar esse código caso os dados sejam de abundância 
batimetria <- ifelse(batimetria > 0, 1, 0)

beta_total <- beta.multi(batimetria, index.family = "sorensen")
beta_total

#Beta diversidade par a par 
beta_par <- beta.pair(batimetria, index.family = "sorensen")
beta_par
#para ver cada tipo de beta separadamente 
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

#################################################### localidades ######################################################
#para ler os próprios dados 
localidades<-read.csv("betadiv_por_reg.csv", header = TRUE, row.names = 1)
localidades<-as.matrix(localidades)

#Beta diversidade total

#usar esse código caso os dados sejam de abundância 
localidades <- ifelse(localidades > 0, 1, 0)

beta_total <- beta.multi(localidades, index.family = "sorensen")
beta_total

#Beta diversidade par a par 
beta_par <- beta.pair(localidades, index.family = "sorensen")
beta_par
#para ver cada tipo de beta separadamente 
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

