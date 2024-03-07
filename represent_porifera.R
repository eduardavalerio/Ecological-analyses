### --------------------------------------------------------------------------------------------------
#  Stack graph, representativeness and distribution of sponges species 
#
#   Author : Eduarda Valério de Jesus 
#   Last updated 2024/02/29
#
#   Script to create stack graphs to represent the distribution and caracterization of the sponges 
#   species from coast of São Paulo
#
#   https://rpubs.com/techanswers88/stackedbarcharts
### --------------------------------------------------------------------------------------------------

# install and load packages 
install.packages("ggplot2")
install.packages("ggthemes")
install.packages("scales")
install.packages("dplyr")
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)

# distribuição de frequência por batimetria (eixo x = batimetria -> intervalo; eixo y = frequência)

# upload planilha de dados (arquivo csv)
data <- read.csv("betadiv_por_bat - NCBI.csv", header = TRUE)

data$Site <- factor(data$Site, levels = c(">200", "180-200", "160-180", "140-160", "120-140", "100-120", "80-100", "60-80", "40-60", "20-40", "0-20"))
total <- sum(data$n)

ggplot(data, aes(x = n, y = Site)) +
  geom_bar(stat = "identity", fill = "dark gray") +  
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = "transparent")) +
  labs(x = "Número de espécies", y = "Batimetria") +
  ggtitle("graphic bar")
