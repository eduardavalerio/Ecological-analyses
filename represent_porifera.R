# REPRESENTATIVIDADE DE PORIFEROS EM BANCOS DE DADOS DE SEQUÊNCIAS GENÉTICAS 

# representação gráfica da representatividade de sequencias de poriferos em 
# bancos de dados para cada marcador molecular utilizado. 

# instalar e carregar pacotes 
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("viridis")
library(ggplot2)
library(tidyverse)
library(viridis)

# importar dados 
data <- read.csv("represent_poriferos.csv")

# Gráfico de pizza 
data %>%
  ggplot(aes(x = "", y = number.of.spp., fill = Order)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  scale_fill_viridis_d() +
  theme_void()

