##################################################
##################################################
############ Dendrogram and NMDS #################
##################################################
##################################################


###### Install packages #####
install.packages('vegan')
install.packages('ggplot2')
install.packages('factoextra')
install.packages("dendextend")
install.packages('igraph')
library (vegan)
library (ggplot2)
library (factoextra)
library (dendextend)
library (igraph)

# load script_estatistica_metabarcoding
load("nmds&dendogram_porifera.R")

################################ Batimetria ##################################


### input file
porifera_bat <- read.csv('nmds_por_bat.csv', header = T, row.names = 1)
str (porifera_bat)
head(porifera_bat)


### transform data into a matrix w/ presence absence information

# bat.matrix<-as.matrix(bat)
com_bat = porifera_bat[,3:ncol(porifera_bat)]
m_com_bat = as.matrix(com_bat)


######################
######################
##### dendrogram #####
######################
######################

clust_bat<-vegdist (m_com_bat, method = "jaccard")
hc_bat <- hclust(clust_bat, method = "average")

#rename labels of dendrogram
labels_bat = c("0-20", "20-40", "40-60", "60-80", "80-100", "100-120", "120-140", "140-160", "160-180", "180-200", "200-250", ">250")
labels_bat

hc_bat$labels_bat <- labels_bat

#raw dendrogram graph
tiff(file = "dendogram_bat.tiff", width=10, heigh=10, 
     unit="in",res = 300, compression="lzw")
fviz_dend(hc_bat, cex = 1.5, main = "", lwd = 1.5, ylab = "")
dev.off()

### input file with Site
porifera_bat <- read.csv('nmds_por_bat.csv', header = T)
str (porifera_bat)
head(porifera_bat)


### transform data into a matrix w/ presence absence information

# bat.matrix<-as.matrix(bat)
com_bat = porifera_bat[,3:ncol(porifera_bat)]
m_com_bat = as.matrix(com_bat)

#customized dendrogram graph
#tiff(file = "nmds_porifera_bat.tiff", width=10, heigh=8, 
#     unit="in",res = 300, compression="lzw")
#fviz_dend(hc_bat, k = 2, # cortando em 2 grupos
#          cex = 1.5, # tamanho do rótulo
#          ylab = "",
#          main = "",
#          k_colors = c("#211B15", "#211B15"),
#          color_labels_by_k = TRUE, # cores por grupo
#          rect = TRUE, # Adicionar retângulo ao redor dos grupos
#          rect_border = c("#C8C5BE", "#C8C5BE"),
#          rect_fill = TRUE)
#dev.off()

#colored dendrogram
#fviz_dend(hc_bat, cex = 1.5, k = 2, # corte em 2 grupos
#         k_colors = "jco", main = "", lwd = 1.5, ylab = "")

#phylogenetic dendrogram
#fviz_dend(hc_bat, k = 2, k_colors = "jco",
#          type = "phylogenic", repel = TRUE, main = "", lwd = 1.5, ylab = "")


###############
###############
##### NMDS ####
###############
###############

# run NMDS analysis
set.seed(123)
nmds = metaMDS(m_com_bat, distance = "jaccard")
nmds

# See data distribution 
plot(nmds)
ordiplot (nmds, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds)$sites)

# add site information
data.scores$Site = porifera_bat$Site
data.scores$Local = porifera_bat$Local
head (data.scores)

#plot NMDS in ggplot2

tiff(file = "nmds_porifera_bat.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
nmds_bat_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 12, aes( colour = Site))+ 
  scale_color_viridis_d()+
 scale_x_continuous(limits=c(-1.90,1))+
 scale_y_continuous(limits=c(-3,1.5)) +
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 20), 
        legend.text = element_text(size = 20, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 20), 
        axis.title.x = element_text(face = "bold", size = 20, colour = "black"), 
        legend.title = element_text(size = 20, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Site", y = "NMDS2") 
nmds_bat_plot
dev.off()

##################
##################
##### ANOSIM #####
##################
##################

#ANOSIM #anosim necessita de replicas por localidades
ano = anosim(m_com_bat, porifera_bat$Local, distance = "jaccard", permutations = 9999)
ano

############################################################### localidades #####################################################################

### input file
porifera_reg <- read.csv('nmds_por_reg.csv', header = T, row.names = 1)
str (porifera_reg)
head(porifera_reg)


### transform data into a matrix w/ presence absence information

# fish.matrix<-as.matrix(fish)
com_reg = porifera_reg[,3:ncol(porifera_reg)]
m_com_reg = as.matrix(com_reg)


######################
######################
##### dendrogram #####
######################
######################

clust_reg<-vegdist (m_com_reg, method = "jaccard")
hc_reg <- hclust(clust_reg, method = "average")

#rename labels of dendrogram
labels_reg = c("AL", "CSS", "IB", "UT", "PL", "PLI", "PLC", "PLU", "SP", "PLB", "GRJ", "PLB","PLBI","OIB")
labels_reg

hc_reg$labels_reg <- labels_reg

#raw dendrogram graph
tiff(file = "dendogram_reg.tiff", width=10, heigh=10, 
     unit="in",res = 300, compression="lzw")
fviz_dend(hc_reg, cex = 1.5, main = "", lwd = 1.5, ylab = "")
dev.off()

##input file with Site
porifera_reg <- read.csv('nmds_por_reg.csv', header = T)
str (porifera_reg)
head(porifera_reg)


### transform data into a matrix w/ presence absence information

# fish.matrix<-as.matrix(fish)
com_reg = porifera_reg[,3:ncol(porifera_reg)]
m_com_reg = as.matrix(com_reg)

#customized dendrogram graph
#tiff(file = "nmds_porifera_reg.tiff", width=10, heigh=8, 
#     unit="in",res = 300, compression="lzw")
#fviz_dend(hc_eggs, k = 2, # cortando em 2 grupos
#          cex = 1.5, # tamanho do rótulo
#          ylab = "",
#          main = "",
#          k_colors = c("#211B15", "#211B15"),
#          color_labels_by_k = TRUE, # cores por grupo
#          rect = TRUE, # Adicionar retângulo ao redor dos grupos
#          rect_border = c("#C8C5BE", "#C8C5BE"),
#          rect_fill = TRUE)
#dev.off()

#colored dendrogram
#fviz_dend(hc_reg, cex = 1.5, k = 2, # corte em 2 grupos
#         k_colors = "jco", main = "", lwd = 1.5, ylab = "")

#phylogenetic dendrogram
#fviz_dend(hc_reg, k = 2, k_colors = "jco",
#          type = "phylogenic", repel = TRUE, main = "", lwd = 1.5, ylab = "")


###############
###############
##### NMDS ####
###############
###############

# run NMDS analysis
set.seed(123)
nmds = metaMDS(m_com_reg, distance = "jaccard")
nmds

# See data distribution 
plot(nmds)
ordiplot (nmds, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds)$sites)

# add site information
data.scores$Site = porifera_reg$Site
data.scores$Local = porifera_reg$Local
head (data.scores)

#plot NMDS in ggplot2

tiff(file = "nmds_porifera_reg.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
nmds_reg_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 12, aes( colour = Site))+ 
  scale_color_viridis_d()+
  scale_x_continuous(limits=c(-0.70,1))+ #se não aparecer nada, mudar a escala
  scale_y_continuous(limits=c(-1,0.5)) +
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 20), 
        legend.text = element_text(size = 20, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 20), 
        axis.title.x = element_text(face = "bold", size = 20, colour = "black"), 
        legend.title = element_text(size = 20, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Sites", y = "NMDS2") 
nmds_reg_plot
dev.off()

##################
##################
##### ANOSIM #####
##################
##################

#ANOSIM #anosim necessita de replicas por localidades
ano = anosim(m_com_reg, porifera_reg$Local, distance = "jaccard", permutations = 9999)
ano