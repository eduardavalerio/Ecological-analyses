### --------------------------------------------------------------------------------------------------
#  Genetic distance and Phylogenetic tree based on the molecular markers COI, 18s and 28s
#
#   Author : Rodrigo Rodrigues Domingues, minor modif by Eduarda Valério de Jesus 
#   Last updated 2024/01/17 
#
#   Script to perform genetic distance and phylogenetic tree between species based on 
#   molecular markers of the species represented in genetic sequence databases of the 
#   phylum Porifera
### --------------------------------------------------------------------------------------------------


#install and load packages
install.packages("adegenet")
install.packages("ape")
install.packages("phangorn")
install.packages("viridis")
library(adegenet)
library(ape)
library(phangorn)
library(viridis)


#-------------------#
#--------COI--------#
#-------------------#

#------------------------------genetic distance---------------------------------------
# A distância genética é baseada na diferença de uma sequência para outra, portanto,
# quanto mais diferenciação entre nucleotídeos, maior será a distância genética entre
# as espécies que estão sendo comparadas.

#inserir os dados -> arquivo deve estar no mesmo diretório que o script
dna_coi <- fasta2DNAbin(file="porifera_coi_alinhado_distgen.fasta") #arquivo fasta alinhado no MEGA sem grupo externo
dna_coi

D_coi <- dist.dna(dna_coi, model = "k80")
length(D_coi) #numero de comparações
D_coi

#transformar o dataframe como matriz 
temp_coi <- as.data.frame(as.matrix(D_coi))

#representar a distância genética em uma tabela (heatmap)
#cores mais claras - menores diferenças; cores mais escuras - maiores diferenças 
table.paint(temp_coi, cleg=0, clabel.row=.6, clabel.col=.6) 


#manipular a matrix - útil para a construção do gráfico abaixo
#esse código converte a matriz D em uma nova matriz chamada temp e inverte a ordem das colunas dessa matriz
temp1_coi <- t(as.matrix(D_coi))
temp1_coi <- temp1_coi[,ncol(temp1_coi):1]

#estilizar gráfico e salvar no diretório local
tiff(file = "distancia_genetica_porifera_coi.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1_coi, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="") #mudar tamanho da imagem em x e y 1:A, sendo A o número de sequências contidas 
axis(side = 2, at = 1:44, lab = rev(rownames(dna_coi)), las = 2, cex.axis = .5) #mudar numero de at para A 
axis(side = 3, at = 1:44, lab = rownames(dna_coi), las = 3, cex.axis = .5) #mudar numero de at para A 
dev.off()


#--------------------------------árvore filogenética---------------------------------- 

############ Maximum Likelihood-based ###############
# O método da Máxima Verossimilhança (Maximum Likelihood, ML) é um princípio 
# estatístico utilizado para estimar os parâmetros de um modelo estatístico. 
# A ideia central é encontrar os valores dos parâmetros que maximizam a probabilidade
# de observar os dados que foram realmente observados.
# resumindo: maximizar a probabilidade de observar os dados, dado os parâmetros do modelo.

#inserir os dados 
dna1_coi <- read.phyDat("porifera_coi_alinhado_arvgen.fasta", format = "fasta")
#arquivo fasta alinhado no MEGA com duas espécies de um grupo externo

# escolha de modelo evolutivo de substituição de base
mt_coi <- modelTest(dna1_coi)
fit_coi <- as.pml(mt_coi, "BIC")

fit_mt_coi <- pml_bb(mt_coi, control = pml.control(trace = 0))
fit_mt_coi
#model: TPM3u+G(4) 
#loglikelihood: -9271.962 
#unconstrained loglikelihood: -2841.588 
#Discrete gamma model
#Number of rate categories: 4 
#Shape parameter: 0.404116 

#Rate matrix:
#  a        c        g        t
#a 0.000000 1.525797 4.475491 1.000000
#c 1.525797 0.000000 1.525797 4.475491
#g 4.475491 1.525797 0.000000 1.000000
#t 1.000000 4.475491 1.000000 0.000000

#Base frequencies:  
#  a         c         g         t 
#0.2495950 0.1327952 0.1721350 0.4454748 

  
fitTIM_coi <- pml_bb(dna1_coi, model="TPM3u+G(4)") #colocar o modelo fornecido

# bootstrap
# O método de Bootstrap é uma técnica estatística que envolve a amostragem 
# repetida com substituição a partir de um conjunto de dados existente. 
# Ele é usado para estimar a variabilidade de uma estimativa estatística ou para
# avaliar a robustez de um procedimento estatístico.

bs_coi <- bootstrap.pml(fit_mt_coi, bs=1000, optNni=TRUE,
                    control = pml.control(trace = 0))
tree_coi <- plotBS(midpoint(fit_mt_coi$tree), bs_coi, p = 70, type="p", main="Phylogenetic tree of sponges COI")


### salvar a árvore
write.tree(fit_mt_coi$tree, "porifera_coi.tree")

#salvar árvore como imagem tiff
tiff("arvore_porifera_coi.tiff", width = 20, height = 10, units = 'in', res = 300)
plot(tree_coi) # Make plot
dev.off()


#-------------------#
#--------18S--------#
#-------------------#

#------------------------------genetic distance---------------------------------------

#inserir os dados -> arquivo deve estar no mesmo diretório que o script
dna_18s <- fasta2DNAbin(file="porifera_18s_alinhado_distgen.fasta") #arquivo fasta alinhado no MEGA sem grupo externo
dna_18s

D_18s <- dist.dna(dna_18s, model = "k80")
length(D_18s) #numero de comparações
D_18s

#transformar o dataframe como matriz 
temp_18s <- as.data.frame(as.matrix(D_18s))

#representar a distância genética em uma tabela (heatmap)
#cores mais claras - menores diferenças; cores mais escuras - maiores diferenças 
table.paint(temp_18s, cleg=0, clabel.row=.6, clabel.col=.6) 


#manipular a matrix - útil para a construção do gráfico abaixo
#esse código converte a matriz D em uma nova matriz chamada temp e inverte a ordem das colunas dessa matriz
temp1_18s <- t(as.matrix(D_18s))
temp1_18s <- temp1_18s[,ncol(temp1_18s):1]

#estilizar gráfico e salvar no diretório local
tiff(file = "distancia_genetica_porifera_18s.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:51, y = 1:51, temp1_18s, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="") #mudar tamanho da imagem em x e y 1:A, sendo A o número de sequências contidas 
axis(side = 2, at = 1:51, lab = rev(rownames(dna_18s)), las = 2, cex.axis = .5) #mudar numero de at para A 
axis(side = 3, at = 1:51, lab = rownames(dna_18s), las = 3, cex.axis = .5) #mudar numero de at para A 
dev.off()


#--------------------------------árvore filogenética---------------------------------- 

#inserir os dados 
dna1_18s <- read.phyDat("porifera_18s_alinhado_arvgen.fasta", format = "fasta")
#arquivo fasta alinhado no MEGA com duas espécies de um grupo externo

# escolha de modelo evolutivo de substituição de base
mt_18s <- modelTest(dna1_18s)
fit_18s <- as.pml(mt_18s, "BIC")

fit_mt_18s <- pml_bb(mt_18s, control = pml.control(trace = 0))
fit_mt_18s

fitTIM_18s <- pml_bb(dna1_18s, model="TrN+G(4)+I") #colocar o modelo fornecido

bs_18s <- bootstrap.pml(fit_mt_18s, bs=1000, optNni=TRUE,
                        control = pml.control(trace = 0))

tree_18s <- plotBS(midpoint(fit_mt_18s$tree), bs_18s, p = 70, type="p", main="Phylogenetic tree of sponges 18s")


### salvar a árvore
write.tree(fit_mt_18s$tree, "porifera_18s.tree")

#salvar árvore como imagem tiff
tiff("arvore_porifera_18s.tiff", width = 20, height = 10, units = 'in', res = 300)
plot(tree_18s) # Make plot
dev.off()

#-------------------#
#--------28S--------#
#-------------------#

#------------------------------genetic distance---------------------------------------

#inserir os dados -> arquivo deve estar no mesmo diretório que o script
dna_28s <- fasta2DNAbin(file="porifera_28s_alinhado_distgen.fasta") #arquivo fasta alinhado no MEGA sem grupo externo
dna_28s

D_28s <- dist.dna(dna_28s, model = "k80")
length(D_28s) #numero de comparações
D_28s

#transformar o dataframe como matriz 
temp_28s <- as.data.frame(as.matrix(D_28s))

#representar a distância genética em uma tabela (heatmap)
#cores mais claras - menores diferenças; cores mais escuras - maiores diferenças 
table.paint(temp_28s, cleg=0, clabel.row=.6, clabel.col=.6) 


#manipular a matrix - útil para a construção do gráfico abaixo
#esse código converte a matriz D em uma nova matriz chamada temp e inverte a ordem das colunas dessa matriz
temp1_28s <- t(as.matrix(D_28s))
temp1_28s <- temp1_28s[,ncol(temp1_28s):1]

#estilizar gráfico e salvar no diretório local
tiff(file = "distancia_genetica_porifera_28s.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="")
axis(side = 2, at = 1:43, lab = rev(rownames(dna)), las = 2, cex.axis = .5)
axis(side = 3, at = 1:43, lab = rownames(dna), las = 3, cex.axis = .5)
dev.off()
#talvez tenha que mudar at dos axis

#--------------------------------árvore filogenética---------------------------------- 

#inserir os dados 
dna1_28s <- read.phyDat("teste.fas", format = "fasta")
#arquivo fasta alinhado no MEGA com duas espécies de um grupo externo

# escolha de modelo evolutivo de substituição de base
mt_28s <- modelTest(dna1_28s)
fit_28s <- as.pml(mt_28s, "BIC")

fit_mt_28s <- pml_bb(mt_28s, control = pml.control(trace = 0))
fit_mt_28s

fitTIM_28s <- pml_bb(dna1_28s, model="TIM2+G(4)+I") #colocar o modelo fornecido

bs_28s <- bootstrap.pml(fit_mt_28s, bs=1000, optNni=TRUE,
                        control = pml.control(trace = 0))

plotBS(midpoint(fit_mt_28s$tree), bs_28s, p = 70, type="p", main="Phylogenetic tree of sponges 28s")


### salvar a árvore
tree_28s <- write.tree(fit_mt_28s$tree, "porifera_28s.tree")

#salvar árvore como imagem tiff
tiff("arvore_porifera_28s.tiff", width = 20, height = 10, units = 'in', res = 300)
plot(tree_28s) # Make plot
dev.off()

##EDIÇÃO DAS ÁRVORES 
#apos salvar as árvores é possível edita-las usando o programa FigTree.
# http://tree.bio.ed.ac.uk/software/figtree/