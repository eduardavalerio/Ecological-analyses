# DISTÂNCIA GENÉTICA E ÁRVORE FILOGENÉTICA 
# Script para performar distância genética e árvore filogenética entre espécies 
# a partir de sequências de marcadores moleculares.

#instalar e carregar os pacotes
install.packges ("adegenet")
install.packges ("ape")
install.packges ("phangorn")
install.packages("viridis")
library(adegenet)
library(ape)
library(phangorn)
library(viridis)

################################### distância genética ################################### 
# A distância genética é baseada na diferença de uma sequência para outra, portanto,
# quanto mais diferenciação entre nucleotídeos, maior será a distância genética entre
# as espécies que estão sendo comparadas.

#inserir os dados -> arquivo deve estar no mesmo diretório que o script
dna <- fasta2DNAbin(file="porifera_coi_alinhado_distgen.fasta") #arquivo fasta alinhado no MEGA sem grupo externo
dna

D <- dist.dna(dna, model = "k80")
length(D) #numero de comparações
D

#transformar o dataframe como matriz 
temp <- as.data.frame(as.matrix(D))

#representar a distância genética em uma tabela (heatmap)
#cores mais claras - menores diferenças; cores mais escuras - maiores diferenças 
table.paint(temp, cleg=0, clabel.row=.6, clabel.col=.6)


#manipular a matrix - útil para a construção do gráfico abaixo
#esse código converte a matriz D em uma nova matriz chamada temp e inverte a ordem das colunas dessa matriz
temp1 <- t(as.matrix(D))
temp1 <- temp1[,ncol(temp1):1]

#estilizar gráfico e salvar no diretório local
tiff(file = "distancia_genetica_porifera_COI.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="")
axis(side = 2, at = 1:44, lab = rev(rownames(dna)), las = 2, cex.axis = .5)
axis(side = 3, at = 1:44, lab = rownames(dna), las = 3, cex.axis = .5)
dev.off()


################################### árvore filogenética ################################### 


############ Maximum Likelihood-based ###############
# O método da Máxima Verossimilhança (Maximum Likelihood, ML) é um princípio 
# estatístico utilizado para estimar os parâmetros de um modelo estatístico. 
# A ideia central é encontrar os valores dos parâmetros que maximizam a probabilidade
# de observar os dados que foram realmente observados.
# resumindo: maximizar a probabilidade de observar os dados, dado os parâmetros do modelo.

#inserir os dados 
dna1 <- read.phyDat("porifera_coi_alinhado_arvgen_grupoexterno.fas", format = "fasta")
#arquivo fasta alinhado no MEGA com duas espécies de um grupo externo

# escolha de modelo evolutivo de substituição de base
mt <- modelTest(dna1)
fit <- as.pml(mt, "BIC")

fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
fit_mt
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

  
fitTIM <- pml_bb(dna1, model="TPM3u+G(4)") #colocar o modelo fornecido

# bootstrap
# O método de Bootstrap é uma técnica estatística que envolve a amostragem 
# repetida com substituição a partir de um conjunto de dados existente. 
# Ele é usado para estimar a variabilidade de uma estimativa estatística ou para
# avaliar a robustez de um procedimento estatístico.

bs <- bootstrap.pml(fit_mt, bs=1000, optNni=TRUE,
                    control = pml.control(trace = 0))

plotBS(midpoint(fit_mt$tree), bs, p = 70, type="p", main="Standard bootstrap")


### salvar a árvore
write.tree(fit_mt$tree, "porifera.tree")

##EDIÇÃO DA ÁRVORE 
#apos salvar a árvore é possível edita-la usando o programa FigTree.
# http://tree.bio.ed.ac.uk/software/figtree/
