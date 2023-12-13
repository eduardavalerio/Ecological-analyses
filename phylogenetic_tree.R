#Disciplina: Genetica aplicada à Oceanografia biologica
#Instituto Oceanografica da USP
#Professores: Dr. Rodrigo Domingues, Dr. Tito Lotufo
#Aula pratica 6 - filogenetica

#instale os pacotes
install.packges ("adegenet")
install.packges ("ape")
install.packges ("phangorn")
install.packages("viridis")

#carregue os pacotes
library(adegenet)
library(ape)
library(phangorn)
library(viridis)

#insira os dados
dna <- fasta2DNAbin(file="porifera_coi_alinhado_distgen.fasta")
dna

#distancia genetica
#qual a diferença de uma seq p/ outra -> quanto diferem de um nucleotídeo para outro
D <- dist.dna(dna, model = "k80")
length(D) #numero de comparacoes
D

#transforme o dataframe como matriz 
temp <- as.data.frame(as.matrix(D))

#represente a distancia genetica em uma tabela (heatmap)
table.paint(temp, cleg=0, clabel.row=.6, clabel.col=.6)
#cores mais claras - menores diferenças; cores mais escuras - maiores diferenças 

#manipule a matrix - util para a construcao do grafico abaixo
# esse código converte a matriz D em uma nova matriz chamada temp e inverte a ordem das colunas dessa matriz
temp1 <- t(as.matrix(D))
temp1 <- temp1[,ncol(temp1):1]

#grafico estilizado
tiff(file = "distancia_genetica_porifera_COI.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="")
axis(side = 2, at = 1:44, lab = rev(rownames(dna)), las = 2, cex.axis = .5)
axis(side = 3, at = 1:44, lab = rownames(dna), las = 3, cex.axis = .5)
dev.off()


################################### arvore filogenetica ########################


#arvore  baseada em distancia genetica 
# Os algoritmos UPGMA (Unweighted Pair Group Method with Arithmetic Mean) e NJ 
#(Neighbor-Joining) são metodos utilizados em analises filogeneticas para construir arvores filogeneticas
#a partir de uma matriz de distancias. Ele sao particularmente utilizado quando as
#distancias representam similaridades ou dissimilaridades entre grupos.

###### UPGMA ########
dm  <- dist.ml(dna)
treeUPGMA  <- upgma(dna)
plot(treeUPGMA, main="UPGMA")

############ Neighbor-Joining ################
treeNJ  <- NJ(dna)
plot(treeNJ, "unrooted", main="NJ")


############ Maximum Likelihood-based ###############
#O metodo da Maxima Verossimilhanca (Maximum Likelihood, ML) e um principio 
#estatistico utilizado para estimar os parametros de um modelo estatistico. 
#a ideia central e encontrar os valores dos parametros que maximizam a probabilidade
#de observar os dados que foram realmente observados.
#resumindo: maximizar a probabilidade de observar os dados, dados os parâmetros do modelo.

#insira o arquivo fasta
dna1 <- read.phyDat("porifera_coi_alinhado_arvgen_grupoexterno.fas", format = "fasta")

# escolha de modelo evolutivo de substituicao de base
mt <- modelTest(dna1)
fit <- as.pml(mt, "BIC")

fit_mt <- pml_bb(mt, control = pml.control(trace = 0))
fit_mt

fitTIM <- pml_bb(dna1, model="TPM3u+G(4)")

#bootstrap
#O metodo de Bootstrap e uma tecnica estatistica que envolve a amostragem 
#repetida com substituicao a partir de um conjunto de dados existente. 
#Ele e usado para estimar a variabilidade de uma estimativa estatistica ou para
#avaliar a robustez de um procedimento estatistico.

bs <- bootstrap.pml(fit_mt, bs=1000, optNni=TRUE,
                    control = pml.control(trace = 0))

plotBS(midpoint(fit_mt$tree), bs, p = 70, type="p", main="Standard bootstrap")


### salvar a arvore
write.tree(fit_mt$tree, "porifera.tree")

##EDIÇÃO DA ÁRVORE 
#apos salvar a arvore você pode edita-la usando o programa FigTree.
# http://tree.bio.ed.ac.uk/software/figtree/


################## parsimony ###############

#O principio da parcimonia, ou metodo da parcimonia, e uma abordagem na 
#construcao de arvores filogeneticas que busca encontrar a arvore mais simples
#e menos complexa possivel, minimizando o numero de mudancas evolutivas 
#(por exemplo, mutacoes geneticas) necessarias para explicar as diferencas 
#observadas entre os organismos. O princípio subjacente e que a explicacao mais
#simples geralmente e a mais plausivel.


parsimony(treeUPGMA, dna1)
parsimony(treeNJ, dna1)


treeRatchet  <- pratchet(dna1, trace = 0, minit=100)
parsimony(treeRatchet, dna1)


treeRatchet  <- acctran(treeRatchet, dna1)

treeRatchet  <- di2multi(treeRatchet)


if(inherits(treeRatchet, "multiPhylo")){
  treeRatchet <- unique(treeRatchet)
}


plotBS(midpoint(treeRatchet), type="phylogram")
add.scale.bar()
