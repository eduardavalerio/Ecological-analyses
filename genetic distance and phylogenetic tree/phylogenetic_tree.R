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
install.packges ("adegenet")
install.packges ("ape")
install.packges ("phangorn")
install.packages("viridis")
library(adegenet)
library(ape)
library(phangorn)
library(viridis)


#-------------------#
#--------COI--------#
#-------------------#

#------------------------------genetic distance---------------------------------------
# The genetic distance is based on the difference from one sequence to another.
# The more differentiation between nucleotides, the greater the genetic distance between
# the species being compared.

#data
dna_coi <- fasta2DNAbin(file="porifera_coi_alinhado_distgen.fasta") #fasta file aligned on MEGA without outgroup
dna_coi

D_coi <- dist.dna(dna_coi, model = "k80")
length(D_coi) 
D_coi

#transform dataframe as matrix
temp_coi <- as.data.frame(as.matrix(D_coi))

#represent the genetic distance in a table (heatmap)
#light colors- minor differences; dark colors - biggest differences 
table.paint(temp_coi, cleg=0, clabel.row=.6, clabel.col=.6) 


temp1_coi <- t(as.matrix(D_coi))
temp1_coi <- temp1_coi[,ncol(temp1_coi):1]

#save graph on local directory
tiff(file = "distancia_genetica_porifera_coi.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1_coi, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="") #change the image size in x and y 1:A, with A being the number of sequences 
axis(side = 2, at = 1:44, lab = rev(rownames(dna_coi)), las = 2, cex.axis = .5) #change the at number for A 
axis(side = 3, at = 1:44, lab = rownames(dna_coi), las = 3, cex.axis = .5) #change the at number for A 
dev.off()


#--------------------------------phylogenetic tree---------------------------------- 

############ Maximum Likelihood-based ###############
# The Maximum Likelihood (ML) method is a statistical principle used to estimate the 
# parameters of a statistical model.
# The central idea is to find the parameter values ​​that maximize the probability
# of observing the data that were actually observed.
# In sume: The probability of maximizing the data, given the model parameters.


#data
dna1_coi <- read.phyDat("porifera_coi_alinhado_arvgen_grupoexterno.fasta", format = "fasta")
#fasta file aligned on MEGA with two species  from an outgroup.

# evolutionary model of base substitution
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

  
fitTIM_coi <- pml_bb(dna1_coi, model="TPM3u+G(4)") 

# bootstrap
# The Bootstrap method is a statistical technique that involves sampling
# repeated with replacement from an existing data set.
# It is used to estimate the variability of a statistical estimate or to
# evaluate the robustness of a statistical procedure.

bs_coi <- bootstrap.pml(fit_mt_coi, bs=1000, optNni=TRUE,
                    control = pml.control(trace = 0))

plot_coi <- plotBS(midpoint(fit_mt_coi$tree), bs_coi, p = 70, type="p", main="Standard bootstrap")


### save tree
write.tree(fit_mt_coi$tree, "porifera_coi.tree")

#save tree as tiff image
tiff("arvore_porifera_coi.tiff", width = 20, height = 10, units = 'in', res = 300)
plot(plot_coi) # Make plot
dev.off()


#-------------------#
#--------18S--------#
#-------------------#

#------------------------------genetic distance---------------------------------------

#data
dna_18s <- fasta2DNAbin(file="porifera_18s_alinhado_distgen.fasta") 
dna_18s

D_18s <- dist.dna(dna_18s, model = "k80")
length(D_18s)
D_18s

temp_18s <- as.data.frame(as.matrix(D_18s))

table.paint(temp_18s, cleg=0, clabel.row=.6, clabel.col=.6) 


temp1_18s <- t(as.matrix(D_18s))
temp1_18s <- temp1_18s[,ncol(temp1_18s):1]

#save graph
tiff(file = "distancia_genetica_porifera_18s.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1_18s, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="") #change the image size in x and y 1:A, with A being the number of sequences
axis(side = 2, at = 1:43, lab = rev(rownames(dna_18s)), las = 2, cex.axis = .5) #change the at number for A 
axis(side = 3, at = 1:43, lab = rownames(dna_18s), las = 3, cex.axis = .5) #change the at number for A
dev.off()


#--------------------------------phylogenetic tree---------------------------------- 

#data
dna1_18s <- read.phyDat("porifera_18s_alinhado_arvgen_grupoexterno.fasta", format = "fasta")

# evolutionary model of base substitution
mt_18s <- modelTest(dna1_18s)
fit_18s <- as.pml(mt_18s, "BIC")

fit_mt_18s <- pml_bb(mt_18s, control = pml.control(trace = 0))
fit_mt_18s

fitTIM_18s <- pml_bb(dna1_18s, model="TrN+G(4)+I")

bs_18s <- bootstrap.pml(fit_mt_18s, bs=1000, optNni=TRUE,
                        control = pml.control(trace = 0))

plotBS(midpoint(fit_mt_18s$tree), bs_18s, p = 70, type="p", main="Standard bootstrap")


### save tree
tree_18s <- write.tree(fit_mt_18s$tree, "porifera_18s.tree")

#save tree as tiff image
tiff("arvore_porifera_18s.tiff", width = 20, height = 10, units = 'in', res = 300)
plot(tree_18s) # Make plot
dev.off()


#-------------------#
#--------28S--------#
#-------------------#

#------------------------------genetic distance---------------------------------------

#data
dna_28s <- fasta2DNAbin(file="porifera_28s_alinhado_distgen.fasta")
dna_28s

D_28s <- dist.dna(dna_28s, model = "k80")
length(D_28s)
D_28s

temp_28s <- as.data.frame(as.matrix(D_28s))

table.paint(temp_28s, cleg=0, clabel.row=.6, clabel.col=.6) 


temp1_28s <- t(as.matrix(D_28s))
temp1_28s <- temp1_28s[,ncol(temp1_28s):1]

#save graph
tiff(file = "distancia_genetica_porifera_28s.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:44, y = 1:44, temp1, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="") #change the image size in x and y 1:A, with A being the number of sequences
axis(side = 2, at = 1:43, lab = rev(rownames(dna)), las = 2, cex.axis = .5) #change the at number for A
axis(side = 3, at = 1:43, lab = rownames(dna), las = 3, cex.axis = .5) #change the at number for A
dev.off()

#--------------------------------phylogenetic tree---------------------------------- 

#data 
dna1_28s <- read.phyDat("porifera_28s_alinhado_arvgen_grupoexterno.fasta", format = "fasta")


# evolutionary model of base substitution
mt_28s <- modelTest(dna1_28s)
fit_28s <- as.pml(mt_28s, "BIC")

fit_mt_28s <- pml_bb(mt_28s, control = pml.control(trace = 0))
fit_mt_28s

fitTIM_28s <- pml_bb(dna1_28s, model="TrN+G(4)+I") 

bs_28s <- bootstrap.pml(fit_mt_28s, bs=1000, optNni=TRUE,
                        control = pml.control(trace = 0))

plotBS(midpoint(fit_mt_28s$tree), bs_28s, p = 70, type="p", main="Standard bootstrap")


###save tree
tree_28s <- write.tree(fit_mt_28s$tree, "porifera_28s.tree")

#save tree as tiff image
tiff("arvore_porifera_28s.tiff", width = 20, height = 10, units = 'in', res = 300)
plot(tree_28s) # Make plot
dev.off()

# It is possible to edit the phylogenetic trees on FigTree
# http://tree.bio.ed.ac.uk/software/figtree/
