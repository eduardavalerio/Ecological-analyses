################################################################################    
#       S1: R codes for phylogenetic analyses by using FASTA format file

#	A workflow with RStudio: phylogenetic analyses and 
#	      visualizations using mtDNA-Cytb sequences
#		
#	      Emine TOPARSLAN, Kemal KARABAG, Ugur BILGE
#			       June 2020

#R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
#Copyright (C) 2020 The R Foundation for Statistical Computing
###############################################################################
#  GENERAL INSTRUCTIONS

#This workflow is planned according to data sets without missing data and gaps.
#Before starting, please check your sequences that have any gaps or missing and fix them.


#The following packages must be installed to perform all the analyses outlined in the article. 
#Use these commands if you have not already installed them. 

#If your samples belonging to populations or groups, 
#you must modify the name of samples using  "_" between population/group name and number.  
#This naming method allows extracting unique names as population names from sample names with the help of a short command. 
#Thus, the name of the population in all the analyzes do not need to be entered again. 

#If you can't change the name of the samples, please run the commands following. 
#After a few steps, you will see the example of name change codes. Please, follow the commands below.
############################################################################################  

#set working directory and input file name below
#setwd("C:/Users/ugur/desktop/ileribiyoinformatik")
setwd("C:/Users/usr/Documents")

sink("18s-tree-anthozoa.txt")
###----------- INSTALL PACKAGES
#If you have already them, you can skip this part (InstallPackages is set to FALSE by default)
InstallPackages = T

if (InstallPackages) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) 
    install.packages("BiocManager")
  BiocManager::install("msa")
  
  install.packages("adegenet")
  install.packages("ape")
  install.packages("ggtree")
  install.packages("ggplot2")
  install.packages("ips")
  install.packages("bios2mds")
  install.packages("haplotypes")
  install.packages("pegas")
  install.packages("phytools")
  install.packages("stats")
  install.packages("treeio")
  install.packages("geiger")
}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
library(adegenet)
library(ape)
library(ggtree)
library(ggplot2)
library(stats)
library(ips)
library(msa)
library(geiger)
library(reshape)
library(pegas)
library(phytools)
library(stats)

##############################################################
# WARNING
#We would like to warn those who will use the Align command that the alignment process can take 
#only a minute or a few hours depending on the length of your sequence and the number of samples or strong of your devices. 
#For example, 740 samples that length is 1226 bp may take approximately 4 hours for alignment.
#We aligned FASTA formatted file using msa() command and ordered the names in the "nbin" object, 
#when coloring haplotype networks and renaming samples for grouping.
#Note: "S2_Appendix.fas" file have been non-aligned before.

#Before starting the analysis, it is necessary to align the FASTA format sequences. 
#If your sequences are already aligned, please proceed from "READING AND PLOTTING OF ALIGNMENT-NJ-MSAPLOT" section.
#Even so, if you want to alignment again and export your FASTA file, you can follow instructions as below.
#If your sequences are not aligned, please run the commands following.
#If your samples already aligned, you can also set "AlignNeeded = FALSE" below.
##############################################################

#########################################################################
###########   INPUT FASTA FILES AND ALIGN THE SEQUENCES #################
#########################################################################

###---------- input file

fname = "anthozoa_18s_same_length_small.fasta" 

###---------- The program reads fasta file and aligns it

AlignNeeded = TRUE 

if (AlignNeeded) {
  
file <- readDNAStringSet(fname)#for reading multiple DNA sequences from msa package
file
  
###---------- multiple sequence alignment from msa package

cb<- msa(file)    
cb
  
###----------- CONVERTING ALIGN FILE TO FASTA FILE    
  
cv<-msaConvert(cb, type=c("bios2mds::align"))
  
###----------- EXPORTING  ALIGNED FASTA FILE              
  
  library(bios2mds) # for exporting fasta file
  
  export.fasta(cv, outfile = "outfile.fas", ncol(cb), open = "w")
}

#After align and export, please check your sequence that has whether gaps or missing in the "outfile.fas".
#If sequences have gaps or are missing, please fix them and call again using "fasta2DNAbin(fname/"outfile.fas")" as below.
#If your sequences have not gaps or are not missing, you can directly contuniue to "nbin<-as.DNAbin(cb)" command as below.

nbin<-as.DNAbin(cb) #read aligned data from cb above


################################################################################
############ READING AND PLOTTING OF ALIGNMENT-NJ-MSAPLOT (kimura tree) ######## 
################################################################################

#If you don't need to align your file or you get repaired "outfile.fas", 
#you can use "readFASTA2DNAbin(fname)" command to call them to RStudio console as below.

nbin<-fasta2DNAbin(fname)# read your previously aligned FASTA file name.


#After align you must use trimEnds() command as below. 
#If you don't need to trim you can follow the next command.

TRIM = FALSE       # Already trimmed sequence is assumed

if (TRIM) {
  nbin<-trimEnds(nbin)#trimming of sequences ends
}

############################################# 
#WARNING    
#Before starting the analysis, you should be sure that the names of samples like "name_1; name_2, ...".
#If you couldn't change it before, you can fastly change them as below.
#Firstly, you can write names of aligned sequences to see the name list using the "labels(nbin)" command, 
#then you can rename them as below;
#A, B, and C letters represent names of populations or groups.
#rownames(nbin)[1:10]=paste("A_",1:10, sep="")
#rownames(nbin)[11:20]=paste("B_",1:10, sep="")
#rownames(nbin)[21:30]=paste("C_",1:10, sep="")....
#After changing the name of samples, you can follow codes as below.
#############################################

###---------- converting DNAbin to alignment format

an<-as.alignment(nbin)  

###---------- converting alignment to matrix

nm<-as.matrix(an) 

###---------- extraction of the sample names

nbinmat<-as.matrix(labels(nbin)) 
nbin
class(nbin)

###---------- computing distance by ape package with K80 model derived by Kimura (1980)

dnbin<-dist.dna(nbin, model = "K80") 
tree<-nj(dnbin)
tiff(file = "tree_antozoa_16s_nj.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
ggt<-ggtree(tree, cex = 0.8, aes(color=branch.length))+
scale_color_continuous(high='grey',low='black')+
geom_tiplab(align=TRUE, size=2)+
geom_treescale(y = - 5, color = "black", fontsize = 4)
dev.off()
ggt

#The shared S2_Appendix.fas dataset is not contained missing value, different letter (N, R, K etc.) or gaps.
#If your sequences have gaps (-) or N, Y, K, R letter (IUPAC nucleotide code), 
#you have to add additional color in "color argument" below.
#For example, "rep("green",1,)" for gaps, "rep("pink",1,)" for N letter etc.
#However, you have to delete or fix the gaps to be able to do other analysis.
#If you have big data set, you can modified "width" and "height" arguments below for better images.

tiff(file = "tree_antozoa_16s_nj.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
njmsaplot<-msaplot(ggt, nbin, offset = 0.2, width=1.5, height = 0.6, 
            color = c(rep("grey", 1), rep("green", 1), rep("blue", 1), rep("black", 1), rep("yellow", 1),
            rep("red", 1)))
njmsaplot
dev.off()


###---------- neighbor-joining tree estimation of Saitou and Nei (1987)

tiff(file = "porifera_nj_coi_distance.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
njdistree<-ggtree(tree,layout = 'circular', branch.length='branch.length', aes(color=branch.length), lwd = 0.5)+xlim(-2, NA)+
  geom_tiplab(names(nbin), size = 3, offset=0.002)+scale_color_continuous(high='lightskyblue1',low='black')+
  geom_treescale(x=-2, color = "black", fontsize = 3, offset = 9) 
dev.off()
njdistree

############################################################################
####################   DISTANCE MATRIX OF HAPLOTYPES    ####################
############################################################################

library(pegas)

D_coi <- dist.dna(nbin, model = "k80")
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
tiff(file = "distancia_genetica_anthozoa_coi.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
par(mar = c(2,7,7,2))
image(x = 1:21, y = 1:21, temp1_coi, keep.dendro=TRUE, symm=TRUE, col = rev(viridisLite::viridis(100)),
      xaxt = "n", yaxt = "n", xlab="",ylab="") #change the image size in x and y 1:A, with A being the number of sequences 
axis(side = 2, at = 1:21, lab = rev(rownames(temp1_coi)), las = 2, cex.axis = .5) #change the at number for A 
axis(side = 3, at = 1:21, lab = rownames(temp1_coi), las = 3, cex.axis = .5) #change the at number for A 
dev.off()

########## calculate the Hamming distance from pegas package
dh<-dist.hamming(nbin) 
dh
dhm<-as.matrix(dh)
write.table(dhm,file="dhm.txt", quote=FALSE, sep="\t")

###---------- heatmap figure 

tiff(file = "heatmap.tiff", width=10, heigh=8, 
     unit="in",res = 600, compression="lzw")
heatmap(dhm,scale="none", col= heat.colors(100),keep.dendro=TRUE, symm=TRUE) #stats package
dev.off()

sink()
