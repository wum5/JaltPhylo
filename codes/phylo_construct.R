setwd("/Users/mengwu/Documents/Jalt_Phylo/results")
library(phytools)


plot_tre <-function(raxml.tree, name=NA){
  raxml.tree <- read.tree(raxml.tree)
  plot(raxml.tree, main=name)
  #nodelabels(raxml.tree$node.label, adj = c(1, 0), frame = "none")
  nodelabels(node=1:raxml.tree$Nnode+Ntip(raxml.tree),
           pie=cbind(as.numeric(raxml.tree$node.label),100-as.numeric(raxml.tree$node.label)),
           piecol=c("black","white"),cex=0.4)
}



pdf("Whole-transcriptome-tree.pdf") 
par(mar=c(0,0,1,0)+0.1)
par(mfrow=c(1,1))
plot_tre("RAxML_bipartitions.transcriptome")
dev.off()


pdf("Chromosomes-trees.pdf") 
par(mar=c(0,0,2,0)+0.1)
par(mfrow=c(4,3))
plot_tre("RAxML_bipartitions.ch1", "Chromosome 1")
plot_tre("RAxML_bipartitions.ch2", "Chromosome 2")
plot_tre("RAxML_bipartitions.ch3", "Chromosome 3")
plot_tre("RAxML_bipartitions.ch4", "Chromosome 4")
plot_tre("RAxML_bipartitions.ch5", "Chromosome 5")
plot_tre("RAxML_bipartitions.ch6", "Chromosome 6")
plot_tre("RAxML_bipartitions.ch7", "Chromosome 7")
plot_tre("RAxML_bipartitions.ch8", "Chromosome 8")
plot_tre("RAxML_bipartitions.ch9", "Chromosome 9")
plot_tre("RAxML_bipartitions.ch10", "Chromosome 10")
plot_tre("RAxML_bipartitions.ch11", "Chromosome 11")
plot_tre("RAxML_bipartitions.ch12", "Chromosome 12")
dev.off()

