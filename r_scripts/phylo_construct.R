setwd("/Users/mengwu/Documents/JaltPhylo/data")
library(phytools)
library(phangorn)


# A function to generate tree
plot_tre <-function(raxml.tree, name=NA, ultrametric=F, scale_axis=F, boostrap=T, rotate_node=NA){
  raxml.tree <- read.tree(raxml.tree)
  if (ultrametric) raxml.tree<-ultrametric(raxml.tree)
  for (x in 1:length(rotate_node)){
    if (!is.na(rotate_node[x])) raxml.tree <- rotate(raxml.tree, rotate_node[x])}
  plot(raxml.tree, main=name)
  if (boostrap) {
  for (i in 1:length(raxml.tree$node.label)){
    if (raxml.tree$node.label[i]=='') raxml.tree$node.label[i]=100}
  for (i in 1:length(raxml.tree$node.label)){
      if (as.numeric(raxml.tree$node.label[i])==100) raxml.tree$node.label[i]=''
      #else raxml.tree$node.label[i]='*'
      }
  #nodelabels(raxml.tree$node.label, adj = c(-0.2, 0.5), frame = "none", cex=0.6)}
  nodelabels(node=1:raxml.tree$Nnode+Ntip(raxml.tree),
           pie=cbind(as.numeric(raxml.tree$node.label),100-as.numeric(raxml.tree$node.label)),
           piecol=c("blue","white"),cex=0.2) }
  if (scale_axis) {
  h<-max(nodeHeights(raxml.tree))
  axis(side = 1, at = seq(0, h+0.003, 0.003))}
  }


# A function to generate ultrametric tree
ultrametric <- function(tree){
  ultrametric.tree <- nnls.tree(cophenetic(tree),tree,rooted=TRUE)
  return(ultrametric.tree)
}



## Plot each chromosome tree
pdf("../results/Chromosomes-trees.pdf", width=10, height=10) 
par(mar=c(0,0,2,0)+0.1)
par(mfrow=c(4,3))
plot_tre("RAxML_bipartitions.ch1", name="Chromosome 1", rotate_node=c(17,19,23,24,28,29))
plot_tre("RAxML_bipartitions.ch2", name="Chromosome 2", rotate_node=c(23,24,27,29))
plot_tre("RAxML_bipartitions.ch3", name="Chromosome 3", rotate_node=c(17,19,25,26,27,28,29))
plot_tre("RAxML_bipartitions.ch4", name="Chromosome 4", rotate_node=c(27,29))
plot_tre("RAxML_bipartitions.ch5", name="Chromosome 5", rotate_node=c(17,18,19,22,23,25,27,28))
plot_tre("RAxML_bipartitions.ch6", name="Chromosome 6", rotate_node=c(17,19,20,21,23,24,25,26,29))
plot_tre("RAxML_bipartitions.ch7", name="Chromosome 7", rotate_node=c(20,23))
plot_tre("RAxML_bipartitions.ch8", name="Chromosome 8", rotate_node=c(17,19,22,23,25))
plot_tre("RAxML_bipartitions.ch9", name="Chromosome 9", rotate_node=c(20,22,23,24,25,26,28,29))
plot_tre("RAxML_bipartitions.ch10", name="Chromosome 10", rotate_node=c(17,18,19,20,22,23,26,28))
plot_tre("RAxML_bipartitions.ch11", name="Chromosome 11", rotate_node=c(17,23,24,25,27,28))
plot_tre("RAxML_bipartitions.ch12", name="Chromosome 12", rotate_node=c(17,25,26,28))
dev.off()



# Plot whole-transcriptome concatenated tree
pdf("../results/Whole-transcriptome-tree.pdf", width=10, height=10) 
par(mar=c(0,0,1,0)+0.1)
par(mfrow=c(2,2))
plot_tre("RAxML_bipartitions.transcriptome", name="Concatenated Tree", rotate_node=c(18,20,21,23,29))
text(0.01, 3, "Bootstrap values\nsmaller than 100\nare shown in\ncorresponding node")

# Plot consensus tree 
tree <- read.tree("consensus_tree2.tre")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 17)
tree <- rotate(tree, 19)
tree <- rotate(tree, 23)
tree <- rotate(tree, 27)
tree <- rotate(tree, 29)
plot(tree, main="Consensus Tree")
tree$node.label[3]='84/74/74'
tree$node.label[4]='39/1/2'
tree$node.label[5]='62/49/44'
tree$node.label[6]='34/33/28'
tree$node.label[7]='14/4/7'
tree$node.label[8]='12/4/4'
tree$node.label[9]='8/0/1'
tree$node.label[10]='9/0/0'
tree$node.label[11]='19/0/9'
tree$node.label[12]='22/25/21'
tree$node.label[13]='22/0/0'
tree$node.label[14]='65/75/75'
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)
text(3, 3, "% supported/IC/TCA\nare shown on each node")

# Plot coalscence tree 
plot_tre("astral_tree.tre", name="Coalescence Tree", rotate_node=c(17,19,23,27,29))
text(3, 3, "Bootstrap values\nsmaller than 100\nare shown on\ncorresponding node")

# Plot BUCKy tree
plot_tre("bucky_tree.tre", name="BUCKy Concordance Tree", rotate_node=c(19,23,24,27,28,29))
text(3, 3, "Concordance factor (CF)\nare shown on each\ninternal node")
dev.off()



# Plot claudogram tree
pdf("../results/Claudogram-trees-50.pdf", width=6, height=6) 
par(mar=c(3,1,1,2)+0.1)
trees <- read.tree("bootstrap_50.tre")
ultrametric.trees<-list()
ultrametric.trees <- lapply(trees, ultrametric)
class(ultrametric.trees)<-"multiPhylo"
densiTree(ultrametric.trees, type = "cladogram", alpha = 0.006, consensus = tree, 
          optim = F, scaleX = T, col = "dark grey", width = 1, cex = 0.8)  ## boostrap_70
dev.off()

pdf("../results/Claudogram-trees-70.pdf", width=6, height=6) 
par(mar=c(3,1,1,2)+0.1)
trees <- read.tree("bootstrap_70.tre")
ultrametric.trees<-list()
ultrametric.trees <- lapply(trees, ultrametric)
class(ultrametric.trees)<-"multiPhylo"
densiTree(ultrametric.trees, type="cladogram", scaleX=T, width=1, cex=0.8, 
          alpha=0.07, col="dark grey", consensus = tree, optim = F)  ## boostrap_50
dev.off()



## Introgression pattern
dat <- read.csv("Dstat.csv", head=TRUE, sep=",")
dat$group <- factor(dat$species, levels = c("AUR","YUN","BIF","SIN","AIJ","UMB","GRA","INC","DEN"))
pdf("../results/ABBA_pattern.pdf", width=5, height=4) 
ggplot(dat, aes(x=group, y=dstat, colour=group))+geom_point(size = 2)+ 
  xlab("") + ylab(substitute(paste(italic(D)," statistic")))+theme(axis.title=element_text(size=16)) +
  theme_bw() + theme(legend.title=element_blank()) + theme(legend.position="none") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
dev.off()
dat.mod = lm(dstat ~ group, data = dat)
summary(dat.mod)

pdf("../results/common_ABBA.pdf") 
dat <- read.csv("ABBA_common.csv", head=TRUE, sep=",")
boxplot(dat, main="Number of common sites support ABBA in pairwise comparisons", 
        ylab="Number of Sites", cex.lab = 1.2, cex = 1.5) 
dev.off()


#### Extral trees

#tree <- read.tree("consensus_tree.tre")
#tree$edge.length<-NULL # get rid of edge lengths
#plot(tree, main="consensus tree")
#tree$node.label[3]='51/45/39'
#tree$node.label[4]='22/26/21'
#tree$node.label[5]='9/3/3'
#tree$node.label[6]='6/0/0'
#tree$node.label[7]='7/0/0'
#tree$node.label[7]='14/8/7'
#tree$node.label[8]='16/0/0'
#tree$node.label[9]='5/0/0'
#tree$node.label[10]='14/1/3'
#tree$node.label[11]='7/0/0'
#tree$node.label[12]='46/63/63'
#tree$node.label[13]='71/60/60'
#tree$node.label[14]='34/0/1'
#nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)

###
cat("owls((J. auriculata:0.00193666032429746886,((J. yungayensis:0.00270300477905288574,(((J. sinuosa:0.00183756164901645921,(J. umbellata:0.00227921246034051145,J. aijana:0.00211501518103496440)100:0.00038070730247107091)100:0.00017859566424201920,J. biflora:0.00212070108423424806)100:0.00021688734778193693,(J. incahausina:0.00168344686867075892,(J. dendroidea:0.00154424769359758573,J. grandibaccata:0.00173981671290508689)93:0.00030253985858484921)100:0.00056474523031534466)100:0.00028782623666485701)100:0.00037677388005408061,(J. calliantha:0.00181474582432030666,J. quipuscoae:0.00172340931249726696)100:0.00132258957044916413)100:0.00088510281238342564)100:0.00253812316389537187,(J. darcyana:0.00446646938706040612,(J. repandidentata:0.00420040000984493647,J. procumbens:0.00314202964408389577)100:0.00041219999319868386)100:0.00349825342819147901);", file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
is.ultrametric(tree.owls)
utree = chronos(tree.owls, lambda = 1, model = "correlated")
write.tree(utree, file = "example.trees")
cat(readLines("example.trees"), sep = "\n")
###

# Plot just the concatenated tree
cat("owls(((J. auriculata:0.00193666032429746886,((J. yungayensis:0.00270300477905288574,(((J. sinuosa:0.00183756164901645921,(J. umbellata:0.00227921246034051145,J. aijana:0.00211501518103496440)100:0.00038070730247107091)100:0.00017859566424201920,J. biflora:0.00212070108423424806)100:0.00021688734778193693,(J. incahausina:0.00168344686867075892,(J. dendroidea:0.00154424769359758573,J. grandiccata:0.00173981671290508689)93:0.00030253985858484921)100:0.00056474523031534466)100:0.00028782623666485701)100:0.00037677388005408061,(J. calliantha:0.00181474582432030666,J. quipuscoae:0.00172340931249726696)100:0.00132258957044916413)100:0.00088510281238342564)100:0.00253812316389537187,(J. darcyana:0.00446646938706040612,(J. repandidentata:0.00420040000984493647,J. procumbens:0.00314202964408389577)100:0.00041219999319868386)100:0.00349825342819147901):0.03542038372018924131,S. lycopersicum:0.03542038372018924131);", file = "ex.tre", sep = "\n")
tree.owls <- read.tree("ex.tre")
is.ultrametric(tree.owls)
utree = chronos(tree.owls, lambda = 1, model = "correlated")
write.tree(utree, file = "example.trees")
cat(readLines("example.trees"), sep = "\n")

# Plot BUCKy tree
pdf("../results/concatenated-tree.pdf", width=15, height=5) 
par(mar=c(3,1,1,2)+0.1)
plot_tre("bucky_tree.tre", ultrametric=T, scale_axis=T, rotate_node=c(19,23,24,27,28,29))
dev.off()
