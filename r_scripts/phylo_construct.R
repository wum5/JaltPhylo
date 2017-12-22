setwd("/Users/mengwu/Documents/Research/JaltPhylo/updated_scripts/r_scripts")
library(ggplot2)
require(gridExtra)
source("base_functions.R")



### Generate Figure 2
## Plot the tree in Figure 2A
pdf("fig2-tree.pdf", width=6, height=8) 
par(mar=c(3,1,1,2)+0.1)
plot_tre("../tree/Fig2A_tree.tre", ultrametric=T, tipdrop="S.lycopersicum", 
         scale_axis=T, rotate_node=c(18,20,21,23,29), boostrap_pie=T)
dev.off()

# Plot claudogram tree in Figure 2B
tree <- read.tree("../tree/RAxML_bipartitions.transcriptome")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 18)
tree <- rotate(tree, 20)
tree <- rotate(tree, 21)
tree <- rotate(tree, 23)
tree <- rotate(tree, 29)

pdf("Claudogram-trees-70.pdf", width=6, height=6) 
par(mar=c(3,1,1,2)+0.1)
trees70 <- read.tree("../tree/trees_bootstrap_70")
ultrametric.trees70 <- list()
ultrametric.trees70 <- lapply(trees70, ultrametric)
class(ultrametric.trees70) <- "multiPhylo"
densiTree(ultrametric.trees70, type="cladogram", scaleX=T, width=1, cex=0.8, 
          alpha=0.05, col="dark grey", consensus = tree, optim = F)  ## boostrap_70

## plot the concatenated tree into the claudogram
tree <- read.tree("../tree/RAxML_bipartitions.transcriptome")
ultrametric.tre <- ultrametric(tree)
ultrametric.tre <- rotate(ultrametric.tre, 18)
ultrametric.tre <- rotate(ultrametric.tre, 20)
ultrametric.tre <- rotate(ultrametric.tre, 21)
ultrametric.tre <- rotate(ultrametric.tre, 23)
ultrametric.tre <- rotate(ultrametric.tre, 29)
par(new=T)
plot(ultrametric.tre, type="cladogram", show.tip.label=FALSE)
dev.off()



### Generate Figure S3
pdf("Figure_S3.pdf", width=10, height=10) 
par(mar=c(0,0,1,0)+0.1)
par(mfrow=c(2,2))

# Plot boostrap50-gene concatenated tree
plot_tre("../tree/RAxML_bipartitions.seqs_bootstrap50", name="RAxML Concatenated Tree", rotate_node=c(23,25))
text(0.01, 3, "Bootstrap values\nsmaller than 100\nare shown in\ncorresponding node")

# Plot consensus tree 
tree <- read.tree("../tree/consense_bootstrap50")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 21)
tree <- rotate(tree, 22)
tree <- rotate(tree, 24)
tree <- rotate(tree, 25)
plot(tree, main="Majority Consensus Tree")
tree$node.label[17]='64/0.50/0.43'
tree$node.label[18]='36/0.32/0.26'
tree$node.label[19]='16/0.06/0.06'
tree$node.label[20]='14/0.04/0.05'
tree$node.label[21]='10/0.00/0.02'
tree$node.label[22]='10/-0.01/-0.02'
tree$node.label[23]='18/0.00/0.08'
tree$node.label[24]='25/0.28/0.22'
tree$node.label[25]='22/0.00/0.00'
tree$node.label[26]='68/0.73/0.73'
tree$node.label[27]='86/0.78/0.78'
tree$node.label[28]='40/0.00/0.03'
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)
text(3, 3, "% supported/IC/TCA\nare shown on each node")

# Plot quartet-based tree 
plot_tre("../tree/astral_tree.tre", name="ASTRAL Quartet-based Tree", 
         outgroup="S.lycopersicum", rotate_node=c(18,22,26,27))
text(3, 3, "Bootstrap values\nsmaller than 100\nare shown on\ncorresponding node")

# Plot BUCKy tree
plot_tre("../tree/bucky_tree.tre", name="BUCKy Concordance Tree", 
         branch_length=F, rotate_node=c(19,23,24,27,28,29))
text(3, 3, "Concordance factor (CF)\nare shown on each\ninternal node")
dev.off()




## Plot each chromosome tree
pdf("Chromosomes-trees.pdf", width=10, height=10) 
par(mar=c(0,0,2,0)+0.1)
par(mfrow=c(4,3))
plot_tre("../tree/RAxML_bipartitions.ch1", name="Chromosome 1", 
         tipdrop="S.lycopersicum", rotate_node=c(21,22,24,25,26,27,29))
plot_tre("../tree/RAxML_bipartitions.ch2", name="Chromosome 2", 
         tipdrop="S.lycopersicum", rotate_node=c(17,19,20,21,23,25,26,27,28))
plot_tre("../tree/RAxML_bipartitions.ch3", name="Chromosome 3", 
         tipdrop="S.lycopersicum", rotate_node=c(23,25,26,27,29))
plot_tre("../tree/RAxML_bipartitions.ch4", name="Chromosome 4", 
         tipdrop="S.lycopersicum", rotate_node=c(17,22,23,24,26))
plot_tre("../tree/RAxML_bipartitions.ch5", name="Chromosome 5", 
         tipdrop="S.lycopersicum", rotate_node=c(22,25,27,28))
plot_tre("../tree/RAxML_bipartitions.ch6", name="Chromosome 6", 
         tipdrop="S.lycopersicum", rotate_node=c(17,19,20,21,23,24,25,26,29))
plot_tre("../tree/RAxML_bipartitions.ch7", name="Chromosome 7", 
         tipdrop="S.lycopersicum", rotate_node=c(20,22,23,24,26,27,29))
plot_tre("../tree/RAxML_bipartitions.ch8", name="Chromosome 8", 
         tipdrop="S.lycopersicum", rotate_node=c(17,19,20,21,23,25,26,27,28))
plot_tre("../tree/RAxML_bipartitions.ch9", name="Chromosome 9", 
         tipdrop="S.lycopersicum", rotate_node=c(22,23,24,25,27))
plot_tre("../tree/RAxML_bipartitions.ch10", name="Chromosome 10", 
         tipdrop="S.lycopersicum", rotate_node=c(21,22,23,24,26,27,29))
plot_tre("../tree/RAxML_bipartitions.ch11", name="Chromosome 11", 
         tipdrop="S.lycopersicum", rotate_node=c(22,24,25,27,29))
plot_tre("../tree/RAxML_bipartitions.ch12", name="Chromosome 12", 
         tipdrop="S.lycopersicum", rotate_node=c(19,20,22,25,27,29))
dev.off()



## Introgression pattern
dat <- read.csv("Dstat.csv", head=TRUE, sep=",")
dat$group <- factor(dat$species, levels = c("AUR","YUN","BIF","SIN","AIJ","UMB","GRA","INC","DEN"))
pdf("ABBA_pattern.pdf", width=5, height=4) 
ggplot(dat, aes(x=group, y=dstat, colour=group))+geom_point(size = 2)+ 
  xlab("") + ylab(substitute(paste(italic(D)," statistic")))+theme(axis.title=element_text(size=16)) +
  theme_bw() + theme(legend.title=element_blank()) + theme(legend.position="none") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_blank())
dev.off()
dat.mod = lm(dstat ~ group, data = dat)
summary(dat.mod)

pdf("common_ABBA.pdf") 
dat <- read.csv("ABBA_common.csv", head=TRUE, sep=",")
boxplot(dat, main="Number of common sites support ABBA in pairwise comparisons", 
        ylab="Number of Sites", cex.lab = 1.2, cex = 1.5) 
dev.off()



#### branch lengths vs. branch concordance
pdf("branch_lengths_vs_concordance.pdf", width=11, height=5) 

mytree <- read.tree("../tree/RAxML_bipartitions.transcriptome")
ultra_tre <- ultrametric(mytree)
rotate_node=c(18,20,23,25,26,28,29)
for (x in 1:length(rotate_node)){
  if (!is.na(rotate_node[x])) ultra_tre <- rotate(ultra_tre, rotate_node[x])}
#plot(ultra_tre)
edgelabels(round(ultra_tre$edge.length, digits = 2),frame="none",bg=NA,adj = c(0.5, 0))
brchLen <- c(0.16,1.21,0.54,1.32,0.3,0.13,0.1,0.08,0.06,0.12,0.12,0.25)
concord <- c(0.395,0.873,0.670,0.698,0.387,0.195,0.159,0.114,0.120,0.201,0.211,0.282)
ica_val <- c(0.03,0.78,0.43,0.73,0.26,0.06,0.05,0.02,-0.02,0.08,0.00,0.22)
subclad <- c("black","black","green","grey","grey", rep("orange",7))
df <- data.frame(brchLen,concord,ica_val,subclad)

cor.test(df$brchLen,df$concord) 
cor.test(df$brchLen,df$ica_val) 

plot1 <- qplot(brchLen, concord, data=df, color=subclad, size=1, xlim = c(0,1.5), ylim = c(0,1), 
      xlab="Internal Branch Length", ylab="Internode Concordance") + theme_bw() + 
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text=element_text(size=12), 
        axis.title=element_text(size=16)) +
  scale_color_manual(breaks = subclad, values= c("purple","dark green","grey","orange")) +
  annotate("text", x = 1.1, y=0.2, label = "Pearspn correlation: r=0.908\nP-value: 4.538e-05")

plot2 <- qplot(brchLen, ica_val, data=df, color=subclad, size=1, xlim = c(0,1.5), ylim = c(-0.1,1), 
      xlab="Internal Branch Length", ylab="Internode Certainty All (ICA)") + theme_bw() + 
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.text=element_text(size=12), 
        axis.title=element_text(size=16)) +
  scale_color_manual(breaks = subclad, values= c("purple","dark green","grey","orange")) +
  annotate("text", x = 1.1, y=0.2, label = "Pearspn correlation: r=0.986\nP-value: 6.327e-08")

grid.arrange(plot1, plot2, ncol=2)

dev.off()



#### comparing two alternative topologies
pdf("alt_topology.pdf", width=10, height=5) 
par(mfrow=c(1,2))
tree1 <- read.tree(text = "((J. umbellata,J. aijana),(J. sinuosa,J. biflora));")
plot(tree1, edge.width = 2, cex=1.5)
tree2 <- read.tree(text = "(((J. umbellata,J. aijana),J. sinuosa),J. biflora);")
plot(tree2, edge.width = 2, cex=1.5)
dev.off()



## plot consense trees with different bootstrap cutoff
pdf("consensus-trees-bootstrap.pdf", width=8, height=8) 
par(mfrow=c(2,2))
par(mar=c(0,0,2,0)+0.1)

tree <- read.tree("../tree/consense_bootstrap50")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 21)
tree <- rotate(tree, 22)
tree <- rotate(tree, 24)
tree <- rotate(tree, 25)
plot(tree, main="50% dataset (1190 gene trees)")
tree$node.label[17]='64/0.50/0.43'
tree$node.label[18]='36/0.32/0.26'
tree$node.label[19]='16/0.06/0.06'
tree$node.label[20]='14/0.04/0.05'
tree$node.label[21]='10/0.00/0.02'
tree$node.label[22]='10/-0.01/-0.02'
tree$node.label[23]='18/0.00/0.08'
tree$node.label[24]='25/0.28/0.22'
tree$node.label[25]='22/0.00/0.00'
tree$node.label[26]='68/0.73/0.73'
tree$node.label[27]='86/0.78/0.78'
tree$node.label[28]='40/0.00/0.03'
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)

tree <- read.tree("../tree/consense_bootstrap60")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 19)
tree <- rotate(tree, 22)
tree <- rotate(tree, 23)
tree <- rotate(tree, 25)
tree <- rotate(tree, 26)
tree$node.label[17]='69/0.55/0.48'
tree$node.label[18]='38/0.32/0.25'
tree$node.label[19]='76/0.79/0.79'
tree$node.label[20]='18/0.05/0.05'
tree$node.label[21]='16/0.10/0.10'
tree$node.label[22]='11/0.01/0.03'
tree$node.label[23]='11/-0.01/-0.02'
tree$node.label[24]='19/0.00/0.07'
tree$node.label[25]='28/0.34/0.29'
tree$node.label[26]='24/0.00/0.00'
tree$node.label[27]='88/0.79/0.79'
tree$node.label[28]='42/0.01/0.03'
plot(tree, main="60% dataset (628 gene trees)")
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)

tree <- read.tree("../tree/consense_bootstrap70")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 20)
tree <- rotate(tree, 21)
tree <- rotate(tree, 22)
tree <- rotate(tree, 23)
tree <- rotate(tree, 25)
tree <- rotate(tree, 28)
tree$node.label[17]='75/0.57/0.54'
tree$node.label[18]='43/0.35/0.31'
tree$node.label[19]='17/0.02/0.02'
tree$node.label[20]='18/0.10/0.07'
tree$node.label[21]='14/0.00/0.05'
tree$node.label[22]='11/-0.03/-0.05'
tree$node.label[23]='25/0.02/0.09'
tree$node.label[24]='31/0.33/0.28'
tree$node.label[25]='25/0.00/0.00'
tree$node.label[26]='83/0.81/0.81'
tree$node.label[27]='93/0.86/0.86'
tree$node.label[28]='43/0.00/0.05'
plot(tree, main="70% dataset (207 gene trees)")
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)

tree <- read.tree("../tree/consense_bootstrap80")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 28)
tree <- rotate(tree, 18)
tree <- rotate(tree, 19)
tree <- rotate(tree, 20)
tree <- rotate(tree, 22)
tree <- rotate(tree, 26)
tree <- rotate(tree, 27)
tree$node.label[17]='86/0.67/0.67'
tree$node.label[18]='42/0.35/0.33'
tree$node.label[19]='94/0.81/0.81'
tree$node.label[20]='14/0.05/0.04'
tree$node.label[21]='14/0.01/0.03'
tree$node.label[22]='14/-0.04/-0.08'
tree$node.label[23]='22/0.00/0.04'
tree$node.label[24]='14/-0.04/-0.04'
tree$node.label[25]='31/0.16/0.14'
tree$node.label[26]='28/0.00/0.00'
tree$node.label[27]='97/0.82/0.82'
tree$node.label[28]='47/0.00/0.09'
plot(tree, main="80% dataset (36 gene trees)")
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)

dev.off()


# Plot whole-transcriptome concatenated tree (50% dataset)
pdf("50bootstrap-transcriptome_tree.pdf", width=10, height=10) 
par(mar=c(0,0,1,0)+0.1)
par(mfrow=c(2,2))
plot_tre("../tree/RAxML_bipartitions.seqs_bootstrap50", 
         name="RAxML Concatenated Tree (50% dataset)", 
         rotate_node=c(23,25))
text(0.01, 3, "Bootstrap values\nsmaller than 100\n
     are shown in\ncorresponding node")

# Plot consensus tree 
tree <- read.tree("../tree/consense_bootstrap50")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 21)
tree <- rotate(tree, 22)
tree <- rotate(tree, 24)
tree <- rotate(tree, 25)
plot(tree, main="50% dataset (1190 gene trees)")
tree$node.label[17]='64/0.50/0.43'
tree$node.label[18]='36/0.32/0.26'
tree$node.label[19]='16/0.06/0.06'
tree$node.label[20]='14/0.04/0.05'
tree$node.label[21]='10/0.00/0.02'
tree$node.label[22]='10/-0.01/-0.02'
tree$node.label[23]='18/0.00/0.08'
tree$node.label[24]='25/0.28/0.22'
tree$node.label[25]='22/0.00/0.00'
tree$node.label[26]='68/0.73/0.73'
tree$node.label[27]='86/0.78/0.78'
tree$node.label[28]='40/0.00/0.03'
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)
text(3, 3, "% supported/IC/TCA\nare shown on each node")

# Plot coalscence tree 
plot_tre("../tree/astral_tree.tre", name="ASTRAL Summary Coalescence Tree", 
         root="S.lycopersicum", rotate_node=c(18,22,26,27))
text(3, 3, "Bootstrap values\nsmaller than 100\nare shown on\ncorresponding node")

# Plot BUCKy tree
plot_tre("../tree/bucky_tree.tre", name="BUCKy Concordance Tree", 
         branch_length=F, rotate_node=c(19,23,24,27,28,29))
text(3, 3, "Concordance factor (CF)\nare shown on each\ninternal node")
dev.off()



### Parse the tree information used in ms simulation
mytree <- read.tree("../tree/RAxML_bipartitions.transcriptome")
ultra_tre <- ultrametric(mytree)
rotate_node=c(18,20,23,25,26,28,29)
for (x in 1:length(rotate_node)){
  if (!is.na(rotate_node[x])) ultra_tre <- rotate(ultra_tre, rotate_node[x])}
ultra_tre <- drop.tip(ultra_tre,"S.lycopersicum")
plot(ultra_tre)
tiplabels()
dist.nodes(ultra_tre)[4,5]/2



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


#pdf("Claudogram-trees-50.pdf", width=6, height=6) 
#par(mar=c(3,1,1,2)+0.1)
#trees50 <- read.tree("../tree/trees_bootstrap_50")
#ultrametric.trees50 <- list()
#ultrametric.trees50  <- lapply(trees50, ultrametric)
#class(ultrametric.trees50)<-"multiPhylo"
#densiTree(ultrametric.trees, type = "cladogram", alpha = 0.006, consensus = tree, 
#          optim = F, scaleX = T, col = "dark grey", width = 1, cex = 0.8)  ## boostrap_50
#dev.off()

