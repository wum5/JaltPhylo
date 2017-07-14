setwd("/Users/mengwu/Documents/Research/JaltPhylo/updated_scripts/r_scripts")
library(ggplot2)


## Tree in Figure 2
pdf("fig2-tree.pdf", width=15, height=5) 
par(mar=c(3,1,1,2)+0.1)
plot_tre("../tree/bucky_tree.tre", ultrametric=T, scale_axis=T, boostrap_pie=T, rotate_node=c(19,23,24,27,28,29))
dev.off()



# Plot whole-transcriptome concatenated tree
pdf("Whole-transcriptome-tree.pdf", width=10, height=10) 
par(mar=c(0,0,1,0)+0.1)
par(mfrow=c(2,2))
plot_tre("../tree/RAxML_bipartitions.transcriptome", name="RAxML Concatenated Tree", rotate_node=c(18,20,21,23,29))
text(0.01, 3, "Bootstrap values\nsmaller than 100\nare shown in\ncorresponding node")

# Plot consensus tree 
tree <- read.tree("../tree/consensus_tree2.tre")
tree$edge.length<-NULL # get rid of edge lengths
tree <- rotate(tree, 27)
tree <- rotate(tree, 24)
plot(tree, main="Consensus Tree")
tree$node.label[3]='50/45/39'
tree$node.label[4]='22/26/21'
tree$node.label[5]='9/3/2'
tree$node.label[6]='6/0/0'
tree$node.label[7]='14/8/7'
tree$node.label[8]='16/0/0'
tree$node.label[9]='5/0/0'
tree$node.label[10]='7/0/0'
tree$node.label[11]='14/1/3'
tree$node.label[12]='46/63/63'
tree$node.label[13]='71/60/60'
tree$node.label[14]='34/0/1'
nodelabels(tree$node.label, adj = c(1.1, 1.5), frame = "none", cex=0.6)
text(3, 3, "% supported/IC/TCA\nare shown on each node")

# Plot coalscence tree 
plot_tre("../tree/astral_tree.tre", name="ASTRAL Summary Coalescence Tree", rotate_node=c(17,19,23,27,29))
text(3, 3, "Bootstrap values\nsmaller than 100\nare shown on\ncorresponding node")

# Plot BUCKy tree
plot_tre("../tree/bucky_tree.tre", name="BUCKy Concordance Tree", branch_length=F, rotate_node=c(19,23,24,27,28,29))
text(3, 3, "Concordance factor (CF)\nare shown on each\ninternal node")
dev.off()



# Plot claudogram tree
pdf("Claudogram-trees-50.pdf", width=6, height=6) 
par(mar=c(3,1,1,2)+0.1)
trees <- read.tree("../tree/bootstrap_50.tre")
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



## Plot each chromosome tree
pdf("Chromosomes-trees.pdf", width=10, height=10) 
par(mar=c(0,0,2,0)+0.1)
par(mfrow=c(4,3))
plot_tre("../tree/RAxML_bipartitions.ch1", name="Chromosome 1", rotate_node=c(17,19,23,24,28,29))
plot_tre("../tree/RAxML_bipartitions.ch2", name="Chromosome 2", rotate_node=c(23,24,27,29))
plot_tre("../tree/RAxML_bipartitions.ch3", name="Chromosome 3", rotate_node=c(17,19,25,26,27,28,29))
plot_tre("../tree/RAxML_bipartitions.ch4", name="Chromosome 4", rotate_node=c(27,29))
plot_tre("../tree/RAxML_bipartitions.ch5", name="Chromosome 5", rotate_node=c(17,18,19,22,23,25,27,28))
plot_tre("../tree/RAxML_bipartitions.ch6", name="Chromosome 6", rotate_node=c(17,19,20,21,23,24,25,26,29))
plot_tre("../tree/RAxML_bipartitions.ch7", name="Chromosome 7", rotate_node=c(20,23))
plot_tre("../tree/RAxML_bipartitions.ch8", name="Chromosome 8", rotate_node=c(17,19,22,23,25))
plot_tre("../tree/RAxML_bipartitions.ch9", name="Chromosome 9", rotate_node=c(20,22,23,24,25,26,28,29))
plot_tre("../tree/RAxML_bipartitions.ch10", name="Chromosome 10", rotate_node=c(17,18,19,20,22,23,26,28))
plot_tre("../tree/RAxML_bipartitions.ch11", name="Chromosome 11", rotate_node=c(17,23,24,25,27,28))
plot_tre("../tree/RAxML_bipartitions.ch12", name="Chromosome 12", rotate_node=c(17,25,26,28))
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
pdf("branch_lengths_vs_concordance.pdf", width=5.2, height=5) 
mytree <- read.tree("../tree/bucky_tree.tre")
ultra_tre <- ultrametric(mytree)
rotate_node=c(19,23,24,27,28,29)
for (x in 1:length(rotate_node)){
  if (!is.na(rotate_node[x])) ultra_tre <- rotate(ultra_tre, rotate_node[x])}
#plot(ultra_tre)
#edgelabels(round(ultra_tre$edge.length, digits = 2),frame="none",bg=NA,adj = c(0.5, 0))
brchLen <- c(0.16,1.22,0.54,1.32,0.3,0.13,0.1,0.08,0.06,0.12,0.12,0.25)
concord <- c(41,87,67,66,37,19,14,11,11,20,20,27)
subclad <- c("black","black","green","grey","grey", rep("orange",7))
df <- data.frame(brchLen,concord,subclad)
cor.test(df$brchLen,df$concord) 
qplot(brchLen, concord, data=df, color=subclad, size=1, xlim = c(0,1.5), ylim = c(0,100), 
      xlab="Internal Branch Length", ylab="Internode Concordance (%)") + theme_bw() + 
      theme(legend.position="none", panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), axis.text=element_text(size=12), 
            axis.title=element_text(size=16)) +
      scale_color_manual(breaks = subclad, values= c("purple","dark green","grey","orange")) +
      annotate("text", x = 1.1, y=20, label = "Pearspn correlation:r=0.89\nP-value: 0.0001")
dev.off()



#### comparing two alternative topologies
pdf("alt_topology.pdf", width=10, height=5) 
par(mfrow=c(1,2))
tree1 <- read.tree(text = "((J. umbellata,J. aijana),(J. sinuosa,J. biflora));")
plot(tree1, edge.width = 2, cex=1.5)
tree2 <- read.tree(text = "(((J. umbellata,J. aijana),J. sinuosa),J. biflora);")
plot(tree2, edge.width = 2, cex=1.5)
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





