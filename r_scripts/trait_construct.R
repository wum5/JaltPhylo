setwd("/Users/mengwu/Documents/Research/JaltPhylo/updated_scripts/r_scripts")
source("base_functions.R")
library(ape)


#### Ancestral state reconstruction
mytree <- read.tree("../tree/bucky_tree.tre")
ultra_tre <- ultrametric(mytree)
rotate_node=c(19,23,24,27,28,29)
for (x in 1:length(rotate_node)){
  if (!is.na(rotate_node[x])) ultra_tre <- rotate(ultra_tre, rotate_node[x])}

fruit_color = c(rep("black",3),"red",rep("green",2),rep("orange",8))
nectar_color = c(rep("clear",4),"red","red","clear","red","clear",rep("red",5))
nectar_volume = as.matrix(c(2,2,2,3,50,51,3,137,7,23,18,38,33,15))
floral_shape = c(rep("rotate",4),"companulate","companulate","rotate","tubular","rotate",
                 "tubular","tubular","companulate","companulate","tubular")
df = data.frame(fruit_color, nectar_color, nectar_volume, floral_shape)
row.names(df) <- c("J.repandidentata","J.procumbens","J.darcyana",
                   "J.auriculata","J.quipuscoae","J.calliantha","J.yungayensis",
                   "J.biflora","J.sinuosa","J.aijana","J.umbellata",
                   "J.grandibaccata","J.dendroidea","J.incahuasina")
new_tree<-drop.tip(ultra_tre,tip="S.lycopersicum")


#### Maximum likelihood Approach ####

pdf("ancestral_states.pdf", width=10, height=10) 
par(mfrow=c(2,2))

## estimate ancestral states under a ER model
## fruit color
fruit_color<-as.matrix(df)[,1]
cols<-setNames(c("purple","green 3","orange","red"),sort(unique(fruit_color)))

fitER<-ace(fruit_color,new_tree,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)
plotTree(new_tree,fsize=1,ftype="i")
nodelabels(node=1:new_tree$Nnode+Ntip(new_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.8)
tiplabels(pie=to.matrix(fruit_color[new_tree$tip.label],
                        seq=sort(unique(fruit_color))),piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.3,y=4,fsize=1)

# nectar color
nectar_color<-as.matrix(df)[,2]
cols<-setNames(c("white","red"),sort(unique(nectar_color)))

fitER<-ace(nectar_color,new_tree,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)
plotTree(new_tree,fsize=1,ftype="i")
nodelabels(node=1:new_tree$Nnode+Ntip(new_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.8)
tiplabels(pie=to.matrix(nectar_color[new_tree$tip.label],
                        seq=sort(unique(nectar_color))),piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.3,
                  y=3,fsize=1)

# nectar volume
nectar_volume<-as.matrix(df)[,3]
nectar_volume<-setNames(sqrt(c(2,2,2,3,50,51,3,137,7,23,18,38,33,15)),
                        c("J.repandidentata","J.procumbens","J.darcyana",
                          "J.auriculata","J.quipuscoae","J.calliantha","J.yungayensis",
                          "J.biflora","J.sinuosa","J.aijana","J.umbellata",
                          "J.grandibaccata","J.dendroidea","J.incahuasina"))
fit<-fastAnc(new_tree,nectar_volume,vars=TRUE,CI=TRUE)

obj<-contMap(new_tree,nectar_volume,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(new_tree)),
     fsize=c(1,0.9))

# floral shape
floral_shape<-as.matrix(df)[,4]
cols<-setNames(c("pink","white","red"),sort(unique(floral_shape)))

fitER<-ace(floral_shape,new_tree,model="ER",type="discrete")
fitER
round(fitER$lik.anc,3)
plotTree(new_tree,fsize=1,ftype="i")
nodelabels(node=1:new_tree$Nnode+Ntip(new_tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.8)
tiplabels(pie=to.matrix(floral_shape[new_tree$tip.label],
                        seq=sort(unique(floral_shape))),piecol=cols,cex=0.5)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.3,
                  y=3,fsize=1)
dev.off()


