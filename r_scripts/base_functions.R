library(phytools)
library(phangorn)


# A function to generate tree
plot_tre <-function(raxml.tree, name=NA, tipdrop = NA, ultrametric=F, 
                    scale_axis=F, boostrap_pie=F, branch_length=T, 
                    outgroup=NA, rotate_node=NA){
  
  raxml.tree <- read.tree(raxml.tree)
  
  if (!is.na(outgroup)) raxml.tree <- root(raxml.tree, outgroup, resolve.root=T)
  if (ultrametric) raxml.tree<-ultrametric(raxml.tree)
  
  for (x in 1:length(rotate_node)){
    if (!is.na(rotate_node[x])) raxml.tree <- rotate(raxml.tree, rotate_node[x])}
  
  if (!is.na(tipdrop)){
    raxml.tree <- drop.tip(raxml.tree, tipdrop)
    raxml.tree$root.edge <- 0.001
    plot(raxml.tree, root.edge = TRUE, use.edge.length = branch_length, main=name)
    add.scale.bar(length=0.001)}
  else plot(raxml.tree, use.edge.length = branch_length, main=name)
  
  for (i in 1:length(raxml.tree$node.label)){
    if (is.null(raxml.tree$node.label[i])) raxml.tree$node.label[i]=''}
  for (i in 1:length(raxml.tree$node.label)){
    if (is.na(raxml.tree$node.label[i])) raxml.tree$node.label[i]=''}
  for (i in 1:length(raxml.tree$node.label)){
    if (raxml.tree$node.label[i]=='') raxml.tree$node.label[i]=100}
  for (i in 1:length(raxml.tree$node.label)){
    if (raxml.tree$node.label[i]=='Root') raxml.tree$node.label[i]=100}
  for (i in 1:length(raxml.tree$node.label)){
    if (as.numeric(raxml.tree$node.label[i])==100) raxml.tree$node.label[i]=''
    #else raxml.tree$node.label[i]='*'
  }
  
  if (boostrap_pie) {
    nodelabels(node=1:raxml.tree$Nnode+Ntip(raxml.tree),
               pie=cbind(as.numeric(raxml.tree$node.label),100-as.numeric(raxml.tree$node.label)),
               piecol=c("black","white"),cex=1)}
  else nodelabels(raxml.tree$node.label, adj = c(-0.2, 0.5), frame = "none", cex=0.6)
  
  if (scale_axis) {
    h<-max(nodeHeights(raxml.tree))
    axis(side = 1, at = seq(0, h, 0.5))}
}


# A function to generate ultrametric tree
ultrametric <- function(mytree){
  calibra <- makeChronosCalib(mytree, node = "root", age.min = 17,
                              age.max = 17, interactive = FALSE, soft.bounds = FALSE)
  chrono_tree <- chronos(mytree, lambda = 1, model = "discrete", quiet = FALSE,
                         calibration = calibra,
                         control = chronos.control())
  return(chrono_tree)
}
