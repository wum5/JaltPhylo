setwd("/Users/mengwu/Documents/Research/JaltPhylo/updated_scripts/r_scripts")
library(tidyverse)
library(pepo)
library(treeplyr)
library(ggtree, quietly=T)


## The function is find the linear correlation relationship between the branch lengths in
## concatenated tree and the coalescence tree generated from MP-EST. The topology of the 
## concatenated and coalescence tree must be exactly the same.
brchL_correlation <- function(mpest_tree, concate_tree){
  x = c()
  y = c()
  for (i in 1:length(mpest_tree$edge.length)) {
    if (mpest_tree$edge.length[i] != 9 & mpest_tree$edge.length[i] != 7) {
      x = c(x,concate_tree$edge.length[i])
      y = c(y,mpest_tree$edge.length[i])}}
  fit <- lm(y~x)
  intercept <- coef(fit)[1] 
  slope <- coef(fit)[2] 
  return(c(slope,intercept))
}


## The function is to transform the branch length '9.0/7.0' to the correlated branch based
## on the linear correlation.
brchL_trans <- function(slope,intercept,mpest_tree,concate_tree){
  for (i in 1:length(mpest_tree$edge.length)) {
    if (mpest_tree$edge.length[i] == 9 | mpest_tree$edge.length[i] == 7) 
      mpest_tree$edge.length[i] <- slope*concate_tree$edge.length[i] + intercept
  }
  return(mpest_tree)
}



# Load the data
mpest_tree <- read.tree("../tree/mpest_tree.tre")
mpest_tree <- rotateNodes(mpest_tree,c(16,18,19,20,21,27))
plot(mpest_tree)

concat_tree <- read.tree("../tree/RAxML_bipartitions.transcriptome")
concat_tree <- rotateNodes(concat_tree,c(18,20,21,23,29))
plot(concat_tree)

# find out the correaltion of branch length in two different trees
my_corr <- brchL_correlation(mpest_tree,concat_tree)

# correct the long branch in MP-EST tree
mpest_tree <- brchL_trans(my_corr[1],my_corr[2],mpest_tree,concat_tree)
mpest_tree$edge.length
plot(mpest_tree)

# start the P(e)/P(o) test
mpest_tree_branches <- prep_branch_lengths(mpest_tree) 
mpest_tree_branches

mpest_tree_hrf <- tree_hrf(mpest_tree_branches)
mpest_tree_hrf

solgg <- to_treedata(mpest_tree, mpest_tree_hrf)

ggtree(solgg, aes(color=hrf), size=2)+geom_tiplab(color='black')+
  scale_color_gradient2(low = "darkblue", mid = "yellow", high="red", midpoint=1, na.value="grey80")+
  theme(legend.position = c(.05, .85))





