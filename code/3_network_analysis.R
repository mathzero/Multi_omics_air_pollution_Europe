### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"

library(varbvs)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)
library(pheatmap)
library(corpcor)
library(ppcor)
library(abind)
library(parallel)
library(RColorBrewer)
library(igraph)
library(psych)
library(dils)
library(ppcor)
library(gplots)

###########################################################################
results.sig <- read.csv("sig_results_all.csv")
results.sig[c(66,67,138:156),]$OMIC.level <- "metabolite"

exposures <- readRDS("data_denoised/exposures.rds")
proteins <- readRDS("data_denoised/proteins_denoised.rds")
metabolites <- readRDS("data_denoised/metabolites_denoised.rds")
covariates <- readRDS("data_denoised/covariates.rds")
transcripts <- readRDS("data_denoised/transcripts_denoised.rds")
methylation <- as.data.frame(readRDS("data_denoised/methylation_denoised.rds"))

###########################################################################
a <- exposures$subjectidp[(exposures$subjectidp  %in% proteins$subjectidp)]
a <- metabolites$subjectidp[(metabolites$subjectidp  %in% a)]
a <- covariates$subjectidp[(covariates$subjectidp  %in% a)]
a <- transcripts$subjectidp[(transcripts$subjectidp  %in% a)]
index <- methylation$subjectidp[(methylation$subjectidp  %in% a)]
index <- droplevels(index)
#creating lists of significant OMICs

sig.univ <- results.sig$OMIC.ID[results.sig$sigcount > 0 & !is.na(results.sig$sigcount) ]
sig.bvs <- results.sig$OMIC.ID[!is.na(results.sig$pm25.posterior) | !is.na(results.sig$pnc.posterior)]
sig.bol <- results.sig$OMIC.ID[!is.na(results.sig$BOL_probpm25) | !is.na(results.sig$BOL_probpnc)]

significant.matrix <- cbind(exposures[index,]$subjectidp,
                            proteins[index,(colnames(proteins) %in% sig.univ)],
                            metabolites[index,(colnames(metabolites) %in% sig.univ)],
                            transcripts[index,(colnames(transcripts) %in% sig.univ)],
                            methylation[index,(colnames(methylation) %in% sig.univ)])

colnames(significant.matrix)[1] <- "subjectidp"

omic.refs <- results.sig[,c(1,16)]
colnames(omic.refs)[1] <- "name"
omic.refs <- omic.refs[!duplicated(omic.refs$name), ]

# cols 2-28 METABOLITES
# cols 29-33 TRANSCRIPTS
###########################################################################

## Univariate networks

heatmap.2(cor(significant.matrix[,2:ncol(significant.matrix)]),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

cor.pws <- corr.test(significant.matrix[,2:ncol(significant.matrix)], use ="pairwise", method = "pearson", adjust = "none", alpha =0.05, ci=FALSE)
pvals <-  as.matrix(cor.pws$p)

heatmap(cor.pws$p)
                            
adj.matrix <- matrix(as.numeric(pvals < (0.05 / 32^2)),nrow =  ncol(pvals),ncol =  ncol(pvals))
rownames(adj.matrix) <- rownames(pvals)
colnames(adj.matrix) <- rownames(pvals)

for (i in 1:nrow(adj.matrix)){
  for (j in 1:ncol(adj.matrix)){
    if (i == j){
      adj.matrix[i,j] <- 0
    }
  }
}

network=graph_from_adjacency_matrix(adj.matrix)
# as_edgelist(network, names=T)
network <- simplify(network, remove.multiple = F, remove.loops = T)

###########################################################################
###########################################################################


png("Marginal_network_univariate_vars.png",width=5000,height=3000, res = 500) 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)

myplot <- delete.vertices(simplify(network), degree(network)==0)
ordered.vertices <-get.data.frame(myplot, what="vertices")
bound.frame <- left_join(ordered.vertices, omic.refs, by = "name")
coul = brewer.pal(nlevels(as.factor(bound.frame$OMIC.level)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(bound.frame$OMIC.level))]

# myplot$layout <- layout_nicely()

coords <- layout.circle(myplot)

plot(myplot, 
     vertex.size=6,
     vertex.color=my_color,
     vertex.label.cex=0.5,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.arrow.size=.003,
     vertex.label.degree = -pi/2,
     edge.width = pvals,
     layout = coords
)

# title and legend
text(-1,-1,"Multi-omic marginal network: univariate",col="white")
legend(x=-1.6, y=1, legend=levels(as.factor(bound.frame$OMIC.level)), col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.8, text.col="white" , horiz = F)
dev.off()

###########################################################################
###########################################################################
###########################################################################
###########################################################################


## BVS networks



significant.matrix <- cbind(exposures[index,]$subjectidp,
                            proteins[index,(colnames(proteins) %in% sig.bvs)],
                            metabolites[index,(colnames(metabolites) %in% sig.bvs)],
                            transcripts[index,(colnames(transcripts) %in% sig.bvs)],
                            methylation[index,(colnames(methylation) %in% sig.bvs)])

colnames(significant.matrix)[1] <- "subjectidp"
colnames(omic.refs)[1] <- "name"

# cols 2-28 METABOLITES
# cols 29-33 TRANSCRIPTS
###########################################################################

heatmap.2(cor(significant.matrix[,2:ncol(significant.matrix)]),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

cor.pws <- corr.test(significant.matrix[,2:ncol(significant.matrix)], use ="pairwise", method = "pearson", adjust = "none", alpha =0.05, ci=FALSE)
pvals <-  as.matrix(cor.pws$p)

heatmap(cor.pws$p)

adj.matrix <- matrix(as.numeric(pvals < (0.05 / ncol(pvals)^2)),nrow =  ncol(pvals),ncol =  ncol(pvals))
rownames(adj.matrix) <- rownames(pvals)
colnames(adj.matrix) <- rownames(pvals)

for (i in 1:nrow(adj.matrix)){
  for (j in 1:ncol(adj.matrix)){
    if (i == j){
      adj.matrix[i,j] <- 0
    }
  }
}

network=graph_from_adjacency_matrix(adj.matrix)
# as_edgelist(network, names=T)
network <- simplify(network, remove.multiple = F, remove.loops = T)

###########################################################################
###########################################################################


png("Marginal_network_bvs_vars.png",width=5000,height=3000, res = 500) 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)

myplot <- delete.vertices(simplify(network), degree(network)==0)
ordered.vertices <-get.data.frame(myplot, what="vertices")
bound.frame <- left_join(ordered.vertices, omic.refs, by = "name")
coul = brewer.pal(nlevels(as.factor(bound.frame$OMIC.level)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(bound.frame$OMIC.level))]

# myplot$layout <- layout_nicely()

coords <- layout_nicely(myplot)

plot(myplot, 
     vertex.size=3,
     vertex.color=my_color,
     vertex.label.cex=0.3,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.arrow.size=.003,
     vertex.label.degree = -pi/2,
     edge.width = pvals,
     layout = coords
)

# title and legend
text(-1,-1,"Multi-omic marginal network: Bayesian variable selection",col="white")
legend(x=-1.6, y=1, legend=levels(as.factor(bound.frame$OMIC.level)), col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.8, text.col="white" , horiz = F)
dev.off()
###########################################################################
###########################################################################
###########################################################################
###########################################################################

## Bolasso networks




significant.matrix <- cbind(exposures[index,]$subjectidp,
                            proteins[index,(colnames(proteins) %in% sig.bol)],
                            metabolites[index,(colnames(metabolites) %in% sig.bol)],
                            transcripts[index,(colnames(transcripts) %in% sig.bol)],
                            methylation[index,(colnames(methylation) %in% sig.bol)])

colnames(significant.matrix)[1] <- "subjectidp"
colnames(omic.refs)[1] <- "name"

# cols 2-28 METABOLITES
# cols 29-33 TRANSCRIPTS
###########################################################################

heatmap.2(cor(significant.matrix[,2:ncol(significant.matrix)]),dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

cor.pws <- corr.test(significant.matrix[,2:ncol(significant.matrix)], use ="pairwise", method = "pearson", adjust = "none", alpha =0.05, ci=FALSE)
pvals <-  as.matrix(cor.pws$p)

heatmap(cor.pws$p)

adj.matrix <- matrix(as.numeric(pvals < (0.05 / ncol(pvals)^2)),nrow =  ncol(pvals),ncol =  ncol(pvals))
rownames(adj.matrix) <- rownames(pvals)
colnames(adj.matrix) <- rownames(pvals)

for (i in 1:nrow(adj.matrix)){
  for (j in 1:ncol(adj.matrix)){
    if (i == j){
      adj.matrix[i,j] <- 0
    }
  }
}

network=graph_from_adjacency_matrix(adj.matrix)
# as_edgelist(network, names=T)
network <- simplify(network, remove.multiple = F, remove.loops = T)

###########################################################################
###########################################################################


png("Marginal_network_bol_vars.png",width=5000,height=3000, res = 500) 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)

myplot <- delete.vertices(simplify(network), degree(network)==0)
ordered.vertices <-get.data.frame(myplot, what="vertices")
bound.frame <- left_join(ordered.vertices, omic.refs, by = "name")
coul = brewer.pal(nlevels(as.factor(bound.frame$OMIC.level)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(bound.frame$OMIC.level))]

# myplot$layout <- layout_nicely()

coords <- layout_nicely(myplot)

plot(myplot, 
     vertex.size=3,
     vertex.color=my_color,
     vertex.label.cex=0.3,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.arrow.size=.003,
     vertex.label.degree = -pi/2,
     edge.width = pvals,
     layout = coords
)

# title and legend
text(-1,-1,"Multi-omic marginal network: Bolasso variable selection",col="white")
legend(x=-1.6, y=1, legend=levels(as.factor(bound.frame$OMIC.level)), col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.8, text.col="white" , horiz = F)
dev.off()




###########################################################################
###########################################################################
###########################################################################
###########################################################################
###########################################################################


## Conditional networks ###




##########################################################################
###########################################################################

a <- exposures$subjectidp[(exposures$subjectidp  %in% proteins$subjectidp)]
a <- metabolites$subjectidp[(metabolites$subjectidp  %in% a)]
a <- covariates$subjectidp[(covariates$subjectidp  %in% a)]
a <- transcripts$subjectidp[(transcripts$subjectidp  %in% a)]
index <- methylation$subjectidp[(methylation$subjectidp  %in% a)]
index <- droplevels(index)
#creating lists of significant OMICs

sig.univ <- results.sig$OMIC.ID[results.sig$sigcount > 0 & !is.na(results.sig$sigcount) ]
sig.bvs <- results.sig$OMIC.ID[!is.na(results.sig$pm25.posterior) | !is.na(results.sig$pnc.posterior)]
sig.bol <- results.sig$OMIC.ID[!is.na(results.sig$BOL_probpm25) | !is.na(results.sig$BOL_probpnc)]

significant.matrix <- cbind(exposures[index,]$subjectidp,
                            proteins[index,(colnames(proteins) %in% sig.univ)],
                            metabolites[index,(colnames(metabolites) %in% sig.univ)],
                            transcripts[index,(colnames(transcripts) %in% sig.univ)],
                            methylation[index,(colnames(methylation) %in% sig.univ)])

colnames(significant.matrix)[1] <- "subjectidp"

omic.refs <- results.sig[,c(1,16)]
colnames(omic.refs)[1] <- "name"
omic.refs <- omic.refs[!duplicated(omic.refs$name), ]

# cols 2-28 METABOLITES
# cols 29-33 TRANSCRIPTS
###########################################################################

## Univariate networks


cor.pws <- pcor.shrink(significant.matrix[,2:ncol(significant.matrix)])
adj.matrix <- matrix(as.numeric(abs(cor.pws) > 0.1),nrow =  ncol(cor.pws), ncol =  ncol(cor.pws))

rownames(adj.matrix) <- colnames(significant.matrix[,2:ncol(significant.matrix)])
colnames(adj.matrix) <- colnames(significant.matrix[,2:ncol(significant.matrix)])


for (i in 1:nrow(adj.matrix)){
  for (j in 1:ncol(adj.matrix)){
    if (i == j){
      adj.matrix[i,j] <- 0
    }
  }
}

network=graph_from_adjacency_matrix(adj.matrix)
# as_edgelist(network, names=T)
network <- simplify(network, remove.multiple = F, remove.loops = T)

###########################################################################
###########################################################################


png("conditional_network_univariate_vars.png",width=5000,height=3000, res = 500) 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)

myplot <- delete.vertices(simplify(network), degree(network)==0)
ordered.vertices <-get.data.frame(myplot, what="vertices")
bound.frame <- left_join(ordered.vertices, omic.refs, by = "name")
coul = brewer.pal(nlevels(as.factor(bound.frame$OMIC.level)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(bound.frame$OMIC.level))]

# myplot$layout <- layout_nicely()

coords <- layout_nicely(myplot)

plot(myplot, 
     vertex.size=6,
     vertex.color=my_color,
     vertex.label.cex=0.5,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.arrow.size=.003,
     vertex.label.degree = -pi/2,
     edge.width = cor.pws,
     layout = coords
)

# title and legend
text(-1,-1,"Multi-omic conditional network: univariate",col="white")
legend(x=-1.6, y=1, legend=levels(as.factor(bound.frame$OMIC.level)), col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.8, text.col="white" , horiz = F)
dev.off()

###########################################################################
###########################################################################
###########################################################################
###########################################################################


## BVS networks



significant.matrix <- cbind(exposures[index,]$subjectidp,
                            proteins[index,(colnames(proteins) %in% sig.bvs)],
                            metabolites[index,(colnames(metabolites) %in% sig.bvs)],
                            transcripts[index,(colnames(transcripts) %in% sig.bvs)],
                            methylation[index,(colnames(methylation) %in% sig.bvs)])

colnames(significant.matrix)[1] <- "subjectidp"
colnames(omic.refs)[1] <- "name"

# cols 2-28 METABOLITES
# cols 29-33 TRANSCRIPTS
###########################################################################

cor.pws <- pcor.shrink(significant.matrix[,2:ncol(significant.matrix)])
adj.matrix <- matrix(as.numeric(abs(cor.pws) > 0.1),nrow =  ncol(cor.pws), ncol =  ncol(cor.pws))
rownames(adj.matrix) <- colnames(significant.matrix[,2:ncol(significant.matrix)])
colnames(adj.matrix) <- colnames(significant.matrix[,2:ncol(significant.matrix)])



for (i in 1:nrow(adj.matrix)){
  for (j in 1:ncol(adj.matrix)){
    if (i == j){
      adj.matrix[i,j] <- 0
    }
  }
}

network=graph_from_adjacency_matrix(adj.matrix)
# as_edgelist(network, names=T)
network <- simplify(network, remove.multiple = F, remove.loops = T)

###########################################################################
###########################################################################


png("conditional_network_bvs_vars.png",width=5000,height=3000, res = 500) 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)

myplot <- delete.vertices(simplify(network), degree(network)==0)
ordered.vertices <-get.data.frame(myplot, what="vertices")
bound.frame <- left_join(ordered.vertices, omic.refs, by = "name")
coul = brewer.pal(nlevels(as.factor(bound.frame$OMIC.level)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(bound.frame$OMIC.level))]

# myplot$layout <- layout_nicely()

coords <- layout_nicely(myplot)

plot(myplot, 
     vertex.size=3,
     vertex.color=my_color,
     vertex.label.cex=0.3,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.arrow.size=.003,
     vertex.label.degree = -pi/2,
     edge.width = cor.pws,
     layout = coords
)

# title and legend
text(-1,-1,"Multi-omic conditional network: Bayesian variable selection",col="white")
legend(x=-1.6, y=1, legend=levels(as.factor(bound.frame$OMIC.level)), col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.8, text.col="white" , horiz = F)
dev.off()
###########################################################################
###########################################################################
###########################################################################
###########################################################################

## Bolasso networks




significant.matrix <- cbind(exposures[index,]$subjectidp,
                            proteins[index,(colnames(proteins) %in% sig.bol)],
                            metabolites[index,(colnames(metabolites) %in% sig.bol)],
                            transcripts[index,(colnames(transcripts) %in% sig.bol)],
                            methylation[index,(colnames(methylation) %in% sig.bol)])

colnames(significant.matrix)[1] <- "subjectidp"
colnames(omic.refs)[1] <- "name"

# cols 2-28 METABOLITES
# cols 29-33 TRANSCRIPTS
###########################################################################

cor.pws <- pcor.shrink(significant.matrix[,2:ncol(significant.matrix)])
adj.matrix <- matrix(as.numeric(abs(cor.pws) > 0.1),nrow =  ncol(cor.pws), ncol =  ncol(cor.pws))
rownames(adj.matrix) <- colnames(significant.matrix[,2:ncol(significant.matrix)])
colnames(adj.matrix) <- colnames(significant.matrix[,2:ncol(significant.matrix)])

for (i in 1:nrow(adj.matrix)){
  for (j in 1:ncol(adj.matrix)){
    if (i == j){
      adj.matrix[i,j] <- 0
    }
  }
}

network=graph_from_adjacency_matrix(adj.matrix)
# as_edgelist(network, names=T)
network <- simplify(network, remove.multiple = F, remove.loops = T)

###########################################################################
###########################################################################


png("conditional_network_bol_vars.png",width=5000,height=3000, res = 500) 
# plot
par(bg="grey13", mar=c(0,0,0,0))
set.seed(4)

myplot <- delete.vertices(simplify(network), degree(network)==0)
ordered.vertices <-get.data.frame(myplot, what="vertices")
bound.frame <- left_join(ordered.vertices, omic.refs, by = "name")
coul = brewer.pal(nlevels(as.factor(bound.frame$OMIC.level)), "Set2")
# Map the color to cylinders
my_color=coul[as.numeric(as.factor(bound.frame$OMIC.level))]

# myplot$layout <- layout_nicely()

coords <- layout_nicely(myplot)

plot(myplot, 
     vertex.size=3,
     vertex.color=my_color,
     vertex.label.cex=0.3,
     vertex.label.color="white",
     vertex.frame.color="transparent",
     edge.arrow.size=.003,
     vertex.label.degree = -pi/2,
     edge.width = pvals,
     layout = coords
)

# title and legend
text(-1,-1,"Multi-omic conditional network: Bolasso variable selection",col="white")
legend(x=-1.6, y=1, legend=levels(as.factor(bound.frame$OMIC.level)), col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 0.8, text.col="white" , horiz = F)
dev.off()


