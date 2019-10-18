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

###########################################################################
###########################################################################
### Bayesian analysis #####################################################
###########################################################################
###########################################################################

results.bayes <- load("results/bayesian_variable_selection/bayesian_variable_selection.rdata")

### PM25 ###

bvs.top.markers.pm25 <- list()
omic.names <- c("proteins","transcripts","metabolites", "methylation")


for (i in 1:4){
  bvs.pm25.inclusion.probs <- as.data.frame(bvs.mods.pm25[[i]]$pip)
  bvs.pm25.inclusion.probs$marker <- rownames(bvs.pm25.inclusion.probs)
  colnames(bvs.pm25.inclusion.probs)[1] <- "probs"
  
  bvs.pm25.inclusion.probs <- bvs.pm25.inclusion.probs %>% arrange(-probs,)  %>% as.data.frame()
  bvs.pm25.inclusion.probs <- bvs.pm25.inclusion.probs[,c(2,1)]
  bvs.top.markers.pm25[[i]] <- bvs.pm25.inclusion.probs[1:20,]
  names(bvs.top.markers.pm25)[[i]] <- omic.names[[i]]
  
 print(ggplot(data = bvs.top.markers.pm25[[i]], aes(x = reorder(marker, -probs),probs)) + geom_col() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    coord_flip() + ggtitle(paste("Bayesian variable selection (PM25):",omic.names[[i]]))) + ylab("Posterior probability") + xlab("Biomarker")
 
  ggsave(filename = paste0(omic.names[[i]],"_bvs_var_select_PM25.png"),units = "cm", device = "png", dpi=300, width = 20, height = 20)
}



### PNC ###

bvs.top.markers.pnc <- list()


for (i in 1:4){
  bvs.pm25.inclusion.probs <- as.data.frame(bvs.mods.pnc[[i]]$pip)
  bvs.pm25.inclusion.probs$marker <- rownames(bvs.pm25.inclusion.probs)
  colnames(bvs.pm25.inclusion.probs)[1] <- "probs"
  
  bvs.pm25.inclusion.probs <- bvs.pm25.inclusion.probs %>% arrange(-probs,)  %>% as.data.frame()
  bvs.pm25.inclusion.probs <- bvs.pm25.inclusion.probs[,c(2,1)]
  bvs.top.markers.pnc[[i]] <- bvs.pm25.inclusion.probs[1:20,]
  names(bvs.top.markers.pnc)[[i]] <- omic.names[[i]]
  
  print(ggplot(data = bvs.top.markers.pnc[[i]], aes(x = reorder(marker, -probs),probs)) + geom_col() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          coord_flip() + ggtitle(paste("Bayesian variable selection (PNC):",omic.names[[i]]))) + ylab("Posterior probability") + xlab("Biomarker")
  
  ggsave(filename = paste0(omic.names[[i]],"_bvs_var_select_PNC.png"),units = "cm", device = "png", dpi=300, width = 20, height = 20)
}

###########################################################################
###########################################################################
### Bolasso analysis ######################################################
###########################################################################
###########################################################################

library(pheatmap)
library(RColorBrewer)
library(viridis)
load("results/bolasso_variable_selection/bolasso_variable_selection.rdata")

### PM25 ###

bol.top.markers.pm25 <- list()


for (i in 1:4){
  bol.pm25.inclusion.probs <- as.data.frame(var.importance.pm25.boll[[i]]$frequency)
  bol.pm25.inclusion.probs$marker <- rownames(bol.pm25.inclusion.probs)
  n <- ncol(bol.pm25.inclusion.probs)
  bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[2:nrow(bol.pm25.inclusion.probs),c(n,1:(n-1))]
  bol.pm25.inclusion.probs$average <- rowMeans(bol.pm25.inclusion.probs[,(ncol(bol.pm25.inclusion.probs)-2):ncol(bol.pm25.inclusion.probs)])
  bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs %>% arrange(-average)  %>% as.data.frame()
  bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[,1:(ncol(bol.pm25.inclusion.probs)-1)]
  # if (i == 1){
  # bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[1:13,]
  # } else {
  #   bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[1:20,]
  # }
  rownames(bol.pm25.inclusion.probs) <- bol.pm25.inclusion.probs$marker
  bol.top.markers.pm25[[i]] <- bol.pm25.inclusion.probs[,]
  names(bol.top.markers.pm25)[[i]] <- omic.names[[i]]
  my_heatmap <- pheatmap(
    mat               = as.matrix(bol.pm25.inclusion.probs[,17:ncol(bol.pm25.inclusion.probs)]),
    color             = inferno(10),
    border_color      = NA,
    drop_levels       = TRUE,
    fontsize          = 3,
    cluster_cols      = FALSE,
    cluster_rows      = FALSE,
    main              = paste("Bolasso variable selection (PM25):",omic.names[[i]]))
  save_pheatmap_png <- function(x, filename, width=2000, height=1800, res = 300) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  save_pheatmap_png(my_heatmap, paste0(omic.names[[i]],"_bol_var_select_pm25_full.png"))
  #ggsave(filename = paste0(omic.names[[i]],"_bol_var_select.png"),units = "cm", device = "png", dpi=300, width = 20, height = 20)
}


### PNC ###

bol.top.markers.pnc <- list()

for (i in 1:4){
  bol.pm25.inclusion.probs <- as.data.frame(var.importance.pnc.boll[[i]]$frequency)
  bol.pm25.inclusion.probs$marker <- rownames(bol.pm25.inclusion.probs)
  n <- ncol(bol.pm25.inclusion.probs)
  bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[2:nrow(bol.pm25.inclusion.probs),c(n,1:(n-1))]
  bol.pm25.inclusion.probs$average <- rowMeans(bol.pm25.inclusion.probs[,(ncol(bol.pm25.inclusion.probs)-2):ncol(bol.pm25.inclusion.probs)])
  bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs %>% arrange(-average)  %>% as.data.frame()
  bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[,1:(ncol(bol.pm25.inclusion.probs)-1)]
  if (i == 1){
    bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[1:13,]
  } else {
    bol.pm25.inclusion.probs <- bol.pm25.inclusion.probs[1:20,]
  }
  rownames(bol.pm25.inclusion.probs) <- bol.pm25.inclusion.probs$marker
  bol.top.markers.pnc[[i]] <- bol.pm25.inclusion.probs[1:20,]
  names(bol.top.markers.pnc)[[i]] <- omic.names[[i]]
  my_heatmap <- pheatmap(
    mat               = as.matrix(bol.pm25.inclusion.probs[,17:ncol(bol.pm25.inclusion.probs)]),
    color             = inferno(10),
    border_color      = NA,
    drop_levels       = TRUE,
    fontsize          = 7,
    cluster_cols      = FALSE,
    cluster_rows      = FALSE,
    main              = paste("Bolasso variable selection (PNC):",omic.names[[i]]))
  save_pheatmap_png <- function(x, filename, width=2000, height=1800, res = 300) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  save_pheatmap_png(my_heatmap, paste0(omic.names[[i]],"_bol_var_select_pnc.png"))
  #ggsave(filename = paste0(omic.names[[i]],"_bol_var_select.png"),units = "cm", device = "png", dpi=300, width = 20, height = 20)
}


###########################################################################
###########################################################################
### Univariate analysis ###################################################
###########################################################################
###########################################################################

load("results/univariate_analysis/univariate_results1-3.RData")

TRAP.names <- c("pm25.o", "pm25.p", "pm25.op", "pnc")

names(proteins.results) <- TRAP.names
names(transcripts.results) <- TRAP.names
names(metabolites.results) <- TRAP.names

for (i in 1:4){
  proteins.results[[i]]$OMIC.ID <- rownames(proteins.results[[i]])
  proteins.results[[i]]$bh_pval <-  p.adjust(proteins.results[[i]]$pval, method = "BH") 
  proteins.results[[i]]$bonf_pval <-  p.adjust(proteins.results[[i]]$pval, method = "bonferroni") 
  proteins.results[[i]] <- proteins.results[[i]] %>% arrange(pval)  %>% as.data.frame()
  transcripts.results[[i]]$OMIC.ID <- rownames(transcripts.results[[i]])
  transcripts.results[[i]]$bh_pval <-  p.adjust(transcripts.results[[i]]$pval, method = "BH") 
  transcripts.results[[i]]$bonf_pval <-  p.adjust(transcripts.results[[i]]$pval, method = "bonferroni") 
  transcripts.results[[i]] <- transcripts.results[[i]] %>% arrange(pval)  %>% as.data.frame()
  metabolites.results[[i]]$OMIC.ID <- rownames(metabolites.results[[i]])
  metabolites.results[[i]]$bh_pval <-  p.adjust(metabolites.results[[i]]$pval, method = "BH") 
  metabolites.results[[i]]$bonf_pval <-  p.adjust(metabolites.results[[i]]$pval, method = "bonferroni") 
  metabolites.results[[i]] <- metabolites.results[[i]] %>% arrange(pval)  %>% as.data.frame()
}

proteins.results[[1]]
plots.list <- list()

for (i in 1:length(proteins.results)){
  plots.list[[i]] <- ggplot(data=proteins.results[[i]], aes(x=coef, y=-log10(pval), col = -log10(bh_pval) > -log10(0.05))) +
    geom_point(alpha=0.7, size=3) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(bh_pval) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.9) +
    geom_hline(yintercept=-log10(0.05 / nrow(proteins.results[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot: univ. model of proteins against", TRAP.names[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(proteins.results[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.9) 
}

png("protein_volcanos.png",width = 8000, height = 2500, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[4]], ncol=3)
dev.off()



plots.list <- list()

for (i in 1:length(transcripts.results)){
  plots.list[[i]] <- ggplot(data=transcripts.results[[i]], aes(x=coef, y=-log10(pval), col = -log10(bh_pval) > -log10(0.05))) +
    geom_point(alpha=0.7, size=2) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(bh_pval) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.9) +
    geom_hline(yintercept=-log10(0.05 / nrow(transcripts.results[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot: univ. model of transcripts against", TRAP.names[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(transcripts.results[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.9) 
}

png("transcripts_volcanos.png",width = 8000, height = 2500, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[4]], ncol=3)
dev.off()



plots.list <- list()

for (i in 1:length(metabolites.results)){
  plots.list[[i]] <- ggplot(data=metabolites.results[[i]], aes(x=coef, y=-log10(pval), col = -log10(bh_pval) > -log10(0.05))) +
    geom_point(alpha=0.7, size=2) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(bh_pval) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.9) +
    geom_hline(yintercept=-log10(0.05 / nrow(metabolites.results[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot: univ. of metabolites against", TRAP.names[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(metabolites.results[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.9) 
}

png("metabolites_volcanos.png",width = 8000, height = 2500, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[4]], ncol=3)
dev.off()


###########################################################################
###########################################################################

results.univ.prot <- left_join(left_join(proteins.results[[1]],proteins.results[[2]], by = "OMIC.ID"),proteins.results[[4]],by = "OMIC.ID")
results.univ.prot <- results.univ.prot[,c(3,1,2,4,6,7,8,10,11,12)]
colnames(results.univ.prot) <- c("OMIC.ID","beta_pm25o","pval_pm25o", "pval_pm25o_BH",
                                  "beta_pm25p","pval_pm25p", "pval_pm25p_BH",
                                  "beta_pnc","pval_pnc", "pval_pnc_BH")
results.univ.prot$OMIC.level <- "protein"

prot.sig <- filter(results.univ.prot, pval_pm25o_BH <= 0.05 |pval_pm25p_BH <= 0.05 |pval_pnc_BH <= 0.05)

###########################################################################

results.univ.trans <- left_join(left_join(transcripts.results[[1]],transcripts.results[[2]], by = "OMIC.ID"),transcripts.results[[4]],by = "OMIC.ID")
results.univ.trans <- results.univ.trans[,c(3,1,2,4,6,7,8,10,11,12)]
colnames(results.univ.trans) <- c("OMIC.ID","beta_pm25o","pval_pm25o", "pval_pm25o_BH",
                                  "beta_pm25p","pval_pm25p", "pval_pm25p_BH",
                                  "beta_pnc","pval_pnc", "pval_pnc_BH")
results.univ.trans$OMIC.level <- "transcript"

trans.sig <-filter(results.univ.trans, pval_pm25o_BH <= 0.05 |pval_pm25p_BH <= 0.05 |pval_pnc_BH <= 0.05)


###########################################################################

results.univ.metab <- left_join(left_join(metabolites.results[[1]],metabolites.results[[2]], by = "OMIC.ID"),metabolites.results[[4]],by = "OMIC.ID")
results.univ.metab <- results.univ.metab[,c(3,1,2,4,6,7,8,10,11,12)]
colnames(results.univ.metab) <- c("OMIC.ID","beta_pm25o","pval_pm25o", "pval_pm25o_BH",
                                  "beta_pm25p","pval_pm25p", "pval_pm25p_BH",
                                  "beta_pnc","pval_pnc", "pval_pnc_BH")
results.univ.metab$OMIC.level <- "metabolite"

metab.sig <-filter(results.univ.metab, pval_pm25o_BH <= 0.05 |pval_pm25p_BH <= 0.05 |pval_pnc_BH <= 0.05)

###########################################################################

results.univ.sig <- rbind(prot.sig,trans.sig,metab.sig)
results.univ.sig$sigcount <- ((results.univ.sig$pval_pm25o_BH <= 0.05) +  
                                (results.univ.sig$pval_pm25p_BH <= 0.05) + 
                                (results.univ.sig$pval_pnc_BH <= 0.05))

write_csv(results.univ.sig, "results_univariate_significant.csv")


###########################################################################
###########################################################################

bvs.top.markers.pm25[[1]]$OMIC.level <- "protein"
bvs.top.markers.pm25[[2]]$OMIC.level <- "transcript"
bvs.top.markers.pm25[[3]]$OMIC.level <- "metabolite"
bvs.top.markers.pm25[[4]]$OMIC.level <- "methylate"
bvs.top.markers.pnc[[1]]$OMIC.level <- "protein"
bvs.top.markers.pnc[[2]]$OMIC.level <- "transcript"
bvs.top.markers.pnc[[3]]$OMIC.level <- "metabolite"
bvs.top.markers.pnc[[4]]$OMIC.level <- "methylate"

bvs.sigs.pm25 <- rbind(bvs.top.markers.pm25[[1]],bvs.top.markers.pm25[[2]],bvs.top.markers.pm25[[3]],bvs.top.markers.pm25[[4]]) 
bvs.sigs.pnc <- rbind(bvs.top.markers.pnc[[1]],bvs.top.markers.pnc[[2]],bvs.top.markers.pnc[[3]],bvs.top.markers.pnc[[4]]) 
bvs.sigs <- full_join(bvs.sigs.pm25,bvs.sigs.pnc, by = "marker")

bvs.sigs[is.na(bvs.sigs$OMIC.level.y), 5] <- bvs.sigs[is.na(bvs.sigs$OMIC.level.y), 3]
bvs.sigs <- bvs.sigs[,c(1,2,4,5)]

colnames(bvs.sigs) <- c("OMIC.ID","pm25.posterior","pnc.posterior","OMIC.level")
bvs.sigs <- bvs.sigs[!is.na(bvs.sigs$OMIC.ID), ]


###########################################################################
###########################################################################

bol.top.markers.pm25[[1]]$OMIC.level <- "protein"
bol.top.markers.pm25[[2]]$OMIC.level <- "transcript"
bol.top.markers.pm25[[3]]$OMIC.level <- "metabolite"
bol.top.markers.pm25[[4]]$OMIC.level <- "methylate"
bol.top.markers.pnc[[1]]$OMIC.level <- "protein"
bol.top.markers.pnc[[2]]$OMIC.level <- "transcript"
bol.top.markers.pnc[[3]]$OMIC.level <- "metabolite"
bol.top.markers.pnc[[4]]$OMIC.level <- "methylate"


bol.sigs.pm25 <- rbind(bol.top.markers.pm25[[1]],bol.top.markers.pm25[[2]],bol.top.markers.pm25[[3]],bol.top.markers.pm25[[4]]) 
bol.sigs.pnc <- rbind(bol.top.markers.pnc[[1]],bol.top.markers.pnc[[2]],bol.top.markers.pnc[[3]],bol.top.markers.pnc[[4]]) 
bol.sigs <- full_join(bol.sigs.pm25,bol.sigs.pnc, by = "marker")
bol.sigs <- bol.sigs[!is.na(bol.sigs$marker), ]

bol.sigs$probpm25 <- apply(bol.sigs[, 2:26], 1, max)
bol.sigs$probpnc <- apply(bol.sigs[, 28:52], 1, max)
bol.sigs[is.na(bol.sigs$OMIC.level.y), 53] <- bol.sigs[is.na(bol.sigs$OMIC.level.y), 27]


bol.sigs <- bol.sigs[,c(1,54,55,53)]
colnames(bol.sigs) <- c("OMIC.ID","BOL_probpm25", "BOL_probpnc","OMIC.level")


###########################################################################
###########################################################################

#Creating over all significant results table ##

results.sig <- full_join(results.univ.sig, bvs.sigs, by = "OMIC.ID")
results.sig <- full_join(results.sig, bol.sigs, by = "OMIC.ID")


results.sig[is.na(results.sig$OMIC.level), 18] <- results.sig[is.na(results.sig$OMIC.level), 11]
results.sig[is.na(results.sig$OMIC.level), 18] <- results.sig[is.na(results.sig$OMIC.level), 15]
results.sig <- results.sig[,c(1:10,12:14,16:ncol(results.sig))]


###########################################################################
###########################################################################

# creating 'count' columns for significance  in each variable selection method

results.sig$uni <- ifelse(is.na(results.sig$sigcount), 0,1)

results.sig$bvs <- ifelse(!is.na(results.sig$pm25.posterior) | !is.na(results.sig$pnc.posterior),1,0)

results.sig$bol <- ifelse(!is.na(results.sig$BOL_probpm25) | !is.na(results.sig$BOL_probpnc),1,0)


write_csv(results.sig[,c(16:19)], "sig_results_for_venn.csv")
results.sig[c(66,67,138:156),]$OMIC.level <- "metabolite"
write_csv(results.sig, "sig_results_all.csv")


###########################################################################
###########################################################################
univ.methy.results.pm25 <- readRDS("results/univariate_analysis/methylation_univariate_results_df_1_pm25.rds")
univ.methy.results.pnc <- readRDS("results/univariate_analysis/methylation_univariate_results_df_1_pnc.rds")


univ.methy.results.pm25$OMIC.ID <- rownames(univ.methy.results.pm25)
univ.methy.results.pnc$OMIC.ID <- rownames(univ.methy.results.pnc)

univ.methy.results.pm25 <- univ.methy.results.pm25[,c(9,2,1)]
univ.methy.results.pnc <- univ.methy.results.pnc[,c(9,2,1)]

univ.methy.results.pm25$bh_pval <-  p.adjust(univ.methy.results.pm25$pval, method = "BH") 
univ.methy.results.pm25$bonf_pval <-  p.adjust(univ.methy.results.pm25$pval, method = "bonferroni") 
univ.methy.results.pnc$bh_pval <-  p.adjust(univ.methy.results.pnc$pval, method = "BH") 
univ.methy.results.pnc$bonf_pval <-  p.adjust(univ.methy.results.pnc$pval, method = "bonferroni") 
colnames(univ.methy.results.pm25) <- c("OMIC.ID","coef","pval","bh_pval","bonf_pval")

ggplot(data=univ.methy.results.pm25, aes(x=coef, y=-log10(pval), col = -log10(bh_pval) > -log10(0.05))) +
  geom_point(alpha=0.7, size=1) +
  theme(legend.position = "none") +
  xlab("Beta") + ylab("-log10 BH adjusted p-value") +
  geom_text_repel(aes(label=ifelse(-log10(bh_pval) > -log10(0.05),as.character(OMIC.ID),'')), 
                  hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.9) +
  geom_hline(yintercept=-log10(0.05 / nrow(univ.methy.results.pm25)), linetype="dashed", color = "black", alpha = 0.6) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
  ggtitle("Volcano plot: univ. model of cPgs against pm25") +
  geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(univ.methy.results.pm25))),label="Significant pval (Bonferroni corrected)",
            vjust=-1, size = 2.5, col ="black")  +
  geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
            vjust=-1, size = 2.5, col ="black", alpha = 0.9) 


ggplot(data=univ.methy.results.pnc, aes(x=coef, y=-log10(pval), col = -log10(bh_pval) > -log10(0.05))) +
  geom_point(alpha=0.7, size=1) +
  theme(legend.position = "none") +
  xlab("Beta") + ylab("-log10 BH adjusted p-value") +
  geom_text_repel(aes(label=ifelse(-log10(bh_pval) > -log10(0.05),as.character(OMIC.ID),'')), 
                  hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.9) +
  geom_hline(yintercept=-log10(0.05 / nrow(univ.methy.results.pnc)), linetype="dashed", color = "black", alpha = 0.6) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
  ggtitle("Volcano plot: univ. model of cPgs against PNC") +
  geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(univ.methy.results.pnc))),label="Significant pval (Bonferroni corrected)",
            vjust=-1, size = 2.5, col ="black")  +
  geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
            vjust=-1, size = 2.5, col ="black", alpha = 0.9) 



###########################################################################
###########################################################################
### Stratified analysis #####################################################
###########################################################################
###########################################################################


load("univariate_results_stratified.RData")

for (i in 1:4){
  pm25_adj_p.results.stratified[[i]]$OMIC.ID <- rownames(pm25_adj_p.results.stratified[[i]])
  pm25_adj_p.results.stratified[[i]]$bh_pval <-  p.adjust(pm25_adj_p.results.stratified[[i]]$pval, method = "BH") 
  pm25_adj_p.results.stratified[[i]]$bonf_pval <-  p.adjust(pm25_adj_p.results.stratified[[i]]$pval, method = "bonferroni") 
  pm25_adj_p.results.stratified[[i]] <- pm25_adj_p.results.stratified[[i]] %>% arrange(pval)  %>% as.data.frame()
}

plots.list <- list()

for (i in 1:length(pm25_adj_p.results.stratified)){
  plots.list[[i]] <- ggplot(data=pm25_adj_p.results.stratified[[i]], aes(x=coef, y=-log10(pval), col = -log10(bh_pval) > -log10(0.05))) +
    geom_point(alpha=0.7, size=3) +
    theme(legend.position = "none") +
    xlab("Beta") + ylab("-log10 BH adjusted p-value") +
    geom_text_repel(aes(label=ifelse(-log10(bh_pval) > -log10(0.05),as.character(OMIC.ID),'')), 
                    hjust="inward", vjust="inward", col = 'black',size=1.5, alpha = 0.9) +
    geom_hline(yintercept=-log10(0.05 / nrow(pm25_adj_p.results.stratified[[i]])), linetype="dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "black", alpha = 0.3) +
    ggtitle(paste("Volcano plot: stratified univ. model of metabolites against PM25 in", locations[i])) +
    geom_text(mapping=aes(x=0,y=-log10(0.05 / nrow(pm25_adj_p.results.stratified[[i]]))),label="Significant pval (Bonferroni corrected)",
              vjust=-1, size = 2.5, col ="black")  +
    geom_text(mapping=aes(x=0,y=-log10(0.05)),label="Pval < 0.05",
              vjust=-1, size = 2.5, col ="black", alpha = 0.9) 
}

png("stratified_volcanos.png",width = 10000, height = 3000, res = 300)
grid.arrange(plots.list[[1]],plots.list[[2]],plots.list[[4]],plots.list[[4]], ncol=4)
dev.off()

pm25_adj_p.results.stratified[4]


univ.methy.results.pm25[1:10,]


