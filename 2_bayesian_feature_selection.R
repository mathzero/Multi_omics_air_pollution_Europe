#### Bayesian Variable selection ###
rm(list = ls())

### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"

library(varbvs)
library(dplyr)
###########################################################################
###########################################################################

# Reading in data #

exposures <- readRDS("data_denoised/exposures.rds")
proteins <- readRDS("data_denoised/proteins_denoised.rds")
metabolites <- readRDS("data_denoised/metabolites_denoised.rds")
covariates <- readRDS("data_denoised/covariates.rds")
transcripts <- readRDS("data_denoised/transcripts_denoised.rds")
methylation <- as.data.frame(readRDS("data_denoised/methylation_denoised.rds"))

# colnames(methylation)[1] <- "subjectidp"
# 
# saveRDS(methylation,"data_denoised/methylation_denoised.rds")

exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) / 2 
exposures$pm25delta <-  exposures$pm25_comb - exposures$modeledpm25abs
exposures$pncdelta <-  exposures$pncmedian - exposures$modeledpnc

omics <- list(proteins,transcripts,metabolites,methylation)
omic.names <- c("proteins","transcripts","metabolites","methylation")

load("results/bayesian_variable_selection/bayesian_variable_selection.rdata")
###########################################################################
###########################################################################
# for (i in 4){
# 
# data <- inner_join(exposures[,c("pm25delta", "subjectidp")], omics[[4]], by = "subjectidp")
# data <- subset(data, select=-subjectidp)
# 
# X <- as.matrix(subset(data, select=-pm25delta))
# y <-  matrix(data$pm25delta)
# 
# bvsmodelpm25 <- varbvs(X, Z=NULL, y, family = c("gaussian"),  update.sigma = TRUE,
#              update.sa = TRUE, optimize.eta = TRUE, initialize.params = FALSE,
#              nr = 100, sa0 = 1, n0 = 10, tol = 1e-4, maxiter = 1e4,
#              verbose = TRUE)
# 
# bvs.mods.pm25[[4]] <- bvsmodelpm25
# print(summary(bvsmodelpm25))
# plot(bvs.mods.pm25[[4]])

data <- inner_join(exposures[,c("pncdelta", "subjectidp")], omics[[4]], by = "subjectidp")
data <- subset(data, select=-subjectidp)

X <- as.matrix(subset(data, select=-pncdelta))
y <-  matrix(data$pncdelta)

bvsmodelpnc <- varbvs(X, Z=NULL, y, family = c("gaussian"),  update.sigma = TRUE,
                       update.sa = TRUE, optimize.eta = TRUE, initialize.params = FALSE,
                       nr = 300, sa0 = 1, n0 = 10, tol = 1e-4, maxiter = 1e4,
                       verbose = TRUE)

bvs.mods.pnc[[4]] <- bvsmodelpnc
print(summary(bvsmodelpnc))
#plot(bvs.mods.pnc[[4]])
# }

# names(bvs.mods.pm25) <- c("proteins","transcripts","metabolites", "methylation")
#names(bvs.mods.pnc) <- c("proteins","transcripts","metabolites", "methylation")

saveRDS(bvs.mods.pm25,bvs.mods.pnc, file = "bayesian_variable_selection.rdata")

# for (i in 1:4){
#   png(filename = paste0(omic.names[[4]],"_bvs_pm25.png"),
#       width = 1000, height = 1000)
#   plot(bvs.mods.pm25[[4]], 
#        draw.threshold = NA, gap = 0,col = "midnightblue", pch = 20, alpha = 0.5,
#        scales = NULL, xlab = "variable", ylab = "posterior probability",
#        main = paste0(omic.names[[4]]," varbvs model: PM25"),
#        abline.args = list(lty = "dotted",col = "orangered"),
#        vars.xyplot.args = list(pch = 20,col = "magenta"),
#        vars.ltext.args = list(col = "black",pos = 4,cex = 0.5),
#        par.settings = list(par.main.text = list(font = 1,cex = 0.8),
#                            layout.heights = list(axis.top = 0, axis.bottom = 0)))
#   dev.off()
#   }




