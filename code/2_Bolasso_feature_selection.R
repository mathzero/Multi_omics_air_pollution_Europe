#### Bolasso Variable selection ###
rm(list = ls())

### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"


library(ggplot2)
library(dplyr)
library(caret)
library(randomForest)
library(mlbench)
library(doBy)
library(parallel)
library(foreach)
library(doParallel)
library(party)
library(mht)
library(plotly)

###########################################################################
###########################################################################

# Reading in data #

exposures <- readRDS("data_denoised/exposures.rds")
proteins <- readRDS("data_denoised/proteins_denoised.rds")
metabolites <- readRDS("data_denoised/metabolites_denoised.rds")
covariates <- readRDS("data_denoised/covariates.rds")
transcripts <- readRDS("data_denoised/transcripts_denoised.rds")
methylation <- readRDS("data_denoised/methylation_denoised.rds")



exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) /2 
exposures$pm25delta <-  exposures$pm25_comb - exposures$modeledpm25abs
exposures$pncdelta <-  exposures$pncmedian - exposures$modeledpnc

omics <- list(proteins,transcripts,metabolites,methylation)

omic.names <- c("proteins","transcripts","metabolites","methylation")

###########################################################################
###########################################################################

var.importance.pm25.boll <- list()
var.importance.pnc.boll <- list()

# 
# t0=Sys.time()
# #setup parallel backend to use many processors
# n_cores=detectCores()
# cl <- makeCluster(n_cores-1) #not to overload your computer
# registerDoParallel(cl)


###Looping for variable selection in PM25 ###
for (i in 1:length(omics)){
  i <- 1
  
  data <- inner_join(exposures[,c("pm25delta", "subjectidp")], omics[[i]], by = "subjectidp")
  data <- subset(data, select=-subjectidp)
  
  X <- as.matrix(data[,2:ncol(data)])
  y <- as.matrix(data[,"pm25delta"]) 
    
  # Create model with default paramters
  seed <- 7
  set.seed(seed)
  mod <- bolasso(X,y,mu=seq(0.001,0.5,0.02))
  png(
    file = paste0(omic.names[[i]], "_varImp_pm25_bolasso.png"),
    width     = 7,
    height    = 4,
    units     = "in",
    res       = 1200,
    pointsize = 4
  )
  plot(mod)
  dev.off()
  var.importance.pm25.boll[[i]] <- mod

  data <- inner_join(exposures[,c("pncdelta", "subjectidp")], omics[[i]], by = "subjectidp")
  data <- subset(data, select=-subjectidp)
  # Create model with default paramters
  
  X <- as.matrix(data[,2:ncol(data)])
  y <- as.matrix(data[,"pncdelta"]) 
  
  seed <- 7
  set.seed(seed)
  mod <- bolasso(X,y,mu=seq(0.001,0.5,0.02))
  png(
    file = paste0(omic.names[[i]], "_varImp_pnc_bolasso.png"),
    width     = 7,
    height    = 4,
    units     = "in",
    res       = 1200,
    pointsize = 4
  )
  plot(mod)
  dev.off()
  var.importance.pm25.boll[[i]] <- mod

}


Sys.setenv("plotly_username"="mathzero")
Sys.setenv("plotly_api_key"="zB4ajmpm27Zau1fuWnXx")

testmatrix <- as.data.frame(mod$frequency)
testmatrix <- testmatrix[2:(nrow(testmatrix)),18:ncol(testmatrix)]
proteins.heatmap <- plot_ly(z = as.matrix(testmatrix), x=colnames(testmatrix), y = rownames(testmatrix), type = "heatmap")

api_create(proteins.heatmap, filename = "proteins_heatmap")



# stopCluster(cl)
# t1=Sys.time()
# print(t1-t0) # around 2 minutes

names(var.importance.pm25.boll) <- c("proteins","transcripts","metabolites","methylation")
names(var.importance.pnc.boll) <- c("proteins","transcripts","metabolites","methylation")

save(var.importance.pm25.boll,var.importance.pnc.boll, file = "bolasso_variable_selection.rdata")

