### Feature selection using random forest ###

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

###########################################################################
###########################################################################

# Reading in data #

exposures <- readRDS("data_denoised/exposures.rds")
proteins <- readRDS("data_denoised/proteins_denoised.rds")
metabolites <- readRDS("data_denoised/metabolites_denoised.rds")
covariates <- readRDS("data_denoised/covariates.rds")
transcripts <- readRDS("data_denoised/transcripts_denoised.rds")
#methylation <- readRDS("data_denoised/methylation_denoised.rds")



exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) /2 
exposures$pm25delta <-  exposures$pm25_comb - exposures$modeledpm25abs
exposures$pncdelta <-  exposures$pncmedian - exposures$modeledpnc

omics <- list(proteins,transcripts,metabolites,methylation)

omic.names <- c("proteins","transcripts","metabolites","methylation")

###########################################################################
###########################################################################

var.importance.pm25 <- list()
var.importance.pnc <- list()


t0=Sys.time()
#setup parallel backend to use many processors
n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)


###Looping for variable selection in PM25 ###
for (i in length(omics)){
  data <- inner_join(exposures[,c("pm25delta", "subjectidp")], omics[[i]], by = "subjectidp")
  data <- subset(data, select=-subjectidp)
  
  # Create model with default paramters
  control <- trainControl(method="repeatedcv", number=10, repeats=3,verboseIter = F
                          , savePredictions = T
                          , returnResamp = "all",
                          allowParallel = T)
  seed <- 7
  set.seed(seed)
  mtry <- round(sqrt(ncol(data)))
  lass_default <- train(pm25delta~., data=data, method="lasso", tuneLength = 4, trControl=control, importance=T)
  print(lass_default)
  png(
    file = paste0(omic.names[[i]], "_variImp_pm25LASS.png"),
    width     = 5,
    height    = 8,
    units     = "in",
    res       = 1200,
    pointsize = 4
  )
  plot(varImp(lass_default))
  dev.off()
  var.importance.pm25[[i]] <- varImp(lass_default)
}

stopCluster(cl)
print(t1-t0) # around 2 minutes


t0=Sys.time()
#setup parallel backend to use many processors
n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)


###Looping for variable selection in PM25 ###
for (i in length(omics)){
  data <- inner_join(exposures[,c("pncdelta", "subjectidp")], omics[[i]], by = "subjectidp")
  data <- subset(data, select=-subjectidp)
  # Create model with default paramters
  control <- trainControl(method="repeatedcv", number=10, repeats=3,verboseIter = F
                          , savePredictions = T
                          , returnResamp = "all",
                          allowParallel = T)
  seed <- 7
  set.seed(seed)
  mtry <- round(sqrt(ncol(data)))
  lass_default <- train(pncdelta~., data=data, method="rf", tuneLength = 4, trControl=control, importance=T)
  print(lass_default)
  png(
    file = paste0(omic.names[[i]], "_variImp_pncLASS.png"),
    width     = 5,
    height    = 8,
    units     = "in",
    res       = 1200,
    pointsize = 4
  )
  plot(varImp(lass_default))
  dev.off()
  var.importance.pm25[[i]] <- varImp(lass_default)

}


stopCluster(cl)
print(t1-t0) # around 2 minutes

names(var.importance.pm25) <- c("proteins","transcripts","metabolites","methylation")
names(var.importance.pnc) <- c("proteins","transcripts","metabolites","methylation")

save(var.importance.pm25,var.importance.pnc, file = "lass_variable_selection.rds")



