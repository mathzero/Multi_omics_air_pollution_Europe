
### Denoising data ###

### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"


exposures <- readRDS("data_imputed/exposures.RDS")
covariates <- readRDS("data_imputed/covariates.RDS")
methylation <- readRDS("data_imputed/methylation_meanimpute.RDS")


####### DENOISING ######


library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(mvtnorm)
library(nlme)
library(mlbench)
library(doBy)
library(parallel)
library(foreach)
library(doParallel)
library(party)
library(omics)



exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) /2 

covExp <- inner_join(covariates,exposures, by = "subjectidp")





### Looping through each OMIC ###

denoised = NULL

data <- inner_join(covExp,methylation,by = "subjectidp")

#setup parallel backend to use many processors
n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)

t0=Sys.time()

X <- data[,27:ncol(data)]
fml <-  "X ~ (1|plate) + temp + relhum + gender + age + session + city + bmi + season + physical_act"
model <- mlmer(as.formula(fml), data=data[,1:ncol(data)], lrt=FALSE, save.residuals=TRUE, save.ranks=FALSE)

denoised <- resid(model)

denoised <- cbind(data[,1],denoised)
rownames(denoised) = rownames(data)
colnames(denoised) = colnames(methylation)


stopCluster(cl)
t1=Sys.time()
print(t1-t0)

denoised <- as.data.frame(denoised)

saveRDS(denoised, file = "methylation_denoised_noID.rds")
