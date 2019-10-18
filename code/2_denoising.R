### Denoising data ###

### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"


exposures <- readRDS("data_imputed/exposures.rds")
proteins <- readRDS("data_imputed/proteins_meanimpute.rds")
metabolites <- readRDS("data_imputed/metabolites_meanimpute.rds")
covariates <- readRDS("data_imputed/covariates.rds")
transcripts <- readRDS("data_imputed/transcripts_meanimpute.rds")
methylation1 <- readRDS("data_imputed/methylation_meanimpute.rds")


####### DENOISING ######


library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)



exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) /2 

covExp <- inner_join(covariates,exposures, by = "subjectidp")
random.effects <- list("(1|plate)","(1|isolation) + (1|labeling) + (1|hybridization)",0, "(1|chip) + (1|position)")
omics <- list(proteins,transcripts,metabolites,methylation)

omic.names <- c("proteins","transcripts","metabolites","methylation")


### Looping through each OMIC ###


for (i in 1:3){
  
  tol = 1e-04
  denoised <-  NULL
  data <- inner_join(covExp,omics[[i]],by = "subjectidp")
  denoised <- data.frame(cbind(denoised, data$subjectidp))
  
  
  ### Ensuring that only complete cases are used for each OMIC (some tech confounders are missing) ###
  
  if (i == 1){
    data <- data[complete.cases(data[,14]),] ### plate
  }
  if (i == 2){
    data <- data[complete.cases(data[,15]),] ### isolation
  }
  if (i == 4){
    data <- data[complete.cases(data[,18]),] ### chip
  }

  
  ### Creating variable X - all the 'outcome' variables (ie the OMICs) ###
  
  
  
    t0=Sys.time()
    for (k in 27:ncol(data)){
      print(k)
    fml <-  paste("data[,k] ~ ", random.effects[[i]])
    model <- lmer(as.formula(fml), data=data, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    t1=Sys.time()
    print(t1-t0) # around 2 minutes
    
    denoised <- data.frame(cbind(denoised,resid(model)))
    }
  rownames(denoised) = rownames(data)
  colnames(denoised) = colnames(omics[[i]])
  
  saveRDS(denoised, file = paste0(omic.names[[i]],"_denoised_RE.rds"))
}







#####################################################################################
#####################################################################################

### Denoising data ###

### Your working directory on the server
work_dir<-"/rds/general/user/mw418/home/comp_epi"

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



exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) /2 

covExp <- inner_join(covariates,exposures, by = "subjectidp")





### Looping through each OMIC ###



tol = 1e-04
denoised = NULL

data <- inner_join(covExp,methylation,by = "subjectidp")

#setup parallel backend to use many processors
n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)

t0=Sys.time()
for (k in 27:52){
  print(k)
  fml <-  "data[,k] ~ (1|id.x) + (1|plate) + temp + relhum + gender + age + session + city + bmi + season + physical_act"
  model <- lmer(as.formula(fml), data=data, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  denoised <- data.frame(cbind(denoised,resid(model)))
}
t1=Sys.time()
print(t1-t0)

stopCluster(cl)
t1=Sys.time()
print(t1-t0)

rownames(denoised) = rownames(data)
colnames(denoised) = colnames(methylation)

saveRDS(denoised, file = "methylation_denoised.rds")


