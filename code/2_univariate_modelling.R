### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"


exposures <- readRDS("data_imputed/exposures.rds")
proteins <- readRDS("data_imputed/proteins.rds")
metabolites <- readRDS("data_imputed/metabolites_meanimpute.rds")
covariates <- readRDS("data_imputed/covariates_meanimpute.rds")
transcripts <- readRDS("data_imputed/transcripts_meanimpute.rds")
methylation <- readRDS("data_imputed/methylation_meanimpute.rds")



####### MODELLING ######

library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)


exposures$pm25_comb <- (exposures$pm25_adj_p+exposures$pm25_adj_o) /2 

covExp <- inner_join(covariates,exposures, by = "subjectidp")
random.effects <- list("(1|plate)","(1|isolation) + (1|labeling) + (1|hybridization)",0, "(1|chip) + (1|position)")

###########################################################################
###########################################################################

### Creating results lists ###

pm25_adj_p.results <- list()
pm25_adj_o.results <- list()
pm25_adj_combined.results <- list()
pncmedian.results <- list()


###########################################################################

#### Initiating loop ###

omics <- list(proteins,transcripts,metabolites,methylation)

for (i in 1:length(omics)){

  i <- 2
  
  data <- inner_join(covExp,omics[[i]],by = "subjectidp")
  if (i == 1){
    data <- data[complete.cases(data[,14]),] ### plate
  }
  if (i == 2){
      data <- data[complete.cases(data[,15]),] ### isolation
    }
  if (i == 4){
      data <- data[complete.cases(data[,18]),] ### chip
    }



  # omic.names <- c("proteins","transcripts","metabolites","methylation")
  
  # plate	Plate number (for proteins)
  # isolation	Technical covariate for transcripts
  # labeling	Technical covariate for transcripts
  # hybridization	Technical covariate for transcripts
  # chip	Chip number (for methylation data)
  # position	Position on the chip (for methylation data)



###########################################################################
###########################################################################

### Univariate modelling ###

### pm25_adj_p ###

univ.fun = function(X) {
  test.form <-  (paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpm25abs")))
  model0 = lmer(as.formula(test.form), 
                           data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                                                                                                                                            tol = 1e-04)))
  model1 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpm25abs + pm25_adj_p"))), 
                           data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                                                                                                                                                         tol = 1e-04)))
  vcov = as.data.frame(VarCorr(model1))$vcov
  res = c(summary(model1)$coefficients["pm25_adj_p",
                                       1:2], anova(model0, model1)$'Pr(>Chisq)'[2],
          vcov[1]/sum(vcov), vcov[2]/sum(vcov))
  names(res) = c("coef", "coef.se", "pval", "subjectid",
                 "plate")
  return(res)
}

results = t(apply(data[,27:ncol(data)], 2, FUN = univ.fun))

pm25_adj_p.results[[i]] <- data.frame(results)[,c(1,3)]


###########################################################################

### pm25_adj_o ###

univ.fun = function(X) {
  model0 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpm25abs"))), 
                data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                               tol = 1e-04)))
  model1 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpm25abs + pm25_adj_o"))), 
                data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                               tol = 1e-04)))
  vcov = as.data.frame(VarCorr(model1))$vcov
  res = c(summary(model1)$coefficients["pm25_adj_o",
                                       1:2], anova(model0, model1)$'Pr(>Chisq)'[2],
          vcov[1]/sum(vcov), vcov[2]/sum(vcov))
  names(res) = c("coef", "coef.se", "pval", "subjectid",
                 "plate")
  return(res)
}

results = t(apply(data[,27:ncol(data)], 2, FUN = univ.fun))

pm25_adj_o.results[[i]] <- data.frame(results)[,c(1,3)]



###########################################################################


### pm25_adj_o  AND  pm25_adj_p ###

univ.fun = function(X) {
  model0 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpm25abs"))), 
                data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                               tol = 1e-04)))
  model1 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpm25abs + pm25_comb"))), 
                data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                               tol = 1e-04)))
  vcov = as.data.frame(VarCorr(model1))$vcov
  res = c(summary(model1)$coefficients["pm25_comb",
                                       1:2], anova(model0, model1)$'Pr(>Chisq)'[2],
          vcov[1]/sum(vcov), vcov[2]/sum(vcov))
  names(res) = c("coef", "coef.se", "pval", "subjectid",
                 "plate")
  return(res)
}

results = t(apply(data[,27:ncol(data)], 2, FUN = univ.fun))

pm25_adj_combined.results[[i]] <- data.frame(results)[,c(1,3)]


###########################################################################



### PNC median ###


univ.fun = function(X) {
  model0 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act + modeledpnc"))), 
                data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                               tol = 1e-04)))
  model1 = lmer(as.formula(paste("X ~ ", paste(random.effects[[i]], "+ (1|id.x) + temp + relhum + gender + age + session + city + bmi + season + physical_act +modeledpnc + pncmedian"))), 
                data = data, REML = FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore",
                                                                                               tol = 1e-04)))
  vcov = as.data.frame(VarCorr(model1))$vcov
  res = c(summary(model1)$coefficients["pncmedian",
                                       1:2], anova(model0, model1)$'Pr(>Chisq)'[2],
          vcov[1]/sum(vcov), vcov[2]/sum(vcov))
  names(res) = c("coef", "coef.se", "pval", "subjectid",
                 "plate")
  return(res)
}

results = t(apply(data[,27:ncol(data)], 2, FUN = univ.fun))

pncmedian.results[[i]] <- data.frame(results)[,c(1,3)]

}



names(pm25_adj_combined.results) <- c("proteins","transcripts","metabolites","methylation")
names(pm25_adj_o.results) <- c("proteins","transcripts","metabolites","methylation")
names(pm25_adj_p.results) <- c("proteins","transcripts","metabolites","methylation")
names(pncmedian.results) <- c("proteins","transcripts","metabolites","methylation")

proteins.results <- list(pm25_adj_combined.results[[1]],pm25_adj_o.results[[1]],pm25_adj_p.results[[1]],pncmedian.results[[1]])
transcripts.results <- list(pm25_adj_combined.results[[2]],pm25_adj_o.results[[2]],pm25_adj_p.results[[2]],pncmedian.results[[2]])
metabolites.results <- list(pm25_adj_combined.results[[3]],pm25_adj_o.results[[3]],pm25_adj_p.results[[3]],pncmedian.results[[3]])
methylation.results <- list(pm25_adj_combined.results[[4]],pm25_adj_o.results[[4]],pm25_adj_p.results[[4]],pncmedian.results[[4]])

save(list=c("proteins.results", "transcripts.results", "metabolites.results","methylation.results"), file="univariate_results.RData")
