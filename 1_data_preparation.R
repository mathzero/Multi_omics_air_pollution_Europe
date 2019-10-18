### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"


library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(mice)

###########################################################################
###########################################################################

### Reading in data ###

exposures <- readRDS("data_clean/Exposures.rds")
proteins <- readRDS("data_clean/Proteins.rds")
metabolites <- readRDS("data_clean/Metabolites.rds")
covariates <- readRDS("data_clean/Covariates.rds")
transcripts <- readRDS("data_clean/Transcripts.rds")
methylation1 <- readRDS("data_clean/Methylation.rds")


###########################################################################
###########################################################################

### creating subjectidp variable ###

# metabolites$subjectidp <- rownames(metabolites)
# metabolites <- metabolites %>%
#   select("subjectidp", everything())
# 
# methylation$subjectidp <- rownames(methylation)
# methylation <- methylation %>%
#   select("subjectidp", everything())
# 
# proteins$subjectidp <- rownames(proteins)
# proteins <- proteins %>%
#   select("subjectidp", everything())
# 
# transcripts$subjectidp <- rownames(transcripts)
# transcripts <- transcripts %>%
#   select("subjectidp", everything())


###########################################################################
###########################################################################

### examining missing data ###
# 
omics <- list(proteins,transcripts,metabolites,methylation)
omic.names <- c("proteins","transcripts","metabolites","methylation")
na_counts <- list()
# 
# for (i in 1:length(omics)){
#   nac <-sapply(omics[[i]], function(y) sum(length(which(is.na(y)))))
#   nac <- data.frame(nac)
#   colnames(nac) <- "count"
#   nac$proportionNA <- round(nac$count / nrow(omics[[i]]),3)
#   nac <- nac[order(nac$count),] 
#   omicname <- omic.names[[i]]
#   png(paste0(omicname,"_missing_data.png", units="px", width=2000, height=1500, res=300))
#   plot(sort(nac[,2], decreasing = TRUE), col ="red", type = 'h', xlab = "index of omic", ylab="proportion of missing data", 
#        main = paste0("Propotion of missing data:",omicname))
#   dev.off()
#   na_counts[[i]] <- nac
# }


###########################################################################
###########################################################################

### Removing > certain proportion of NAs (10%) ####
imputed.omics <- list()

for (i in 1:length(omics)){
  data <- omics[[i]]
  data <- Filter(function(x) sum(is.na(x)) < length(x)/10, data)
  for(j in 1:ncol(data)){
    data[is.na(data[,j]), j] <- mean(data[,j], na.rm = TRUE)
  }
  imputed.omics[[i]] <- data
  saveRDS(data, file = paste0(omic.names[[i]],"_meanimpute.rds"))
}

proteins <- imputed.omics[[1]]
transcripts <- imputed.omics[[2]]
metabolites <- imputed.omics[[3]]
methylation <- imputed.omics[[4]]

### Removing observaions with lots of missing data (exposures) ###

exposures <- exposures[-which(rowMeans(is.na(exposures)) > 0.5), ]

sum(is.na(transcripts))
###########################################################################
###########################################################################

####### Imputation of covariates and exposures ####

mice_exposures_trimmed1=mice(exposures[,2:ncol(exposures)],method="cart",m=1,maxit=10)

exposures_imputed_mice=complete(mice_exposures_trimmed1)

###### Subsetting our dataset to get the final cleaned exposures

exposures = data.frame(subjectidp=exposures$subjectidp,exposures_imputed_mice[,c("id", "pm25_adj_p",
                                                                                 "pm25_adj_o", "pncmedian", "modeledpnc", "modeledpm25abs")])

# Imputing covariates
covariates_mice=mice(covariates[,2:14],method="cart",m=1,maxit=10)
covariates[,2:14]=complete(covariates_mice)

###########################################################################
###########################################################################

# log transform proteins
proteins[,2:ncol(proteins)] <-  log(proteins[,2:ncol(proteins)])

#### Scale all data ###
covariates[,c("temp","relhum","age","bmi","physical_act","education_adj")] <- scale(covariates[,c("temp","relhum","age","bmi","physical_act","education_adj")])
exposures[,c(3:ncol(exposures))] <- scale(exposures[,c(3:ncol(exposures))])

metabolites[,2:ncol(metabolites)] <- scale(metabolites[,2:ncol(metabolites)])
transcripts[,2:ncol(transcripts)] <- scale(transcripts[,2:ncol(transcripts)])
proteins[,2:ncol(proteins)] <- scale(proteins[,2:ncol(proteins)])
methylation[,2:ncol(methylation)] <- scale(methylation[,2:ncol(methylation)])


###########################################################################
###########################################################################

saveRDS(exposures,"exposures.rds")
saveRDS(covariates,"covariates.rds")
saveRDS(transcripts,"transcripts.rds")
saveRDS(metabolites,"metabolites.rds")
saveRDS(methylation,"methylation.rds")
saveRDS(proteins,"proteins.rds")



