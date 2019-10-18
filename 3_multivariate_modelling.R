### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"

library(tidyverse)
library(ggplot2)



exposures <- readRDS("data_denoised/exposures.rds")
proteins <- readRDS("data_denoised/proteins_denoised.rds")
metabolites <- readRDS("data_denoised/metabolites_denoised.rds")
covariates <- readRDS("data_denoised/covariates.rds")
transcripts <- readRDS("data_denoised/transcripts_denoised.rds")
methylation <- readRDS("data_denoised/methylation_denoised.rds")

###########################################################################
###########################################################################

data.prot <- inner_join(exposures, proteins, by = "subjectidp")

mod.prot <- lm()