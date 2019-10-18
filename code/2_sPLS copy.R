library(mixOmics)
library(RColorBrewer)
library('ggplot2')
library("corrplot")
library(abind)
library(mlbench)
library(doBy)
library(parallel)
library(foreach)
library(doParallel)
library(party)


### Your working directory on the server
work_dir<-"~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project"

setwd(work_dir)
#("~/Google Drive/Imperial/2 Computational epidemiology/Comp_Epi_Project")
#"/rds/general/user/mw418/home/comp_epi"
source('scripts/functions_iterations.R')

#####################################################################################
#####################################################################################
### PROTEINS ########################################################################
#####################################################################################
#####################################################################################


Y_exposure <- readRDS("data_denoised/exposures.rds")
X_protein <- readRDS("data_denoised/proteins_denoised_noID.rds")

#making ID vector
X_protein$id <- substr(X_protein$subjectidp,1,3)
table(as.factor(X_protein$id))
X_protein <- X_protein[,c(15,1:14)]

#remove samples with only 1 measurement
X_protein <- X_protein[duplicated(X_protein$id) | duplicated(X_protein$id, fromLast = TRUE),]
table(as.factor(X_protein$id))

#Matching Exposure[Y] to Metabolites[X]
Y_exposure <- Y_exposure[Y_exposure[,1] %in% X_protein[,2], ]

#Making the rownames
Y_exposure$subjectidp <- paste0(substr(Y_exposure$subjectidp,1,3), "_", 
                                substr(Y_exposure$subjectidp,4,4), sep="")
rownames(Y_exposure) <- Y_exposure$subjectidp

X_protein$subjectidp <- paste0(substr(X_protein$subjectidp, 1, 3), "_",
                               substr(X_protein$subjectidp, 4, 4), sep = "")
rownames(X_protein) <- X_protein$subjectidp

#Design matrix
design_protein <- as.matrix(substr(Y_exposure$subjectidp,1,3))

#Removing unwanted columns
Y_exposure$pnc_diff <- Y_exposure$modeledpnc - Y_exposure$pncmedian
Y_exposure <- Y_exposure[,c(3:5,8)]

X_protein <- X_protein[,3:ncol(X_protein)]

####################### sPLS MODEL ########################

#[Y] pm25_p Model
Summ_Xprot.Ypm25_p <- PerfSparseOnX(X_protein, data.frame(Y_exposure[,1]), design_protein, niter=100)
KeepXFull <- Summ_Xprot.Ypm25_p[Summ_Xprot.Ypm25_p$best==1, 'nVar']
sPLS_Xprot.Ypm25_p <- getSPLSModelonX(X_protein, data.frame(Y_exposure[,1]), design_protein, NCompX = 1, KeepXFull)

#[Y] pm25_o Model
Summ_Xprot.Ypm25_o <- PerfSparseOnX(X_protein, data.frame(Y_exposure[,2]), design_protein, niter=100)
KeepXFull <- Summ_Xprot.Ypm25_o[Summ_Xprot.Ypm25_o$best==1, 'nVar']
sPLS_Xprot.Ypm25_o <- getSPLSModelonX(X_protein, data.frame(Y_exposure[,2]), design_protein, NCompX = 1, KeepXFull)

#[Y] pnc_median Model
Summ_Xprot.Ypnc <- PerfSparseOnX(X_protein, data.frame(Y_exposure[,3]), design_protein, niter=100)
KeepXFull <- Summ_Xprot.Ypnc[Summ_Xprot.Ypnc$best==1, 'nVar']
sPLS_Xprot.Ypnc <- getSPLSModelonX(X_protein, data.frame(Y_exposure[,3]), design_protein, NCompX = 1, KeepXFull)

#[Y] pnc_diff Model
Summ_Xprot.Ypnc_diff <- PerfSparseOnX(X_protein, data.frame(Y_exposure[,4]), design_protein, niter=100)
KeepXFull <- Summ_Xprot.Ypnc_diff[Summ_Xprot.Ypnc_diff$best==1, 'nVar']
sPLS_Xprot.Ypnc_diff <- getSPLSModelonX(X_protein, data.frame(Y_exposure[,4]), design_protein, NCompX = 1, KeepXFull)

Summ_protein <- list(pm25_p = Summ_Xprot.Ypm25_p, pm25_o = Summ_Xprot.Ypm25_o, pnc = Summ_Xprot.Ypnc, pnc_diff = Summ_Xprot.Ypnc_diff)
sPLS_protein <- list(pm25_p=sPLS_Xprot.Ypm25_p, pm25_o=sPLS_Xprot.Ypm25_o, pnc=sPLS_Xprot.Ypnc, pnc_diff=sPLS_Xprot.Ypnc_diff)

save(Summ_protein, sPLS_protein, file = "protein_PLS.RData")


#####################################################################################
#####################################################################################
### METABOLITES #####################################################################
#####################################################################################
#####################################################################################


Y_exposure <- readRDS("data_denoised/exposures.rds")
X_metab <- readRDS("data_denoised/metabolites_denoised_noID.rds")

#making ID vector
id <- X_metab$subjectidp
id <- substr(id, 1,3)
table(as.factor(id))

#adding ID vector back to X_metab
X_metab$id <- id
X_metab <- X_metab[,c(3631,1:3630)]
View(X_metab[,1:2])      

#remove samples with only 1 measurement
X_metab <- X_metab[duplicated(X_metab$id) | duplicated(X_metab$id, fromLast = TRUE),]
table(as.factor(X_metab$id))

#Matching Exposure[Y] to Metabolites[X]
Y_exposure <- Y_exposure[Y_exposure[,1] %in% X_metab[,2], ]

#Making the rownames
Y_exposure$subjectidp <- paste0(substr(Y_exposure$subjectidp,1,3), "_", 
                                substr(Y_exposure$subjectidp,4,4), sep="")
rownames(Y_exposure) <- Y_exposure$subjectidp

X_metab$subjectidp <- paste0(substr(X_metab$subjectidp, 1, 3), "_",
                             substr(X_metab$subjectidp, 4, 4), sep = "")
rownames(X_metab) <- X_metab$subjectidp

#Design matrix
design_metab <- as.matrix(substr(Y_exposure$subjectidp,1,3))

#Removing unwanted columns
Y_exposure$pnc_diff <- Y_exposure$modeledpnc - Y_exposure$pncmedian
Y_exposure <- Y_exposure[,c(3:5,8)]

X_metab <- X_metab[,3:ncol(X_metab)]

########## sPLS MODEL ##############

n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)

t0=Sys.time()


#[Y] pm25_p Model
Summ_Xmetab.Ypm25_P <- PerfSparseOnX(X_metab, data.frame(Y_exposure[,1]), design_metab,
                                     niter=100)
KeepXFull <- Summ_Xmetab.Ypm25_P[Summ_Xmetab.Ypm25_P$best==1, 'nVar']
sPLS_Xmetab.Ypm25_P <- getSPLSModelonX(X_metab, data.frame(Y_exposure[,1]), design_metab,
                                       NCompX = 1, KeepXFull)

#[Y] pm25_o Model
Summ_Xmetab.Ypm25_O <- PerfSparseOnX(X_metab, data.frame(Y_exposure[,2]), design_metab,
                                     niter=100)
KeepXFull <- Summ_Xmetab.Ypm25_O[Summ_Xmetab.Ypm25_O$best==1, 'nVar']
sPLS_Xmetab.Ypm25_O <- getSPLSModelonX(X_metab, data.frame(Y_exposure[,2]), design_metab,
                                       NCompX = 1, KeepXFull)

#[Y] pnc_median Model
Summ_Xmetab.Ypnc <- PerfSparseOnX(X_metab, data.frame(Y_exposure[,3]), design_metab,
                                  niter=100)
KeepXFull <- Summ_Xmetab.Ypnc[Summ_Xmetab.Ypnc$best==1, 'nVar']
sPLS_Xmetab.Ypnc <- getSPLSModelonX(X_metab, data.frame(Y_exposure[,3]), design_metab,
                                    NCompX = 1, KeepXFull)

#[Y] pnc_diff Model
Summ_Xmetab.Ypnc_diff <- PerfSparseOnX(X_metab, data.frame(Y_exposure[,4]), design_metab,
                                       niter=100)
KeepXFull <- Summ_Xmetab.Ypnc_diff[Summ_Xmetab.Ypnc_diff$best==1, 'nVar']
sPLS_Xmetab.Ypnc_diff <- getSPLSModelonX(X_metab, data.frame(Y_exposure[,4]), design_metab,
                                         NCompX = 1, KeepXFull)

stopCluster(cl)
t1=Sys.time()
print(t1-t0)

Summ_metabolite <- list(pm25_P = Summ_Xmetab.Ypm25_P, pm25_O = Summ_Xmetab.Ypm25_O, pnc = Summ_Xmetab.Ypnc, pnc_diff = Summ_Xmetab.Ypnc_diff)
sPLS_metabolite <- list(pm25_P = sPLS_Xmetab.Ypm25_P, pm25_O = sPLS_Xmetab.Ypm25_O, pnc = sPLS_Xmetab.Ypnc, pnc_diff = sPLS_Xmetab.Ypnc_diff)

save(Summ_metabolite, sPLS_metabolite, file = "metabolite_PLS.RData")
############################################


#####################################################################################
#####################################################################################
### TRANSCRIPTS #####################################################################
#####################################################################################
#####################################################################################


Y_exposure <- readRDS("data_denoised/exposures.rds")
X_trans <- readRDS("data_denoised/transcripts_denoised_noID.rds")

#making ID vector
id <- X_trans$subjectidp
id <- substr(id,1,3)
table(as.factor(id))

#adding ID vector back to X_trans
X_trans$id <- id
X_trans <- X_trans[,c(23559, 1:23558)]
View(X_trans[,1:2]) 

#remove samples w/ only 1 measurement
X_trans <- X_trans[duplicated(X_trans$id) | duplicated(X_trans$id, fromLast = TRUE),]
table(as.factor(X_trans$id))

#Matching Exposure[Y] to transcripts[X]
Y_exposure <- Y_exposure[Y_exposure[,1] %in% X_trans[,2], ]

#Making the rownames
Y_exposure$subjectidp <- paste0(substr(Y_exposure$subjectidp,1,3), "_", 
                                substr(Y_exposure$subjectidp,4,4), sep="")
rownames(Y_exposure) <- Y_exposure$subjectidp

X_trans$subjectidp <- paste0(substr(X_trans$subjectidp, 1, 3), "_",
                             substr(X_trans$subjectidp, 4, 4), sep = "")
rownames(X_trans) <- X_trans$subjectidp

#Design matrix
design_trans <- as.matrix(substr(Y_exposure$subjectidp,1,3))

#Removing unwanted columns
Y_exposure$pnc_diff <- Y_exposure$modeledpnc - Y_exposure$pncmedian
Y_exposure <- Y_exposure[,c(3:5,8)]

X_trans <- X_trans[,3:ncol(X_trans)]

########## sPLS MODEL ##############


n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)

t0=Sys.time()

#[Y] pm25_p Model
Summ_Xtrans.Ypm25_P <- PerfSparseOnX(X_trans, data.frame(Y_exposure[,1]), design_trans,
                                     niter=100)
KeepXFull <- Summ_Xtrans.Ypm25_P[Summ_Xtrans.Ypm25_P$best==1, 'nVar']
sPLS_Xtrans.Ypm25_P <- getSPLSModelonX(X_trans, data.frame(Y_exposure[,1]), design_trans,
                                       NCompX = 1, KeepXFull)

#[Y] pm25_o Model
Summ_Xtrans.Ypm25_O <- PerfSparseOnX(X_trans, data.frame(Y_exposure[,2]), design_trans,
                                     niter=100)
KeepXFull <- Summ_Xtrans.Ypm25_O[Summ_Xtrans.Ypm25_O$best==1, 'nVar']
sPLS_Xtrans.Ypm25_O <- getSPLSModelonX(X_trans, data.frame(Y_exposure[,2]), design_trans,
                                       NCompX = 1, KeepXFull)

#[Y] pnc_median Model
Summ_Xtrans.Ypnc <- PerfSparseOnX(X_trans, data.frame(Y_exposure[,3]), design_trans, niter=100)
KeepXFull <- Summ_Xtrans.Ypnc[Summ_Xtrans.Ypnc$best==1, 'nVar']
sPLS_Xtrans.Ypnc <- getSPLSModelonX(X_trans, data.frame(Y_exposure[,3]), design_trans, NCompX = 1, KeepXFull)

#[Y] pnc_diff Model
Summ_Xtrans.Ypnc_diff <- PerfSparseOnX(X_trans, data.frame(Y_exposure[,4]), design_trans, niter=100)
KeepXFull <- Summ_Xtrans.Ypnc_diff[Summ_Xtrans.Ypnc_diff$best==1, 'nVar']
sPLS_Xtrans.Ypnc_diff <- getSPLSModelonX(X_trans, data.frame(Y_exposure[,4]), design_trans, NCompX = 1, KeepXFull)

stopCluster(cl)
t1=Sys.time()
print(t1-t0)

Summ_transcript <- list(pm25_P = Summ_Xtrans.Ypm25_P, pm25_O = Summ_Xtrans.Ypm25_O, pnc = Summ_Xtrans.Ypnc, pnc_diff = Summ_Xtrans.Ypnc_diff)
sPLS_transcript <- list(pm25_P = sPLS_Xtrans.Ypm25_P, pm25_O = sPLS_Xtrans.Ypm25_O, pnc = sPLS_Xtrans.Ypnc, pnc_diff = sPLS_Xtrans.Ypnc_diff)

save(Summ_transcript, sPLS_transcript, file = "transcript_PLS.RData")

############################################

#####################################################################################
#####################################################################################
### METHYLATION #####################################################################
#####################################################################################
#####################################################################################



Y_exposure <- readRDS("data_denoised/exposures.rds")
X_methyl <- readRDS("data_denoised/methylation_denoised_noID.rds")

#making ID Vector
id <- X_methyl$subjectidp
id <- substr(id,1,3)
table(as.factor(id))

#adding ID vector back to X_methyl
X_methyl$id <- id
X_methyl <- X_methyl[,c(448855, 1:448854)]

#remove samples w/ only 1 measurement
X_methyl <- X_methyl[duplicated(X_methyl$id) | duplicated(X_methyl$id, fromLast = TRUE),]
table(as.factor(X_methyl$id))

#Matching Exposure[Y] to methylation[X]
Y_exposure <- Y_exposure[Y_exposure[,1] %in% X_methyl[,2],]

#Making the rownames
Y_exposure$subjectidp <- paste0(substr(Y_exposure$subjectidp,1,3), "_",
                                substr(Y_exposure$subjectidp,4,4), sep="")
rownames(Y_exposure) <- Y_exposure$subjectidp

X_methyl$subjectidp <- paste0(substr(X_methyl$subjectidp,1,3), "_",
                              substr(X_methyl$subjectidp,4,4), sep="")
rownames(X_methyl) <- X_methyl$subjectidp

#Design matrix
design_methyl <- as.matrix(substr(Y_exposure$subjectidp,1,3))

#removing unwanted columns
Y_exposure$pnc_diff <- Y_exposure$modeledpnc - Y_exposure$pncmedian
Y_exposure <- Y_exposure[,c(3:5,8)]

X_methyl <- X_methyl[,3:ncol(X_methyl)]

########## sPLS MODEL ##############

n_cores=detectCores()
cl <- makeCluster(n_cores-1) #not to overload your computer
registerDoParallel(cl)

t0=Sys.time()


#[Y] pm25_p Model
Summ_Xmethyl.Ypm25_P <- PerfSparseOnX(X_methyl, data.frame(Y_exposure[,1]), design_methyl,
                                      niter=100)
KeepXFull <- Summ_Xmethyl.Ypm25_P[Summ_Xmethyl.Ypm25_P$best==1, 'nVar']
sPLS_Xmethyl.Ypm25_P <- getSPLSModelonX(X_methyl, data.frame(Y_exposure[,1]), design_methyl,
                                        NCompX = 1, KeepXFull)

#[Y] pm25_o Model
Summ_Xmethyl.Ypm25_O <- PerfSparseOnX(X_methyl, data.frame(Y_exposure[,2]), design_methyl,
                                      niter=100)
KeepXFull <- Summ_Xmethyl.Ypm25_O[Summ_Xmethyl.Ypm25_O$best==1, 'nVar']
sPLS_Xmethyl.Ypm25_O <- getSPLSModelonX(X_methyl, data.frame(Y_exposure[,2]), design_methyl,
                                        NCompX = 1, KeepXFull)

#[Y] pnc_median Model
Summ_Xmethyl.Ypnc <- PerfSparseOnX(X_methyl, data.frame(Y_exposure[,3]), design_methyl, niter=100)
KeepXFull <- Summ_Xmethyl.Ypnc[Summ_Xmethyl.Ypnc$best==1, 'nVar']
sPLS_Xmethyl.Ypnc <- getSPLSModelonX(X_methyl, data.frame(Y_exposure[,3]), design_methyl, NCompX = 1, KeepXFull)

#[Y] pnc_diff Model
Summ_Xmethyl.Ypnc_diff <- PerfSparseOnX(X_methyl, data.frame(Y_exposure[,4]), design_methyl, niter=100)
KeepXFull <- Summ_Xmethyl.Ypnc_diff[Summ_Xmethyl.Ypnc_diff$best==1, 'nVar']
sPLS_Xmethyl.Ypnc_diff <- getSPLSModelonX(X_methyl, data.frame(Y_exposure[,4]), design_methyl, NCompX = 1, KeepXFull)


stopCluster(cl)
t1=Sys.time()
print(t1-t0)

Summ_methyl <- list(pm25_P = Summ_Xmethyl.Ypm25_P, pm25_O = Summ_Xmethyl.Ypm25_O, pnc = Summ_Xmethyl.Ypnc, pnc_diff = Summ_Xmethyl.Ypnc_diff)
sPLS_methyl <- list(pm25_P = sPLS_Xmethyl.Ypm25_P, pm25_O = sPLS_Xmethyl.Ypm25_O, pnc = sPLS_Xmethyl.Ypnc, pnc_diff = sPLS_Xmethyl.Ypnc_diff)

save(Summ_methyl, sPLS_methyl, file = "methylation_PLS.RData")