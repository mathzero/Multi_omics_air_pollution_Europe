getPLSModel<-function(X, Y, designId,modelName)
{ 
  
  modelList <- list()
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  mPLS <- pls(X = withinVarX, Y = withinVarY, ncomp=min(dim(X)[2],dim(Y)[2]))
  VIP <- vip(mPLS)
  modelList[[1]]<-withinVarX 
  modelList[[2]]<-withinVarY
  modelList[[3]]<-mPLS
  
  modelList[[4]]<-quality.var(mPLS, withinVarX, withinVarY)
  modelList[[5]] <- vip(mPLS)
  names(modelList)<- c("withinVarX", "withinVarY", modelName,"quality.var","VIP")   
  return(modelList)
}



quality.var <- function(model,X,Y){
  p1 <- predict(model,X) 
  matmean <- matrix(colMeans(Y),ncol=dim(Y)[2],nrow=dim(Y)[1],byrow=TRUE)
  sumYT <- colSums((Y-matmean)**2)
  cumP <- rep(0,model$ncomp)
  cumPY <- matrix(0,nrow=model$ncomp,ncol=dim(Y)[2])
  for (k in 1:model$ncomp){
    selectedY <- which(model$loadings$Y[,k]!=0)
    sumPLS <- colSums(as.data.frame(p1$predict[,selectedY,k]-matmean[selectedY])**2)
    cumPY[k,selectedY] <- sumPLS/sumYT[selectedY]
    cumP[k] <- mean(sumPLS/sumYT[selectedY]) 
  }
  IndP <- diff(c(0,cumP))
  IndPY <- apply(cumPY,2,FUN=function(x) diff(c(0,x)))
  if(dim(cumPY)[1]>1){
    colnames(cumPY) <- colnames(IndPY) <- colnames(Y)
    rownames(cumPY) <- rownames(IndPY) <- paste("comp",1:model$ncomp,sep="")
    names(cumP) <- names(IndP) <- paste("comp",1:model$ncomp,sep="")
  }else{
    colnames(cumPY) <- names(IndPY) <- colnames(Y)
    names(cumP) <- names(IndP) <- paste("comp",1:model$ncomp,sep="")
  }
  res <- list(CumPerExplain=cumP,PerExplain=IndP,CumPerExplainY=cumPY,PerExplainY=IndPY)
}


PerfSparseOnX<-function(X, Y, designId, niter=10)
{
  SelectedX<-NULL
  SelectedY<-NULL
  tunedX<- NULL
  Summary <- NULL
  MaxVarX <- dim(X)[2]
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  nComp <- min(dim(X)[2],dim(Y)[2])
  for(comp in c(1:nComp)){
    TmpSummary <- NULL
    for(i in 1:MaxVarX){
      TmpKeepX <- c(SelectedX,i)
      TmpKeepY <- rep(dim(Y)[2],comp)
      TmpsPLS <- spls(withinVarX,withinVarY,keepX=TmpKeepX,keepY=TmpKeepY,ncomp=comp,mode='regression')
      Q2_iter=MSEP_iter=NULL
      cat("\n")
      cat(paste0("\n", comp, " comp., ", i, " var.", "\n"))
      pb=txtProgressBar(style=3)
      for (k in 1:niter){
        setTxtProgressBar(pb, k/niter)
        TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=1, progressBar = FALSE)
        Q2_iter=cbind(Q2_iter, TmpPerf$Q2.total)
        MSEP_iter=abind(MSEP_iter, TmpPerf$MSEP, along=3)
      }
      Q2.total=matrix(apply(Q2_iter, 1, mean), ncol=comp)
      MSEP=apply(matrix(apply(MSEP_iter, 3, colSums), ncol=comp), 2, mean)
      TmpSummary <- rbind(TmpSummary,c(comp,i,c(Q2.total,rep(NA,nComp-comp)),c(MSEP,rep(NA,nComp-comp))))
    }
    TmpTune=tune.splslevel(X = X, Y = Y, multilevel = as.matrix(as.numeric(as.factor(designId), ncol = 1)), 
                           ncomp = comp, mode = "regression", already.tested.X = SelectedX, already.tested.Y = SelectedY,
                           test.keepX = seq(1, ncol(X)), test.keepY = ncol(Y))
    NVarIn <- which.max(TmpTune$cor.value)
    TmpSel <- rep(0,MaxVarX)
    TmpSel[NVarIn] <- 1
    TmpSummary <- as.data.frame(cbind(TmpSummary,TmpTune$cor.value,TmpSel))
    SelectedX <- c(SelectedX,as.numeric(NVarIn))
    SelectedY <- c(SelectedY,dim(Y)[2])
    colnames(TmpSummary) <- c('nComp','nVar',paste('Q2_C',c(1:nComp),sep=''),paste('sum_MSEP',c(1:nComp)),'max_Corr','best')
    Summary <- rbind(Summary,TmpSummary)
    if(comp==1){
      colnames(Summary) <- c('nComp','nVar',paste('Q2_C',c(1:nComp),sep=''),paste('sum_MSEP',c(1:nComp)),'max_Corr','best')
    }
  }
  
  return(Summary)
}

getSPLSModelonX<-function(X, Y, designId,NCompX,NVarX)
{ 
  
  modelList <- list()
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  mPLS <- spls(X = withinVarX, Y = withinVarY, ncomp=NCompX,keepX=NVarX)
  VIP <- vip(mPLS)
  modelList[[1]]<-withinVarX 
  modelList[[2]]<-withinVarY
  modelList[[3]]<-mPLS
  
  modelList[[4]] <- VIP
  modelList[[5]]<-quality.var(mPLS, withinVarX, withinVarY)
    
  return(modelList)
}

PerfSparseOnY<-function(X, Y, designId, niter=10)
{
  SelectedX<-NULL
  SelectedY<-NULL
  tunedY<- NULL
  Summary <- NULL
  MaxVarX <- dim(X)[2]
  MaxVarY <- dim(Y)[2]
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  nComp <- min(dim(X)[2],dim(Y)[2])
  for(comp in c(1:nComp)){
    TmpSummary <- NULL
    for(i in 1:MaxVarY){
      TmpKeepX <- rep(dim(X)[2], comp)
      TmpKeepY <- c(SelectedY, i)
      TmpsPLS <- spls(withinVarX,withinVarY,keepX=TmpKeepX,keepY=TmpKeepY,ncomp=comp,mode='regression')
      Q2_iter=MSEP_iter=NULL
      cat("\n")
      cat(paste0("\n", comp, " comp., ", i, " var.", "\n"))
      pb=txtProgressBar(style=3)
      for (k in 1:niter){
        setTxtProgressBar(pb, k/niter)
        TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=1, progressBar = FALSE)
        Q2_iter=cbind(Q2_iter, TmpPerf$Q2.total)
        MSEP_iter=abind(MSEP_iter, TmpPerf$MSEP, along=3)
      }
      Q2.total=matrix(apply(Q2_iter, 1, mean), ncol=comp)
      MSEP=apply(matrix(apply(MSEP_iter, 3, colSums), ncol=comp), 2, mean)
      TmpSummary <- rbind(TmpSummary,c(comp,i,c(Q2.total,rep(NA,nComp-comp)),c(MSEP,rep(NA,nComp-comp))))
    }
    TmpTune=tune.splslevel(X = X, Y = Y, multilevel = as.matrix(as.numeric(as.factor(designId), ncol = 1)), 
                           ncomp = comp, mode = "regression", already.tested.X = SelectedX, already.tested.Y = SelectedY,
                           test.keepX = ncol(X), test.keepY = seq(1, ncol(Y)))
    NVarIn <- which.max(TmpTune$cor.value)
    TmpSel <- rep(0,MaxVarY)
    TmpSel[NVarIn] <- 1
    TmpSummary <- as.data.frame(cbind(TmpSummary,t(TmpTune$cor.value),TmpSel))
    SelectedX <- c(SelectedX,dim(X)[2])
    SelectedY <- c(SelectedY,as.numeric(NVarIn))
    colnames(TmpSummary) <- c('nComp','nVar',paste('Q2_C',c(1:nComp),sep=''),paste('sum_MSEP',c(1:nComp)),'max_Corr','best')
    Summary <- rbind(Summary,TmpSummary)
    if(comp==1){
      colnames(Summary) <- c('nComp','nVar',paste('Q2_C',c(1:nComp),sep=''),paste('sum_MSEP',c(1:nComp)),'max_Corr','best')
    }
  }
  
  return(Summary)
}


getSPLSModelonY<-function(X, Y, designId,NCompY,NVarY)
{ 
  
  modelList <- list()
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  mPLS <- spls(X = withinVarX, Y = withinVarY, ncomp=NCompY,keepY = NVarY)
  VIP <- vip(mPLS)
  modelList[[1]]<-withinVarX 
  modelList[[2]]<-withinVarY
  modelList[[3]]<-mPLS
  
  modelList[[4]]<-quality.var(mPLS, withinVarX, withinVarY)
  modelList[[5]] <- vip(mPLS)
  
  return(modelList)
}

PerfSparseOnXAndY<-function(X, Y, designId, niter=10)
{
  SelectedX<-NULL
  SelectedY<-NULL
  tunedX<- NULL
  tunedY<- NULL
  Summary <- NULL
  MaxVarX <- dim(X)[2]
  MaxVarY <- dim(Y)[2]
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  nComp <- min(dim(X)[2],dim(Y)[2])
  for(comp in c(1:nComp)){
    TmpSummary <- NULL
    for(i in 1:MaxVarX){
      TmpKeepX <- c(SelectedX,i)
      for (j in 1:MaxVarY){
        TmpKeepY <- c(SelectedY, i)
        TmpsPLS <- spls(withinVarX,withinVarY,keepX=TmpKeepX,keepY=TmpKeepY,ncomp=comp,mode='regression')
        Q2_iter=MSEP_iter=NULL
        cat("\n")
        cat(paste0("\n", comp, " comp., ", i, ";", j, " var.", "\n"))
        pb=txtProgressBar(style=3)
        for (k in 1:niter){
          setTxtProgressBar(pb, k/niter)
          TmpPerf <- perf(TmpsPLS,validation='Mfold',folds=5,nrepeat=1, progressBar = FALSE)
          Q2_iter=cbind(Q2_iter, TmpPerf$Q2.total)
          MSEP_iter=abind(MSEP_iter, TmpPerf$MSEP, along=3)
        }
        Q2.total=matrix(apply(Q2_iter, 1, mean), ncol=comp)
        MSEP=apply(matrix(apply(MSEP_iter, 3, colSums), ncol=comp), 2, mean)
        TmpSummary <- rbind(TmpSummary,c(comp,i,j,c(Q2.total,rep(NA,nComp-comp)),c(MSEP,rep(NA,nComp-comp))))
      }
    }
    TmpTune=tune.splslevel(X = X, Y = Y, multilevel = as.matrix(as.numeric(as.factor(designId), ncol = 1)), 
                           ncomp = comp, mode = "regression", already.tested.X = SelectedX, already.tested.Y = SelectedY,
                           test.keepX = seq(1, ncol(X)), test.keepY = seq(1, ncol(Y)))
    NVarIn <- which.max(as.vector(t(TmpTune$cor.value)))
    TmpSel <- rep(0,MaxVarX*MaxVarY)
    TmpSel[NVarIn] <- 1
    TmpSummary <- as.data.frame(cbind(TmpSummary,as.vector(t(TmpTune$cor.value)),TmpSel))
    colnames(TmpSummary) <- c('nComp','nVarX', 'nVarY', paste('Q2_C',c(1:nComp),sep=''),paste('sum_MSEP',c(1:nComp)),'max_Corr','best')
    SelectedX <- c(SelectedX, TmpSummary[NVarIn,]$nVarX)
    SelectedY <- c(SelectedY, TmpSummary[NVarIn,]$nVarY)
    Summary <- rbind(Summary,TmpSummary)
    if(comp==1){
      colnames(Summary) <- c('nComp','nVarX', 'nVarY',paste('Q2_C',c(1:nComp),sep=''),paste('sum_MSEP',c(1:nComp)),'max_Corr','best')
    }
  }
  
  return(Summary)
}

getSPLSModelonXAndY<-function(X, Y, designId, NComp, NVarX, NVarY)
{ 
  
  modelList <- list()
  withinVarX <- withinVariation(X,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  withinVarY <- withinVariation(Y,design=as.matrix(as.numeric(as.factor(designId),ncol=1)))
  mPLS <- spls(X = withinVarX, Y = withinVarY, ncomp=NComp, keepX=NVarX, keepY = NVarY)
  VIP <- vip(mPLS)
  modelList[[1]]<-withinVarX 
  modelList[[2]]<-withinVarY
  modelList[[3]]<-mPLS
  
  modelList[[4]]<-quality.var(mPLS, withinVarX, withinVarY)
  modelList[[5]] <- vip(mPLS)
  
  return(modelList)
}