#Step 0: input
filepath <- paste0("/home/proj/jeffwu/isye-ae/")         # set working directory if executed from somewhere else besides directory with RawData
#filepath <- ""
#setwd(filepath)
num.case <- 30
available_file.tmp <- 0:998
TT <- length(available_file.tmp)
numPart <- 10
dimPart <- num.case/numPart * TT
indvec = 7
check.tm.num = 20
numModes <- 300
cutoff.numModes = 300
POD.case = 1

#load libraries
library(Rcpp)
library(bigmemory)
library(BH)
library(calibrator)
library(bigpca)
library(rARPACK)
library(parallel)
library(colorRamps)
library(shape)

sourceCpp(code='
          #include <Rcpp.h>
          using namespace Rcpp;
          
          #include <algorithm>
          // [[Rcpp::depends(BH, bigmemory)]]
          #include <bigmemory/MatrixAccessor.hpp>
          #include <numeric>
          #include <iostream>
          #include <omp.h>
          
          //  function for bigInnerProd
          //
          // [[Rcpp::plugins(openmp)]]
          template <typename T>
          NumericMatrix bigInnerProd(XPtr<BigMatrix> pMat, MatrixAccessor<T> mat) {
          //sym - flag for whether matrix is symmetric
          NumericMatrix ret(pMat->ncol(),pMat->ncol());
          #pragma omp parallel for
          for (size_t i=0; i < pMat->ncol(); ++i){
          for (size_t j=i; j < pMat->ncol(); ++j){
          //      cout << "(i,j): " << i << ", " << j << endl;
          //      cout << mat[i] << endl;
          //      cout << mat[j] << endl;
          //      cout << std::inner_product(mat[i], mat[i]+pMat->nrow(), mat[j], 0.0) << endl;
          ret(i,j) = std::inner_product(mat[i], mat[i]+pMat->nrow(), mat[j], 0.0);
          }
          }
          return (ret);
          }
          
          // Dispatch function for bigInnerProd
          //
          
          // [[Rcpp::export]]
          NumericVector bigInnerProd(SEXP pBigMat) {
          // First we have to tell Rcpp what class to use for big.matrix objects.
          // This object stores the attributes of the big.matrix object passed to it
          // by R.
          XPtr<BigMatrix> xpMat(pBigMat);
          
          switch(xpMat->matrix_type()) {
          case 1:
          return bigInnerProd(xpMat, MatrixAccessor<char>(*xpMat));
          case 2:
          return bigInnerProd(xpMat, MatrixAccessor<short>(*xpMat));
          case 4:
          return bigInnerProd(xpMat, MatrixAccessor<int>(*xpMat));
          case 8:
          return bigInnerProd(xpMat, MatrixAccessor<double>(*xpMat));
          default:
          // This case should never be encountered unless the implementation of
          // big.matrix changes, but is necessary to implement shut up compiler
          // warnings.
          throw Rcpp::exception("unknown type detected for big.matrix object!");
          }
          }')

sourceCpp(code='
          #include <Rcpp.h>
          using namespace Rcpp;
          
          #include <algorithm>
          // [[Rcpp::depends(BH, bigmemory)]]
          #include <bigmemory/MatrixAccessor.hpp>
          #include <numeric>
          #include <iostream>
          #include <omp.h>
          
          //  function for bigInnerProd
          //
          // [[Rcpp::plugins(openmp)]]
          template <typename T>
          NumericMatrix bigInnerProd2(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat2, MatrixAccessor<T> mat, MatrixAccessor<T> mat2) {
          NumericMatrix ret(pMat->ncol(),pMat2->ncol());
          #pragma omp parallel for
          for (size_t i=0; i < pMat->ncol(); ++i){
          for (size_t j=0; j < pMat2->ncol(); ++j){
          //      cout << "(i,j): " << i << ", " << j << endl;
          //      cout << mat[i] << endl;
          //      cout << mat[j] << endl;
          //      cout << std::inner_product(mat[i], mat[i]+pMat->nrow(), mat[j], 0.0) << endl;
          ret(i,j) = std::inner_product(mat[i], mat[i]+pMat->nrow(), mat2[j], 0.0);
          }
          }
          return (ret);
          }
          
          // Dispatch function for bigInnerProd
          //
          
          // [[Rcpp::export]]
          NumericMatrix bigInnerProd2(SEXP pBigMat, SEXP pBigMat2) {
          // First we have to tell Rcpp what class to use for big.matrix objects.
          // This object stores the attributes of the big.matrix object passed to it
          // by R.
          XPtr<BigMatrix> xpMat(pBigMat);
          XPtr<BigMatrix> xpMat2(pBigMat2);
          
          switch(xpMat->matrix_type()) {
          case 1:
          return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<char>(*xpMat), MatrixAccessor<char>(*xpMat2));
          case 2:
          return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<short>(*xpMat), MatrixAccessor<short>(*xpMat2));
          case 4:
          return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<int>(*xpMat), MatrixAccessor<int>(*xpMat2));
          case 8:
          return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<double>(*xpMat), MatrixAccessor<double>(*xpMat2));
          default:
          // This case should never be encountered unless the implementation of
          // big.matrix changes, but is necessary to implement shut up compiler
          // warnings.
          throw Rcpp::exception("unknown type detected for big.matrix object!");
          }
          }')

#THE FOLLOWING IS JUST A TEMPLATE. THE RUN FILES ARE THREADED ON SCRIPTS

#K11.R
############################################################################
#Set working directory where data is
for(ind in indvec){
  cat("ind =", ind, "\n")
  setwd(paste0(filepath,"AE test for ICA/", ind))
  for (i in 1:numPart){
    for (j in i:numPart){
	cat("i = ",i, "j = ", j, "\n")
      #Loading libraries
      FTtmp1 <- attach.big.matrix(paste0('FTtmp',i,'.dsc'))
      FTtmp2 <- attach.big.matrix(paste0('FTtmp',j,'.dsc'))
      
      #Warmup
      #mat1 <- as.big.matrix(FTtmp1[,1:5])
      #tmp <- bigInnerProd2(mat1@address,FTtmp1@address)
      #mat2 <- as.big.matrix(FTtmp2[,1:5])
      #tmp <- bigInnerProd2(mat2@address,FTtmp2@address)
      
      # (3.7 sec with warmup, 7.5 with no warmup)
      ptm <- proc.time()
      Kpart <- bigInnerProd2(FTtmp1@address,FTtmp2@address)
      save(Kpart,file=paste0('K',i,j,'.RData'))
    }
  }
  K = matrix(0,nrow=TT*num.case,ncol=TT*num.case)
  for (i in 1:numPart){
    for (j in i:numPart){
      load(paste0("K",i,j,".RData"))
      #Kpart is the name
      K[ ((i-1)*dimPart+1):(i*dimPart), ((j-1)*dimPart+1):(j*dimPart) ] = Kpart
      if(i != j) K[ ((j-1)*dimPart+1):(j*dimPart), ((i-1)*dimPart+1):(i*dimPart)] <- t(K[ ((i-1)*dimPart+1):(i*dimPart), ((j-1)*dimPart+1):(j*dimPart)])
    }
  }

  
  save(K,file="K.RData")
  proc.time() - ptm
  
  
  # eigen.R
  # Reconstructing POD modes and coefficients
  # Doing eigendecomposition
  load("K.RData")
  eigdat = eigs_sym(K, (numModes))
  save(eigdat, file = "EigenSystem.RData")
  
  # Reconstructing POD modes and coefficients
  ptm <- proc.time()
  FTlist <- vector("list",numPart)
  tFTlist <- vector("list",numPart)
  for (i in 1:numPart){
    FTlist[[i]] <- attach.big.matrix(paste0("FTtmp",i,".dsc"))
    tFTlist[[i]] <- big.t(FTlist[[i]])
  }
  proc.time() - ptm
  
  load("EigenSystem.RData")
  eiglist <- vector("list",numPart)
  for (i in 1:numPart){
    eiglist[[i]] <- as.big.matrix(eigdat$vector[((i-1)*dimPart+1):(i*dimPart),])
  }
  
  # Computing POD modes and coefficients
  ptm <- proc.time()
  
  #POD modes
  PsiTlist <- vector("list",numPart)
  for (i in 1:numPart){
    PsiTlist[[i]] <- bigInnerProd2(tFTlist[[i]]@address,eiglist[[i]]@address)
  }
  PsiT <- matrix(0,nrow=nrow(FTlist[[1]]),ncol=numModes)
  for (i in 1:numPart){
    PsiT = PsiT + PsiTlist[[i]]
  }
  PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
  #POD coefficients
  PsiTbig <- as.big.matrix(PsiT, backingfile = "PsiT.bck", backingpath = getwd(), descriptorfile = "PsiT.dsc")
  BT <- matrix(0,nrow=0,ncol=numModes)
  for (i in 1:numPart){
    BT = rbind( BT, bigInnerProd2(FTlist[[i]]@address,PsiTbig@address) )
  }
  
  proc.time() - ptm
  
  #Save POD information
  PsiTsvg = PsiT
  BTsvg = BT
  eigval = eigdat$values
  eigvec = eigdat$vectors
  save(eigval,eigvec,PsiTsvg,BTsvg, file = paste0("PODdat.RData"))
  
}



#################################################################################
#Test reconstruction
#################################################################################

#Save reconstruction
cas = 1
#   FT <- attach.big.matrix(paste0("FT_",cas,".dsc"))
#   slc = ( (999*(cas-1)+1):(999*cas) ) - (999*3)
for (k in 1:length(indvec)){
  setwd(paste0(filepath,"bigmemory_bk/", indvec[k]))
  numCol <- 1000 #numCol colors in legend
  colors <- matlab.like(numCol+1)
  
  #Load common coordinate indices and data
  load("../../CommonCoord.RData")
  load("FTmean.vector.RData")
  load("PODdat.RData")
  load("na.index.RData")
  #na.index <- c(192, 385, 578, 771, 964)
  
  #Assign range for plot
  if (indvec[k]==4){
    Range = c(-50, 50)
  }else if(indvec[k]==5){
    Range = c(-30, 30)
  }else if(indvec[k]==6){
    Range = c(-30, 0)
  }else if(indvec[k]==7){
    Range = c(-400,400)
  }else if(indvec[k]==8){
    Range = c(100, 350)
  }else if(indvec[k]==9){
    Range = c(0, 1000)
  }
  
  for (j in cas:cas){
    Coordinate = rbind(CommonCoord[[j]][[1]],CommonCoord[[j]][[2]],CommonCoord[[j]][[3]],CommonCoord[[j]][[4]])
    
    #Plot reconstruction
    for (i in 1:check.tm.num){
      print(paste0("case ",j,", timestep ",i))  
      zcolor <- colors[ pmax( (FTmean.vector + PsiTsvg[,1:cutoff.numModes]%*%matrix(BTsvg[(j-1)*TT+i,1:cutoff.numModes],ncol=1) - Range[1])/(Range[2]-Range[1])*numCol, 0) + 1 ]
      #        zcolor <- colors[ pmax( (FTmean.vector + FTtmp[,i] - min(FT[,i]))/diff(range(FT[,i]))*numCol, 0) + 1] 
      png(paste0("../../check/recon-",indvec[k],"-",i,"_",cutoff.numModes,".png"), width = 600, height = 400)
      plot(Coordinate[-na.index,1],Coordinate[-na.index,2],col= zcolor,pch=15,cex=0.5,
           xlab = 'x (m)', ylab = 'y (m)', main = paste0("Reconstructed at T = ", i),
           xlim = c(0,0.04), ylim = c(0,0.008) ) 
      colorlegend(col=colors,zlim=Range,left=T)
      dev.off()
    }
    
    #Plot actual flow
    FT.val <- attach.big.matrix(paste0("FT_",cas,".dsc"))
    for (i in 1:check.tm.num){
      print(paste0("case ",j,", timestep ",i))  
      zcolor <- colors[ pmax( (FT.val[,i] - Range[1])/(Range[2]-Range[1])*numCol, 0) + 1 ]
      #        zcolor <- colors[ pmax( (FTmean.vector + FTtmp[,i] - min(FT[,i]))/diff(range(FT[,i]))*numCol, 0) + 1] 
      png(paste0("../../check/true-",indvec[k],"-",i,".png"), width = 600, height = 400)
      plot(Coordinate[,1],Coordinate[,2],col= zcolor,pch=15,cex=0.5,
           xlab = 'x (m)', ylab = 'y (m)', main = paste0("True at T = ", i),
           xlim = c(0,0.04), ylim = c(0,0.008) ) 
      colorlegend(col=colors,zlim=Range,left=T)
      dev.off()
    }
    
  }
}

#Coef X POD mode and output dat data
load("../../CommonCoord.RData")
for(ind in indvec){
  setwd(paste0(filepath,"bigmemory_bk/", ind))
  numCol <- 1000 #numCol colors in legend
  colors <- matlab.like(numCol+1)
  
  load("FTmean.vector.RData")
  load("PODdat.RData")
  load("na.index.RData")
  #na.index <- c(192, 385, 578, 771, 964)
  
  for (j in POD.case){
    Coordinate = rbind(CommonCoord[[j]][[1]],CommonCoord[[j]][[2]],CommonCoord[[j]][[3]],CommonCoord[[j]][[4]])
    for(k in 1:11){
      for (i in 1:1){
        print(paste0("mode", k, " case ",j,", timestep ",i))  
        if(k == 1){
          CoefTimePOD <- FTmean.vector
        }else{
          CoefTimePOD <- PsiTsvg[,k-1] * matrix(BTsvg[(j-1)*TT+i,k-1],ncol=1)
        }
        
        Range <- range(CoefTimePOD)
        zcolor <- colors[ pmax( (CoefTimePOD - Range[1])/(Range[2]-Range[1])*numCol, 0) + 1 ]
        #        zcolor <- colors[ pmax( (FTmean.vector + FTtmp[,i] - min(FT[,i]))/diff(range(FT[,i]))*numCol, 0) + 1] 
        png(paste0("../../check/mode", k-1, "-", ind,"-",i,".png"), width = 600, height = 400)
        plot(Coordinate[-na.index,1],Coordinate[-na.index,2],col= zcolor,pch=15,cex=0.5,
             xlab = 'x (m)', ylab = 'y (m)', main = paste0("POD = ", k-1, " times coef at T = ", i, " of ind = ", ind)) 
        colorlegend(col=colors,zlim=Range,left=T)
        dev.off()
        OutputDat <- cbind(Coordinate[-na.index,], CoefTimePOD)
        colnames(OutputDat) <- c("x", "y", "Mode")
        write.table(OutputDat, file=paste0("../../POD_result/dat/Mode", k - 1, "-Ind", ind,"-Case",i,".dat"), sep = '\t', row.names = F)
      }
    }
  }
}
