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


# #Test
# # set up big.matrix
# # nrows <- 90000
# # setwd("/tmp")
# bkFile <- "bigmat.bk"
# descFile <- "bigmatk.desc"
# bkFile2 <- "bigmat2.bk"
# descFile2 <- "bigmatk2.desc"
# numcol = 3
# nrows = 1000
# bigmat <- filebacked.big.matrix(nrow=nrows, ncol=numcol, type="double",
#                                 backingfile=bkFile, backingpath=".",
#                                 descriptorfile=descFile,
#                                 dimnames=c(NULL,NULL))
# bigmat2 <- filebacked.big.matrix(nrow=nrows, ncol=numcol*2, type="double",
#                                 backingfile=bkFile2, backingpath=".",
#                                 descriptorfile=descFile2,
#                                 dimnames=c(NULL,NULL))
# set.seed(123)
# for (i in 1:numcol) bigmat[,i] = rnorm(nrows)*i
# for (i in 1:(2*numcol)) bigmat2[,i] = rnorm(nrows)*i
# 
# smallmat <- matrix(rnorm(numcol*5),nrow=numcol,ncol=5)
# 
# ptm <- proc.time()
# K = bigRecon(bigmat@address,smallmat)
# proc.time() - ptm
# 
# #Check
# mat = as.matrix(bigmat)
# t(mat)%*%mat

############################################################################
#Set working directory where data is
ind = 6
setwd(paste0("/home/proj/jeffwu/isye-ae/bigmemory_bk/", ind))

#Loading libraries
FTtmp1 <- attach.big.matrix(paste0("FTtmp1.dsc"))
FTtmp2 <- attach.big.matrix(paste0("FTtmp2.dsc"))

# #Only for case 25
# cas = 26
# slc = ( (999*(cas-1)+1):(999*cas) ) - 14985
# # #check
# # ( slc + (999*30)/2 )%%999
# FT25slice <- as.big.matrix(FTtmp2[,slc], backingfile = "FT25.bck", backingpath = getwd(), descriptorfile = "FT25.dsc")
# ptm <- proc.time()
# K25slice1 = bigInnerProd2(FT25slice@address,FTtmp1@address)
# K25slice2 = bigInnerProd2(FT25slice@address,FTtmp2@address)
# proc.time() - ptm
# save(K25slice1,file="K25slice1.RData")
# save(K25slice2,file="K25slice2.RData")
# 
# K25slice = cbind(K25slice1,K25slice2)
# 
# #Replace each slice
# load("Old stuff/K.RData")
# K[slc,] = K25slice
# K[,slc] = t(K25slice)
# save(K,file="K.RData")
#
# #Compute K11
# ptm <- proc.time()
# K11 = bigInnerProd(FTtmp1@address)
# K11 = symmetrize(K11)
# proc.time() - ptm
# save(K11,file="K11.RData")
# 
# #Compute K12
# ptm <- proc.time()
# K12 = bigInnerProd2(FTtmp1@address,FTtmp2@address)
# proc.time() - ptm
# save(K12,file="K12.RData")
# 
# #Compute K22
# ptm <- proc.time()
# K22 = bigInnerProd(FTtmp2@address)
# K22 = symmetrize(K22)
# proc.time() - ptm
# save(K22,file="K22.RData")
# 
# # Constructing and saving K
# load("K11.RData")
# load("K12.RData")
# load("K22.RData")
# K = rbind(cbind(K11,K12),cbind(t(K12),K22))
# save(K,file="K.RData")

# Doing eigendecomposition
numModes = 300
# load("K.RData")
eigdat = eigs_sym(K, (numModes))
save(eigdat, file = "EigenSystem.RData")

# Reconstructing POD modes and coefficients
ptm <- proc.time()
tFTtmp1 <- big.t(FTtmp1,name="tFTtmp1")
tFTtmp2 <- big.t(FTtmp2,name="tFTtmp2")

proc.time() - ptm

load("EigenSystem.RData")
eig1 <- as.big.matrix(eigdat$vector[1:nrow(tFTtmp1),], backingfile = "eig1.bck", backingpath = getwd(), descriptorfile = "eig1.dsc")
eig2 <- as.big.matrix(eigdat$vector[(nrow(tFTtmp1)+1):nrow(eigdat$vector),], backingfile = "eig2.bck", backingpath = getwd(), descriptorfile = "eig2.dsc")

ptm <- proc.time()
#POD modes
PsiT1 = bigInnerProd2(tFTtmp1@address,eig1@address)
PsiT2 = bigInnerProd2(tFTtmp2@address,eig2@address)
PsiT = PsiT1 + PsiT2
PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
#POD coefficients
PsiTbig <- as.big.matrix(PsiT, backingfile = "PsiT.bck", backingpath = getwd(), descriptorfile = "PsiT.dsc")
BT1 <- bigInnerProd2(FTtmp1@address,PsiTbig@address)
BT2 <- bigInnerProd2(FTtmp2@address,PsiTbig@address)
BT <- rbind(BT1, BT2)
proc.time() - ptm

#Save POD information
PsiTsvg = PsiT
BTsvg = BT
eigval = eigdat$values
eigvec = eigdat$vectors
save(eigval,eigvec,PsiTsvg,BTsvg, file = paste0("PODdat.RData"))

#Test reconstruction
ind = 4
setwd(paste0("/home/proj/jeffwu/isye-ae/bigmemory_bk/", ind))
numCol <- 1000 #numCol colors in legend
colors <- matlab.like(numCol+1)
if (ind==4){
  Range = c(-101, 69)
}else if(ind==5){
  Range = c(-46, 56)
}else if(ind==6){
  Range = c(-65, 12)
}else if(ind==7){
  Range = c(-5e5,1180000)
}else if(ind==8){
  Range = c(120, 310)
}else if(ind==9){
  Range = c(125, 1010)
}

#Load common coordinate indices and data
load("../../CommonCoord.RData")
load("FTmean.vector.RData")
load("PODdat.RData")
na.index <- c(385,578)
# nrow(CommonCoord[[1]][[1]]) + nrow(CommonCoord[[1]][[2]]) + nrow(CommonCoord[[1]][[3]]) + nrow(CommonCoord[[1]][[4]])

#Save reconstruction
cas = 25
for (j in cas:cas){
  Coordinate = rbind(CommonCoord[[j]][[1]],CommonCoord[[j]][[2]],CommonCoord[[j]][[3]],CommonCoord[[j]][[4]])
  for (i in 1:50){
    print(paste0("case ",j,", timestep ",i))  
#     zcolor <- pmin( pmax( (FTmean.vector + PsiTsvg%*%matrix(BTsvg[(j-1)*999+i,],ncol=1) - Range[1])/(Range[2]-Range[1])*numCol , 0 ), numCol ) + 1
    zcolor <- pmin( pmax( (FTmean.vector + FTtmp2[,slc[i]] - Range[1])/(Range[2]-Range[1])*numCol , 0 ), numCol ) + 1
    png(paste0("../check/",ind,"-",i,".png"), width = 600, height = 400)
    plot(Coordinate[-na.index,1],Coordinate[-na.index,2],col=colors[zcolor],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Reconstructed at T = ", i),
         xlim = c(0,0.10), ylim = c(0,0.013) ) 
    colorlegend(col=colors,zlim=Range,left=T)
    dev.off()
  }
}
