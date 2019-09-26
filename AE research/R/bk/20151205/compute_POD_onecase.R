#Set working directory where data is
setwd("../../proj/jeffwu/isye-ae/RawData/")

#Loading libraries
library(foreach)
library(doSNOW)
library(parallel)
ptm <- proc.time()

#Set up clusters
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

#Choose flow property
ind = 8 #Flow temperature

#Subsampling size
TT = 1000 #Dimensions of Q
#qq = 976 #Dimension of subsample
POD_mode <- 10

#Nystrom
available_file.tmp <- (1:999)[-c(1,853)] 

FT = foreach(i = 1:TT, .combine=cbind) %dopar% {
  if(any(i == available_file.tmp)){
    #Read snapshot at i
    snapT <- read.table(paste("1/Week1_1_",i,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
    return(snapT)
  }else {
    return(0)
  }
}
stopCluster(cl)

TT <- length(available_file.tmp)
FT <- FT[,available_file.tmp]
FT.centered <- t(scale(t(FT),center=T,scale=F))

K = crossprod(FT.centered)
#dotMean = (1/TT^2) * (sum(K)) #Dot product of mean snapshot
#rowMean = rowMeans(K) #Matrix of row means
#colMean = rowMean #Matrix of column means

#K = t(t(K - rowMean) - colMean) + dotMean #Normalize
for(qq in c(50, 100, 200, 500, 997)){
  cat("qq =", qq, "\n")  
  ptm <- proc.time()
  spac = floor(TT/qq)
  KTq <- K[,ceiling(spac*1:qq-spac/2)]
  Kq <- K[ceiling(spac*1:qq-spac/2),ceiling(spac*1:qq-spac/2)]
  
  EigSystem <- eigen(Kq)
  Uq <- EigSystem$vector
  lambda_q <- EigSystem$value
  lambda_T <- TT/qq * EigSystem$value #Estimated eigenvalues of Q
  
  # UT <- sqrt(qq/TT) / lambda_q * KTq %*% Uq #Estimated eigenvectors of Q
  # UT <- sqrt(qq/TT) / lambda_q * K %*% Uq #Estimated eigenvectors of Q
  UT <- sqrt(qq/TT) * KTq %*% Uq[,1:POD_mode] %*% diag(1/lambda_q[1:POD_mode])#Estimated eigenvectors of Q
  PsiT <- FT.centered %*% UT #POD modes
  PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
  
  # BT <- diag(lambda_T) %*% t(UT) #POD coefficients
  BT <- t(PsiT) %*% FT.centered #POD coefficients
  
  PsiTsvg = PsiT
  BTsvg = BT
  save(PsiTsvg, file = paste0("../", POD_mode,"mode/Psi_", qq,".RData"))
  save(BTsvg, file = paste0("../", POD_mode,"mode/BT_", qq,".RData"))
  
  print(proc.time() - ptm)
}


## Timing O(q^3) 
#  Load  : 117s
#  q = 10: 6.2s
#  q = 20: 7.9s
#  q = 50: 15.2s
#  q = 200: 45.4s
#  q = 976: 201.2s