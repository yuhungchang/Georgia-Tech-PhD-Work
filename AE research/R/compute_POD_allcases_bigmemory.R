#Set working directory where data is
setwd("/home/proj/jeffwu/isye-ae/RawData/")

#Loading libraries
# install.packages("bigmemory",repos="http://cran.us.r-project.org")
library(foreach)
library(doSNOW)
library(parallel)
library(bigmemory)
library(bigalgebra)
library(bigpca)

#Choose flow property
ind = c(4:9) 
cat("ind = ", ind, "...\n")
#if(!dir.exists(paste0("../bigmemory_bk/", ind))) dir.create(paste0("../bigmemory_bk/", ind))
  
#Number of modes desired
POD_mode <- 500

#Timesteps available
# available_file.tmp <- (1:999)[-c(1,853)] 
#numTime <- 600
#available_file.tmp <- 3:(3+numTime)
#available_file.tmp <- setdiff(available_file.tmp,seq(from=13,by=11,to=1000)) #Remove jumps
available_file.tmp <- 0:998
TT = length(available_file.tmp)
#wonkyInd = c(1,3,13,14,15,16,19,20:30)
wonkyInd = 1:30
rescaleInd_by_1000 = 7
rescaleInd_by_million = 13

#Load common coordinate indices
load("../CommonCoord.MinDistIndex.RData")
load("../CutCoordInd.RData")

#Set up clusters
cl <- makeCluster(detectCores(),type="SOCK")
registerDoSNOW(cl)

#Compute FT
ptm <- proc.time()
#FT = as.big.matrix(foreach(j = 1:30, .combine=cbind) %dopar% {
for(j in 1:30){
  print(j)
  #j indexes case
  # ptm <- proc.time()  
  coord1 = CommonCoord.MinDistIndex[[j]][[1]]
  coord2 = CommonCoord.MinDistIndex[[j]][[2]]
  coord3 = CommonCoord.MinDistIndex[[j]][[3]]
  coord4 = CommonCoord.MinDistIndex[[j]][[4]]
  
  cutcoord1 = CutCoordInd[[j]][[1]]
  cutcoord2 = CutCoordInd[[j]][[2]]
  cutcoord3 = CutCoordInd[[j]][[3]]
  cutcoord4 = CutCoordInd[[j]][[4]]
  
  len = length(coord1) + length(coord2) + length(coord3) + length(coord4)
  week = ceiling(j/6)
  if(j == 14) week = 1 
  
  #Preallocate
  #allSnaps = matrix(0.0,nrow=len,ncol=TT)
  #for (i in 1:TT){
  FT = as.big.matrix(foreach(i = available_file.tmp, .combine=cbind) %dopar% {
    #     print(i)
    tm = i
    #Read snapshot for case j at timestep tm
    if (j %in% wonkyInd){
      snapT <- read.table(paste0(j,"/Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
      #snapT <- read.big.matrix(paste0(j,"/Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', type = "double", backingfile = "test.bck", backingpath = getwd(), descriptorfile = "test.dsc")[,ind]
    }
    else{
      snapT <- read.table(paste0(j,"/Ins_All_Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
      #snapT <- read.big.matrix(paste0(j,"/Ins_All_Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', type = "double", backingfile = "test.bck", backingpath = getwd(), descriptorfile = "test.dsc")[,ind]
    }
    
    #Extract data at common coordinates
    #allSnaps[,i] = c(snapT[cutcoord1[coord1]],snapT[cutcoord2[coord2]],snapT[cutcoord3[coord3]],snapT[cutcoord4[coord4]])
    snapT <- snapT[c(cutcoord1[coord1],cutcoord2[coord2],
                     cutcoord3[coord3],cutcoord4[coord4]), ]
    snapT[,ind == rescaleInd_by_1000] <- snapT[,ind == rescaleInd_by_1000]/1000
    snapT[,ind == rescaleInd_by_million] <- snapT[,ind == rescaleInd_by_million]/1000000
    
    return(snapT)
  }, backingfile = paste0("FT_tmp.bck"), backingpath = paste0("../bigmemory_bk"), descriptorfile = paste0("FT_tmp.dsc"))
  
  for(kk in 1:length(ind)){
    column.index <- kk + 0:(length(available_file.tmp)-1) * length(ind)
    FT.tmp <- as.big.matrix(FT[, column.index], backingfile = paste0("FT_",j,".bck"), backingpath = paste0("../bigmemory_bk/", ind[kk]), descriptorfile = paste0("FT_",j,".dsc"))
  }
  #     proc.time() - ptm
  
#   if(j == 1){
#     FT = allSnaps
#   }else{
#     FT = as.big.matrix(cbind(FT[], allSnaps[]), backingfile = paste0("FT_",ind,".bck"), backingpath = getwd(), descriptorfile = paste0("FT_",ind,".dsc"))
#   }
}

Fmean = bmcapply(FT, 1, mean,n.cores=8)
##

write.big.matrix(FT, filename=paste0("test/rARPACK_output/FT_",numTime,"t_",ind,"f.RData"))
save(Fmean, file=paste0("test/rARPACK_output/Fmean_",numTime,"t_",ind,"f.RData"))

proc.time() - ptm
stopCluster(cl)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Timing:
#  TT = 11  -> 67s
#  TT = 101 -> 723s
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#WARNING: THIS TAKES VERY LONG!
#Center FT and compute Gram matrix
ptm <- proc.time()
FT.t.centered <- as.big.matrix(bmcapply(FT, 1,scale, center=T, scale=F, n.cores=8), backingfile = "FT.centered.t.bck", backingpath = getwd(), descriptorfile = "FT.centered.t.dsc")
FT.centered <- big.t(FT.t.centered, verbose=TRUE)
#K <- crossprod(FT.centered)
K <- FT.t.centered %*% FT.centered
K <- K[]
save(K,file=paste0("test/rARPACK_output/K_",numTime,".RData"))

proc.time() - ptm

#... now move to wren to use packages

######################################################################
#   #Nystrom sampling
######################################################################
#   ptm <- proc.time()
#   spac = floor(TT/qq)
#   KTq <- K[,ceiling(spac*1:qq-spac/2)]
#   Kq <- K[ceiling(spac*1:qq-spac/2),ceiling(spac*1:qq-spac/2)]
#   
#   EigSystem <- eigen(Kq)
#   Uq <- EigSystem$vector
#   lambda_q <- EigSystem$value
#   lambda_T <- TT/qq * EigSystem$value #Estimated eigenvalues of Q
#   
#   # UT <- sqrt(qq/TT) / lambda_q * KTq %*% Uq #Estimated eigenvectors of Q
#   # UT <- sqrt(qq/TT) / lambda_q * K %*% Uq #Estimated eigenvectors of Q
#   UT <- sqrt(qq/TT) * KTq %*% Uq[,1:POD_mode] %*% diag(1/lambda_q[1:POD_mode])#Estimated eigenvectors of Q
#   PsiT <- FT.centered %*% UT #POD modes
#   PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
#   
#   # BT <- diag(lambda_T) %*% t(UT) #POD coefficients
#   BT <- t(PsiT) %*% FT.centered #POD coefficients

######################################################################
#   #rARPACK package
######################################################################
#Now compute the eigenvectors (qq modes)
library(rARPACK)
numModes = 500
numTimes = 200

#Set working directory where data is
setwd("/home/proj/jeffwu/isye-ae/RawData/")
load(paste0("test/rARPACK_output/K_",numTimes,".RData"))
load(paste0("test/rARPACK_output/FT_",numTimes,".RData"))
load(paste0("test/rARPACK_output/Fmean_",numTimes,".RData"))

eigdat = eigs_sym(K, (numModes))
eigval = eigdat$values
eigvec = eigdat$vectors

FT.centered <- (FT - Fmean)
PsiT <- FT.centered %*% eigvec #POD modes
PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
BT <- t(PsiT) %*% FT.centered #POD coefficients

PsiTsvg = PsiT
BTsvg = BT

save(eigval,eigvec,PsiTsvg,BTsvg, file = paste0("test/rARPACK_output/POD_",numTime,"_",numModes,"m.RData"))


