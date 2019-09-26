#Set working directory where data is
setwd("/home/proj/jeffwu/isye-ae/RawData/")

#Loading libraries
library(foreach)
library(doSNOW)
library(parallel)

#Choose flow property
ind = 7 #Flow temperature

#Number of modes desired
POD_mode <- 500

#Timesteps available
# available_file.tmp <- (1:999)[-c(1,853)] 
numTime <- 200
available_file.tmp <- 3:(3+numTime)
available_file.tmp <- setdiff(available_file.tmp,seq(from=13,by=11,to=1000)) #Remove jumps
TT = length(available_file.tmp)
wonkyInd = c(1,3,13,14,15,16,19,20:30)

#Load common coordinate indices
load("../CommonCoord.MinDistIndex.RData")
load("../CutCoordInd.RData")

#Set up clusters
cl <- makeCluster(detectCores(),type="SOCK")
registerDoSNOW(cl)

#Compute FT
ptm <- proc.time()
FT = foreach(j = 1:30, .combine=cbind) %dopar% {
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
  
  #Preallocate
  allSnaps = matrix(0.0,nrow=len,ncol=TT)
  for (i in 1:TT){
    #     print(i)
    tm = available_file.tmp[i]
    #Read snapshot for case j at timestep tm
    week = ceiling(j/6)
    if (j %in% wonkyInd){
      snapT <- read.table(paste0(j,"/Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
    }
    else{
      snapT <- read.table(paste0(j,"/Ins_All_Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
    }
    
    #Extract data at common coordinates
    allSnaps[,i] = c(snapT[cutcoord1[coord1]],snapT[cutcoord2[coord2]],snapT[cutcoord3[coord3]],snapT[cutcoord4[coord4]])
  }
#     proc.time() - ptm
  
  return(allSnaps)
}

Fmean = rowMeans(FT)

save(FT, file=paste0("../FT_",numTime,"t_",ind,"f.RData"))
save(Fmean, file=paste0("../Fmean_",numTime,"t_",ind,"f.RData"))
#
proc.time() - ptm

# Timing: numTime = 50: 284s (should scale linearly)

# stopCluster(cl)

#WARNING: THIS TAKES VERY LONG!
#Center FT and compute Gram matrix
ptm <- proc.time()

FT.centered <- t(scale(t(FT),center=T,scale=F))
K <- crossprod(FT.centered)
save(K,file=paste0("../K_",numTime,"t_",ind,"f.RData"))
#
proc.time() - ptm

# Timing: ???

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
numTime = 200
ind = 7

#Set working directory where data is
setwd("/home/proj/jeffwu/isye-ae/RawData/")
load(paste0("../rARPACK_output/K_",numTime,"t_",ind,"f.RData"))
load(paste0("../rARPACK_output/FT_",numTime,"t_",ind,"f.RData"))
load(paste0("../rARPACK_output/Fmean_",numTime,"t_",ind,"f.RData"))

eigdat = eigs_sym(K, (numModes))
eigval = eigdat$values
eigvec = eigdat$vectors

save(eigval,eigvec, file = paste0("../rARPACK_output/eigen_",numTime,"t_",numModes,"m_", ind, "f.RData"))

#Switch back to castle for computations
setwd("/home/proj/jeffwu/isye-ae/RawData/")
numModes = 500
numTime = 200
ind = 7

load(paste0("../rARPACK_output/eigen_",numTime,"t_",numModes,"m_", ind, "f.RData"))
load(paste0("../rARPACK_output/FT_",numTime,"t_",ind,"f.RData"))
load(paste0("../rARPACK_output/Fmean_",numTime,"t_",ind,"f.RData"))

FT.centered <- (FT - Fmean)
PsiT <- FT.centered %*% eigvec #POD modes
PsiT <- t(t(PsiT) / apply(PsiT, 2, FUN = function(x) sqrt(sum(x^2)))) # Normalize
BT <- t(PsiT) %*% FT.centered #POD coefficients

PsiTsvg = PsiT
BTsvg = BT

save(eigval,eigvec,PsiTsvg,BTsvg, file = paste0("../rARPACK_output/POD_",numTime,"t_",numModes,"m_", ind, "f.RData"))


