##########################################################
#         *** Generate Big FT ***                        #
#    MADE BY CHIH-LI and SIMON ON 12/6/2015              #
#    DESCRIPTION: before doing inner product,            #
#                 obtain the Centered FT                 #
#    NOTE: PLEASE READ THE MANUSCRIPT                    #
#    Reports & Presentations/JASA paper/paper.pdf        # 
##########################################################

#Step 0: input
setwd("RawData/")                    # set working directory
num.case <- 30
available_file.tmp <- 0:998
TT <- length(available_file.tmp)
split.num <- 10
indvec <- c(4:9)
if(num.case %% split.num != 0) stop("num.case %% split.num != 0")

#Loading libraries
# install.packages("bigmemory",repos="http://cran.us.r-project.org")
library(foreach)
library(doSNOW)
library(parallel)
library(bigmemory)
library(bigalgebra)
library(bigpca)

ptm <- proc.time()
load("../CommonCoord.RData")
grid.num <- nrow(CommonCoord[[1]][[1]]) + nrow(CommonCoord[[1]][[2]]) +
  nrow(CommonCoord[[1]][[3]]) + nrow(CommonCoord[[1]][[4]])

#cl <- makeCluster(detectCores(),type="SOCK")
#registerDoSNOW(cl)
for(ind in indvec){
  #foreach(ind = c(4:9, 13), .packages = "bigmemory") %dopar% {
  cat("ind =", ind, "\n")
  FT <- filebacked.big.matrix(grid.num, TT * num.case,                      
                              backingfile = "FT_All.bck",  type='double',
                              backingpath = paste0("../bigmemory_bk/", ind), 
                              descriptorfile = "FT_All.dsc")
  
  for(i in 1:num.case){
    cat("Making FT_All step",i,"\n")
    FT_part.desc <- attach.big.matrix(paste0("../bigmemory_bk/", ind, "/", "FT_",i,".dsc"))
    column.index <- (1 + (i-1) * TT) : (TT * i) 
    FT[,column.index] = FT_part.desc[]
  }
  
  
  FT <- attach.big.matrix(paste0("../bigmemory_bk/", ind, "/", "FT_All.dsc"))
  
  ##### Compute mean at each grid
  na.index <- c(192, 385, 578, 771, 964)
  FTmean.vector <- rep(NA, length(available_file.tmp))
  for(i in 1:length(available_file.tmp)){
    FTmean.vector[i] <- mean(FT[-na.index,i + length(available_file.tmp) * 0:29])
  }
  
  #FTmean.vector <- bmcapply(FT, 1, mean, dir = paste0("../bigmemory_bk/", ind))
  #load(paste0("../bigmemory_bk/", ind, "/FTmean.vector.RData"))
  #na.index.fg <- is.na(FTmean.vector)
  #na.index <- which(na.index.fg)
  #na.index <- c(385,578)
  #na.index <- c(192, 385, 578, 771, 964)
  #print(which(na.index.fg))
  #save(na.index, file = paste0("../bigmemory_bk/", ind, "/na.index.RData"))
  #FTmean.vector <- FTmean.vector[!na.index.fg]
  save(FTmean.vector, file = paste0("../AE test for ICA/", ind, "/FTmean.vector.RData"))
  
  FT.centered <- filebacked.big.matrix(nrow(FT) - length(na.index), TT * num.case,                      
                                       backingfile = "FT.centered.bck", 
                                       backingpath = paste0("../AE test for ICA/", ind), 
                                       descriptorfile = "FT.centered.dsc")
  FTmean.vector.mx <- matrix(rep(FTmean.vector, nrow(FT) - length(na.index)),
                             nrow = nrow(FT) - length(na.index), byrow = TRUE)  
  for(ii in 1:num.case){
    cat("Centering FT.centered step", ii,"\n")
    col.index <- (1 + (ii - 1) * TT) : (ii * TT)
    FT.centered[,col.index] <- FT[-na.index, col.index] - FTmean.vector.mx
  }
  
  for(ii in 1:split.num){
    col.index <- (1 + (ii - 1) * (num.case/split.num) * TT) : (ii * (num.case/split.num) * TT)
    FTtmp <- filebacked.big.matrix(nrow(FT) - length(na.index), length(col.index),                      
                                   backingfile = paste0("FTtmp",ii,".bck"),  type='double',
                                   backingpath = paste0("../AE test for ICA/", ind), 
                                   descriptorfile = paste0("FTtmp",ii,".dsc"))
    
    for(kk in 1:(num.case/split.num)){
      cat("Centering FTtmp", ii, "step", kk,"\n")
      col.tmp.index <- col.index[(1 + (kk - 1) * TT) : (kk * TT)]
      FTtmp[,(1 + (kk - 1) * TT) : (kk * TT)] <- FT.centered[,col.tmp.index]
    }
  }
}
#stopCluster(cl)
print(proc.time() - ptm)
