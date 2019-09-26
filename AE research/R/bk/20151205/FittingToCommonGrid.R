#library(fBasics)
#ans = linearInterp(data.df[,1], data.df[,2], data.df[,3],  gridPoints = nrow(CommonCoord.1),
#                  xo = CommonCoord.1[,1], yo = CommonCoord.1[,2])
setwd("../../proj/jeffwu/isye-ae/")
load("CommonCoord.RData")
load("CutCoordInd.RData")
library(foreach)
library(doSNOW)
library(parallel)
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

CommonCoord.MinDistIndex <- vector("list", 30)
for(i in 1:30){
  CommonCoord.MinDistIndex[[i]] <- vector("list", 4)
}
#wonkyInd = c(1,3,13,14,15,16,19,20:30)
wonkyInd = 1:30
wrongInd = 25

for (i in 1:30){
  print(i)
  
  #Rename files
  week = ceiling(i/6)
  if (i == 14) week = 1
#   if (i %in% wrongInd){
#     RealCoord <- read.table(paste0("RawData/", i,"/Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,1:2]
#     RealCoord[,1] <- RealCoord[,1] - min(RealCoord[,1])
#     RealCoord[,2] <- RealCoord[,2] - min(RealCoord[,2])
#   }else 
  if (i %in% wonkyInd){
    RealCoord <- read.table(paste0("RawData/", i,"/Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,1:2]
  }
  else{
    RealCoord <- read.table(paste0("RawData/",i,"/Ins_All_Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,1:2]
  }
  
  CommonCoord.i <- CommonCoord[[i]]
  CutCoordInd.i <- CutCoordInd[[i]]
  
  ptm <- proc.time()
  for(j in 1:4){
      CommonCoord.i.j <- CommonCoord.i[[j]]
      CutCoordInd.i.j <- CutCoordInd.i[[j]]
      #coords <- t(RealCoord[CutCoordInd.i.j,1:2])
      coords <- t(RealCoord[,1:2])
      
      CommonCoord.MinDistIndex[[i]][[j]] <- foreach(k = 1:nrow(CommonCoord.i.j), .combine = c) %dopar% {
#           x <- as.matrix(CommonCoord.i.j[k,])[,1]
#           SumSquare <- apply(sweep(RealCoord[CutCoordInd.i.j,1:2], 2, x, FUN = "-"), 1, FUN = function(y) sum(y^2))
#           return(which.min(SumSquare))
        which.min(colSums((coords - CommonCoord.i.j[k,])^2))    
      }
    
    names(CommonCoord.MinDistIndex[[i]][[j]]) <- NULL
  }
  proc.time() - ptm
}

save(CommonCoord.MinDistIndex, file = "CommonCoord.MinDistIndex.RData")

stopCluster(cl)

# 
# 
# ## plot
# library(colorRamps)
# library(shape)
# load("CommonCoord.MinDistIndex.RData")
# Range = c(120, 309.8351) # Same range with true FT
# numCol <- 1000
# colors <- matlab.like(numCol+1)
# wonkyInd = c(1,3,13,14,15,16,19,20:30)
# 
# for (i in 1:30){
#   print(i)
#   
#   #Rename files
#   week = ceiling(i/6)
#   if (i %in% wonkyInd){
#     RealCoord <- read.table(paste0("RawData/", i,"/Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,c(1:2,8)]
#   }
#   else{
#     RealCoord <- read.table(paste0("RawData/",i,"/Ins_All_Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,c(1:2,8)]
#   }
#   png(paste0("checkCommonGrid.True.", i , ".png"), width = 600, height = 400)
#   zcolor <- colors[(RealCoord[,3] - min(Range))/diff(range(Range))*numCol + 1]
#   plot(RealCoord[,1], RealCoord[,2],col=zcolor,pch=15,cex=1.2,
#        xlab = 'x (m)', ylab = 'y (m)') 
#   colorlegend(col=colors,zlim=range(Range),left=T)
#   dev.off()
#   
#   
#   for(j in 1:4){
#     CommCoord.AvgResponse <- RealCoord[CutCoordInd[[i]][[j]],3][CommonCoord.MinDistIndex[[i]][[j]]]
#     zcolor <- colors[(CommCoord.AvgResponse - min(Range))/diff(range(Range))*numCol + 1]
#     png(paste0("checkCommonGrid.", i, ".", j, ".png"), width = 600, height = 400)
#     
#     plot(CommonCoord[[i]][[j]][,1], CommonCoord[[i]][[j]][,2],col=zcolor,pch=15,cex=1.2,
#          xlab = 'x (m)', ylab = 'y (m)') 
#     colorlegend(col=colors,zlim=range(Range),left=T)
#     dev.off()
#   }
# }  
# 
# plot(CommonCoord[[1]][[1]][,1], CommonCoord[[1]][[1]][,2],pch=15,cex=1.2,
#      xlab = 'x (m)', ylab = 'y (m)') 
# dev.off()
# 
# ###### CANNOT WORK ######
# d <- 10
# CommonCoord.DistanceAndIndex <- foreach(i = 1:nrow(CommonCoord.1.1), .combine = rbind) %dopar% {
#   x <- as.matrix(CommonCoord.1.1[i,])[,1]
#   RootSumSquare <- apply(sweep(data.df[,1:2], 2, x, FUN = "-"), 1, FUN = function(y) sqrt(sum(y^2)))
#   SortRootSumSquare <- sort.int(RootSumSquare, decreasing = FALSE, index.return = TRUE)
#   #return(c(SortRootSumSquare$x[1:d], SortRootSumSquare$ix[1:d]))
#   return(c(RootSumSquare[SortRootSumSquare$ix[1:d]], SortRootSumSquare$ix[1:d]))
# }
# 
# CommCoord.AvgResponse <- foreach(i = 1:nrow(CommonCoord.1.1), .combine = c) %dopar% {
#   #InverseDistance <- 1/CommonCoord.DistanceAndIndex[i,1:d][CommonCoord.DistanceAndIndex[i,1:d] > 0]
#   InverseDistance <- CommonCoord.DistanceAndIndex[i,1:d][CommonCoord.DistanceAndIndex[i,1:d] > 0]
#   AvgIndex <- CommonCoord.DistanceAndIndex[i,(d+1):(2*d)][CommonCoord.DistanceAndIndex[i,1:d] > 0]
#   sum((InverseDistance/sum(InverseDistance)) * data.df[AvgIndex, 3])
#   #sum(data.df[AvgIndex, 3])/length(AvgIndex)
# }
# 
# CommCoord.AvgResponse <- foreach(i = 1:nrow(CommonCoord.1.1), .combine = c) %dopar% {
#   InverseDistance <- CommonCoord.DistanceAndIndex[i,1:d]
#   AvgIndex <- CommonCoord.DistanceAndIndex[i,(d+1):(2*d)][CommonCoord.DistanceAndIndex[i,1:d] > 0]
#   InverseDistance <- InverseDistance[InverseDistance > 0]
#   data.df[AvgIndex[which.min(InverseDistance)], 3]
# }
