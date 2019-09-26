library(colorRamps)
library(shape)
library(foreach)
library(doSNOW)
library(parallel)

#Plotting true flow
load("../Coordinate.RData")
available_file.tmp <- (1:999)[-c(1,853)] #Files with the same format
cl <- makeCluster(detectCores())
registerDoSNOW(cl)
# Range = foreach(i = available_file.tmp, .combine=c) %dopar% {
#   #Read snapshot at i
#   range(read.table(paste("1/Ins_All_Week1_1_",i,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[index,ind])
# }

Range = c(120, 309.8351) # Same range with true FT

ind = 4 #Flow temperature
numCol <- 1000
colors <- matlab.like(numCol+1)

Plot = foreach(i = available_file.tmp, .packages = c("colorRamps", "shape")) %dopar% {
  png(paste0("../True_Snapshot_4th_Run/snapshot_", i, ".png"), width = 600, height = 400)
  snapT <- read.table(paste("4/Ins_All_Week1_4_",i,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,c(1,2,ind)]
  zcolor <- colors[(snapT[,3] - min(Range))/diff(range(Range))*numCol + 1] 
  Coordinate <- snapT[,1:2]
  #index <- Coordinate[,1] < 0.04 & Coordinate[,2] < 0.01
  index <- rep(TRUE, nrow(Coordinate))
  plot(Coordinate[index,1],Coordinate[index,2],col=zcolor,pch=15,cex=1.2,
       xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", i)) 
  colorlegend(col=colors,zlim=range(Range),left=T)
  
  dev.off()
}
stopCluster(cl)


for(ind in 4){
  j = 26
  for(i in 1:10){
    png(paste0("../check/", ind, "_snapshot__", i, ".png"), width = 600, height = 400)
    week = ceiling(j/6)
    if(j == 14) week = 1 
    snapT <- read.table(paste0(j,"/Week",week,"_",j,"_",i,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,c(1,2,ind)]
    zcolor <- colors[(snapT[,3] - min(snapT[,3]))/diff(range(snapT[,3]))*numCol + 1] 
    Coordinate <- snapT[,1:2]

    #index <- Coordinate[,1] < 0.04 & Coordinate[,2] < 0.01
    index <- rep(TRUE, nrow(Coordinate))
    plot(Coordinate[index,1],Coordinate[index,2],col=zcolor,pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", i)) 
    colorlegend(col=colors,zlim=range(snapT[,3]),left=T)
    
    dev.off()
  }
}

library(colorRamps)
library(shape)
numCol <- 1000
colors <- matlab.like(numCol+1)

for(ind in 4){
  j = 25
  snapT <- attach.big.matrix(paste0("../bigmemory_bk/", ind, "/FT_",j,".dsc"))
  for(i in 21:30){
    png(paste0("../check/", ind, "_newsnapshot__", i, ".png"), width = 600, height = 400)
    week = ceiling(j/6)
    if(j == 14) week = 1 
    zcolor <- colors[(snapT[,i] - min(snapT[,i]))/diff(range(snapT[,i]))*numCol + 1] 
    
    Coordinate = rbind(CommonCoord[[j]][[1]],CommonCoord[[j]][[2]],CommonCoord[[j]][[3]],CommonCoord[[j]][[4]])
    #na.index <- c(385,578)
    plot(Coordinate[,1],Coordinate[,2],col=zcolor,pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", i)) 
    colorlegend(col=colors,zlim=range(snapT[,3]),left=T)
    
    dev.off()
  }
}
