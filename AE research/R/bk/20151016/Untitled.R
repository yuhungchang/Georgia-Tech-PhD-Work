#Plotting true flow
load("../Coordinate.RData")
available_file.tmp <- 1:999 #Files with the same format
ind = 8 #Flow temperature
#index <- Coordinate[,1] < 0.04 & Coordinate[,2] < 0.01
index <- rep(TRUE, nrow(Coordinate))
cl <- makeCluster(detectCores())
registerDoSNOW(cl)
Range = foreach(i = available_file.tmp, .combine=c) %dopar% {
  #Read snapshot at i
  range(read.table(paste("17/Ins_All_Week3_17_",i,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[index,ind])
}

# Range = 120.0000 309.8351

library(colorRamps)
library(shape)

numCol <- 1000
colors <- matlab.like(numCol+1)
load("../Coordinate.RData")

Plot = foreach(i = available_file.tmp, .packages = c("colorRamps", "shape")) %dopar% {
  png(paste0("../True_Snapshot_New/snapshot_", i, ".png"), width = 600, height = 400)
  snapT <- read.table(paste("19/Week4_19_",i,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[index,ind]
  zcolor <- colors[(snapT - min(Range))/diff(range(Range))*numCol + 1] 
  plot(Coordinate[index,1],Coordinate[index,2],col=zcolor,pch=15,cex=1.2,
       xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", i)) 
  colorlegend(col=colors,zlim=range(Range),left=T)
  
  dev.off()
}
stopCluster(cl)
