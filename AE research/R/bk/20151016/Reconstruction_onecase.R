#######################################################################
# Plot Reconstruction
#######################################################################

library(foreach)
library(doSNOW)
library(parallel)
library(colorRamps)
library(shape)

ind = 8 #Flow temperature
POD_mode = 50 #Number of POD modes used
qqVec = c(50, 100, 200, 500, 997)

setwd("../../proj/jeffwu/isye-ae/RawData/")
load("../Coordinate.RData")
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

#Compute average snapshot
TT = 1000;
available_file.tmp <- (1:999)[-c(1,853)] #Files with the same format
FT = foreach(i = 1:TT, .combine=cbind) %dopar% {
  if(any(i == available_file.tmp)){
    #Read snapshot at i
    snapT <- read.table(paste("1/Week1_1_",i,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
    return(snapT)
  }else {
    return(0)
  }
}
FT = FT[,available_file.tmp]


#######################################################################
# Computing reconstruction error
#######################################################################

rErr = rep(0,length(qqVec)) #To store RMSE

for(i in 1:length(qqVec)){
  cat("qq =", qqVec[i], "\n")
  load(paste0("../",POD_mode,"mode(",ind,")/BT_", qqVec[i],".RData"))
  load(paste0("../",POD_mode,"mode(",ind,")/Psi_", qqVec[i], ".RData"))
  
  TT <- length(available_file.tmp)
  Fq = PsiTsvg %*% BTsvg
  Fmean = rowMeans(FT)
  Fq = Fq + Fmean
  
#   varExpl = sum((Fq-Fmean)^2)
#   varTot = sum((FT-Fmean)^2)
  rErr[i] = sqrt(mean((FT-Fq)^2))
#   cat("Total variation explained: ", varExpl/varTot)
#   cat("Average error at given point and time: ", rErr)
}
save(rErr,file=paste0("../rErr",POD_mode,".RData"))

#######################################################################
# Reconstructing and saving flow/modes
#######################################################################

for(qq in qqVec){
  
  cat("qq =", qq, "\n")
  load(paste0("../",POD_mode,"mode(",ind,")/BT_", qq,".RData"))
  load(paste0("../",POD_mode,"mode(",ind,")/Psi_", qq, ".RData"))
  
  TT <- length(available_file.tmp)
  Fq = PsiTsvg %*% BTsvg
  Fmean = rowMeans(FT)
  Fq = Fq + Fmean
  
  #######################################################################
  # Saving images - Flow
  #######################################################################
    
  numCol <- 1000 #numCol colors in legend
  colors <- matlab.like(numCol+1)
  Range = c(120, 309.8351) # Same range with true FT
  
  Plot = foreach(i = 1:TT, .packages = c("colorRamps", "shape")) %dopar% {
    if(!dir.exists(paste0("../",POD_mode,"mode(",ind,")/q", qq, "_Snapshot"))) dir.create(paste0("../",POD_mode,"mode(",ind,")/q", qq, "_Snapshot"))
    
    png(paste0("../",POD_mode,"mode(",ind,")/q", qq, "_Snapshot/snapshot_", available_file.tmp[i], ".png"), width = 600, height = 400)
    index <- Coordinate[,1] < 0.04 & Coordinate[,2] < 0.01 #Cut off region
    zcolor <- colors[pmin(pmax((Fq[index,i] - Range[1])/(Range[2]-Range[1])*numCol,0),numCol) + 1] 
    plot(Coordinate[index,1],Coordinate[index,2],col=zcolor,pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0(paste("q =", qq,"Snapshot at T = "), available_file.tmp[i])) 
    colorlegend(col=colors,zlim=Range,left=T)
    
    dev.off()
  }
  
  
#   #######################################################################
#   # Saving images - POD modes
#   #######################################################################
#   numCol <- 1000 #numCol colors in legend
#   colors <- matlab.like(numCol+1)
#   Plot = foreach(i = 1:ncol(PsiTsvg), .packages = c("colorRamps", "shape")) %dopar% {
#     if(!dir.exists(paste0("../q", qq, "_mode"))) dir.create(paste0("../q", qq, "_mode"))
#     
#     png(paste0("../q", qq, "_mode/mode_", i, ".png"), width = 600, height = 400)
#     index <- Coordinate[,1] < 0.04 & Coordinate[,2] < 0.01
#     zcolor <- colors[(PsiTsvg[index,i] + 0.03)/0.06*numCol + 1] 
#     plot(Coordinate[index,1], Coordinate[index,2],col=zcolor,pch=15,cex=1.2,
#          xlab = 'x (m)', ylab = 'y (m)', main = paste("q =", qq, "Mode", i)) 
#     colorlegend(col=colors,zlim=c(-0.03, 0.03),left=T, digit = 3)
#     dev.off()
#   }
}

stopCluster(cl)