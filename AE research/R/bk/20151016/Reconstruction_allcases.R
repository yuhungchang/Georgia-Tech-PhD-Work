#######################################################################
# Plot Reconstruction
#######################################################################

library(foreach)
library(doSNOW)
library(parallel)
library(colorRamps)
library(shape)

#Set working directory
setwd("../../proj/jeffwu/isye-ae/RawData/")

#Start cluster
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

#Parameters
numTimes = 200
numModes = 500
ind = 8 #Flow temperature

# #######################################################################
# #Compute average snapshot
# #######################################################################
# wonkyInd = c(1,3,13,14,15,16,19,20:30)
# 
# ptm <- proc.time()
# FT = foreach(j = 1:2, .combine='+') %dopar% {
#   #j indexes case
#   coord1 = CommonCoord.MinDistIndex[[j]][[1]]
#   coord2 = CommonCoord.MinDistIndex[[j]][[2]]
#   coord3 = CommonCoord.MinDistIndex[[j]][[3]]
#   len = length(coord1) + length(coord2) + length(coord3)
#   
#   runsum = rep(0,len)
#   for (i in 1:TT){
#     tm = available_file.tmp[i]
#     
#     #Read snapshot for case j at timestep tm
#     week = ceiling(j/6)
#     if (j %in% wonkyInd){
#       snapT <- read.table(paste0(j,"/Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
#     }
#     else{
#       snapT <- read.table(paste0(j,"/Ins_All_Week",week,"_",j,"_",tm,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
#     }
#     
#     runsum = runsum + c(snapT[coord1],snapT[coord2],snapT[coord3])
#   }
#   return(runsum)
# }
# proc.time() - ptm
# 
# Fmean = (FT)/(30*TT)
# 
# #Timing:
# #TT = 100 -> 126 s

#######################################################################
# Computing reconstruction error
#######################################################################

#Parameters
available_file.tmp <- 3:(3+numTimes) #Files with the same format
TT = length(available_file.tmp)

#Testing different number of POD modes
qqVec = c(25,50,75,100,200,300,400,500)
# rErr = rep(0,length(qqVec)) #To store RMSE

#Loading eigenvectors and true data
load(paste0("../rARPACK_output/POD_",numTimes,"_",numModes,"m.RData"))
load(paste0("../rARPACK_output/Fmean_",numTimes,".RData"))
load(paste0("../rARPACK_output/FT_",numTimes,".RData"))

#Start cluster
cl <- makeCluster(detectCores())
# registerDoSNOW(cl)

rErr = rep(0,length(qqVec))

for (j in 1:length(qqVec)){
  Fq = PsiTsvg[,1:qqVec[j]] %*% BTsvg[1:qqVec[j],]
  Fq = Fq + Fmean
  #   varExpl = sum((Fq-Fmean)^2)
  #   varTot = sum((FT-Fmean)^2)
  rErr[j] = (sqrt(mean((FT-Fq)^2)))
  #   cat("Total variation explained: ", varExpl/varTot)
  #   cat("Average error at given point and time: ", rErr)
}

# save(Fq,file=paste0("../rARPACK_output/Frecon_",numTimes,"_",numModes,"m.RData"))
save(rErr,file=paste0("../rARPACK_output/rErr_",numTimes,"_",numModes,"m.RData"))

#######################################################################
# Reconstructing and saving flow/modes
#######################################################################

#######################################################################
# Saving images - Flow
#######################################################################

#Plot
numCol <- 1000 #numCol colors in legend
colors <- matlab.like(numCol+1)
Range = c(120, 309.8351) # Same range with true FT

#Load common coordinate indices and data
load("../CommonCoord.RData")
load(paste0("../rARPACK_output/Frecon_",numTimes,"_",numModes,"m.RData"))
zcolor <- pmin(pmax((Fq - Range[1])/(Range[2]-Range[1])*numCol,0),numCol) + 1
  
#Start cluster
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

#Create directories and 
dir.create(paste0("../rARPACK_output/",numModes,"m_recon (",ind,")/"))
#   Plot = foreach(j = 1:1, .packages = c("colorRamps", "shape")) %dopar% {
for (j in 2:2){
  Coordinate = rbind(CommonCoord[[j]][[1]],CommonCoord[[j]][[2]],CommonCoord[[j]][[3]])
  #       if(!dir.exists(paste0("../rARPACK_output/",numModes,"m_recon (",ind,")/", j, "/"))){
  # dir.create(paste0("../rARPACK_output/",numModes,"m_recon (",ind,")/", j, "/"))
  #       }
  
  for (i in 1:TT){
    print(paste0("case ",j,", timestep ",i))      
    png(paste0("../rARPACK_output/",numModes,"m_recon (",ind,")/", j, "/", j, "_Snapshot_recon_", available_file.tmp[i], ".png"), width = 600, height = 400)
    plot(Coordinate[,1],Coordinate[,2],col=colors[zcolor[,TT*(j-1)+i]],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Reconstructed at T = ", available_file.tmp[i]),
         xlim = c(0,0.10), ylim = c(0,0.012) ) 
    colorlegend(col=colors,zlim=Range,left=T)
    dev.off()
  }    
}
  
#######################################################################
# Saving images - POD modes
#######################################################################

#Parameters
numTimes = 200
numModes = 500
ind = 8 #Flow temperature
available_file.tmp <- 3:(3+numTimes) #Files with the same format
TT = length(available_file.tmp)

load("../CommonCoord.RData")
load(paste0("../rARPACK_output/POD_",numTimes,"_",numModes,"m.RData"))

numCol <- 1000 #numCol colors in legend
colors <- matlab.like(numCol+1)
dir.create(paste0("../rARPACK_output/",numModes,"m_POD (",ind,")/"))
#   Plot = foreach(j = 1:1, .packages = c("colorRamps", "shape")) %dopar% {
Range = c(-0.035,0.035)

#Use case 1 as basis case
Coordinate = rbind(CommonCoord[[1]][[1]],CommonCoord[[1]][[2]],CommonCoord[[1]][[3]])

for (j in 1:100){  
  zcolor <- pmin(pmax((PsiTsvg[,j] - Range[1])/(Range[2]-Range[1])*numCol,0),numCol) + 1
  
  print(paste0("POD mode ",j))      
  png(paste0("../rARPACK_output/",numModes,"m_POD (",ind,")/", numModes, "m_POD",j, ".png"), width = 600, height = 400)
  plot(Coordinate[,1],Coordinate[,2],col=colors[zcolor],pch=15,cex=1.2,
       xlab = 'x (m)', ylab = 'y (m)', main = paste0("POD mode ", j)) 
  colorlegend(col=colors,zlim=Range,left=T)
  dev.off()
}

stopCluster(cl)