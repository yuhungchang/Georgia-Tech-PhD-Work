#######################################################################
# Plot Reconstruction
#######################################################################

library(foreach)
library(doSNOW)
library(parallel)
library(colorRamps)
library(shape)

#Set working directory
setwd("/home/proj/jeffwu/isye-ae/RawData/")

#Start cluster
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

#Parameters
numTimes = 200
numModes = 500
ind = 7 #Flow temperature
available_file.tmp <- 3:(3+numTimes)
available_file.tmp <- setdiff(available_file.tmp,seq(from=13,by=11,to=1000)) #Remove jumps
TT = length(available_file.tmp)

#######################################################################
# Computing reconstructed image
#######################################################################

#Loading eigenvectors and true data
load(paste0("../rARPACK_output/POD_",numTimes,"t_",numModes,"m_", ind, "f.RData"))
load(paste0("../rARPACK_output/FT_",numTimes,"t_",ind,"f.RData"))
load(paste0("../rARPACK_output/Fmean_",numTimes,"t_",ind,"f.RData"))

#Start cluster
cl <- makeCluster(detectCores())
# registerDoSNOW(cl)

qqVec = c(200)
# rErr = rep(0,length(qqVec))
for (j in 1:length(qqVec)){
  Fq = PsiTsvg[,1:qqVec[j]] %*% BTsvg[1:qqVec[j],]
  Fq = Fq + Fmean
  #   varExpl = sum((Fq-Fmean)^2)
  #   varTot = sum((FT-Fmean)^2)
  # rErr[j] = (sqrt(mean((FT-Fq)^2)))
  #   cat("Total variation explained: ", varExpl/varTot)
  #   cat("Average error at given point and time: ", rErr)
}

save(Fq,file=paste0("../rARPACK_output/Frecon_",numTimes,"_200m_",ind,"f.RData"))
# save(rErr,file=paste0("../rARPACK_output/rErr_",numTimes,"_",numModes,"m.RData"))

#######################################################################
# Saving images - Reconstructed flow
#######################################################################

#Plot
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
load("../CommonCoord.RData")

#Save reconstruction
dat <- read.table(paste0("../../RawData/", 1,"/Week",1,"_",1,"_",100,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")
for (j in 1:1){
  Coordinate = rbind(CommonCoord[[j]][[1]],CommonCoord[[j]][[2]],CommonCoord[[j]][[3]],CommonCoord[[j]][[4]])
#   Coordinate = dat[,1:2]
  
  for (i in 1:10){
#   for (i in 1:1){
    print(paste0("case ",j,", timestep ",i))  
    zcolor <- pmin( pmax( (FTmean.vector + PsiTsvg%*%matrix(BTsvg[i,],ncol=1) - Range[1])/(Range[2]-Range[1])*numCol , 0 ), numCol ) + 1
#     zcolor <- pmin( pmax( (dat[,8] - Range[1])/(Range[2]-Range[1])*numCol , 0 ), numCol ) + 1
    png(paste0("../../",ind,"-",i,".png"), width = 600, height = 400)
#     png(paste0("../../test.png"), width = 600, height = 400)
    plot(Coordinate[,1],Coordinate[,2],col=colors[zcolor],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Reconstructed at T = ", i),
         xlim = c(0,0.10), ylim = c(0,0.013) ) 
    colorlegend(col=colors,zlim=Range,left=T)
    dev.off()
  }    
}

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
