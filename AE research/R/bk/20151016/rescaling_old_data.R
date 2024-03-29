#Set working directory where data is
setwd("../../proj/jeffwu/isye-ae/RawData/")

#Loading libraries
library(colorRamps)
library(shape)
library(foreach)
library(doSNOW)
library(parallel)

ptm <- proc.time()

#Set up clusters
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

#################################################################
#Exporting unnormalized coordinates for all runs
#################################################################
Coordinate = vector("list",30)
wonkyInd = c(1,3,13,14,15,16,19,20:30)

for (i in 1:30){
  #Rename files
  week = ceiling(i/6)
  if (i %in% wonkyInd){
    Coordinate[[i]] <- read.table(paste0(i,"/Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,1:2]
  }
  else{
    Coordinate[[i]] <- read.table(paste0(i,"/Ins_All_Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,1:2]
  }
}
save(Coordinate, file = "../Coordinate.RData")

#Check dimensions
minGrid = 1e9
for (i in 1:30){
  print(dim(Coordinate[[i]]))
  minGrid = min(minGrid,nrow(Coordinate[[i]]))
}
print(minGrid)

#################################################################
# Cutting simulation data into interested window
#################################################################

numCol = 1000
case = 22 # Baseline case
ind = 8 #Flow temperature
Range = c(120, 309.8351) # Same range with true FT

# Load in factor settings 
load("../Coordinate.RData")
des = read.csv("../30mppts_by_week.csv")
Lval = des[,2]/1000
Rval = des[,3]/1000
maxL = max(Lval)
maxR = max(Rval)

# How big is the simulation window?
maxX = rep(0,30)
maxY = rep(0,30)
for (i in 1:30){
  maxX[i] = max(Coordinate[[i]][,1])
  maxY[i] = max(Coordinate[[i]][,2])
}
max(maxX) #0.18 - Maximum X dim.
max(maxY) #0.019924 - Maximum Y dim.
cutX = 40/1000 #x-dim. outside injector to cut
cutY = 8/1000 #y-dim outside injector to cut

#Normalize - Divide into three parts
CutCoordInd = vector("list",30)
for (i in 1:30){
  tmpList = vector("list",3)
  tmpList[[1]] = which((Coordinate[[i]][,1]<=Lval[i]) & (Coordinate[[i]][,2]<=Rval[i]))
  tmpList[[2]] = which((Coordinate[[i]][,1]>Lval[i]) & (Coordinate[[i]][,1]<(Lval[i]+cutX)) & (Coordinate[[i]][,2]<=Rval[i]))
  tmpList[[3]] = which((Coordinate[[i]][,1]>Lval[i]) & (Coordinate[[i]][,1]<(Lval[i]+cutX)) & (Coordinate[[i]][,2]>Rval[i]) & (Coordinate[[i]][,2]<(Rval[i]+cutY)))
  CutCoordInd[[i]] = tmpList
}
save(CutCoordInd,file="../CutCoordInd.RData")

#Baseline injector
blInj = Coordinate[[case]][CutCoordInd[[case]][[1]],]
blBC = Coordinate[[case]][CutCoordInd[[case]][[2]],]
blTC = Coordinate[[case]][CutCoordInd[[case]][[3]],]

plot(blInj[,1],blInj[,2],pch=4,cex=0.3)
plot(blBC[,1],blBC[,2],pch=4,cex=0.3)
plot(blTC[,1],blTC[,2],pch=4,cex=0.3)

# #Plot
# cutInd = list(injecInd[[case]],botchamInd[[case]],topchamInd[[case]])
# snapT <- read.table(paste0("1/Ins_All_Week1_1_2.dat"), header = FALSE, skip = 1, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[cutInd,c(1,2,ind)]
# zcolor <- colors[(snapT[,3] - min(Range))/diff(range(Range))*numCol + 1] 
# #index <- Coordinate[,1] < 0.04 & Coordinate[,2] < 0.01
# plot(snapT[,1],snapT[,2],col=zcolor,pch=15,cex=1.2,
#      xlab = 'x (m)', ylab = 'y (m)', main = "True Snapshot")) 
# colorlegend(col=colors,zlim=range(Range),left=T)
# 
# dev.off()

#################################################################
# Deriving common grid
#################################################################

#Loading in factor settings
des = read.csv("../30mppts_by_week.csv")
Lval = des[,2]/1000
Rval = des[,3]/1000

#Cut lengths of chamber
cutX = 40/1000 #x-dim. outside injector to cut
cutY = 8/1000 #y-dim outside injector to cut

# Compute grid for each of the three pieces
numGrid = 60000
numPerPart = rep(numGrid/3,3)

#Output common grid
CommonCoord = vector("list",30)
for (i in 1:30){
  tmpList = vector("list",30)
  
  #First part: inside injector (aspect ratio = 10)
  tmpList[[1]] = cbind(Lval[i]*blInj[,1]/Lval[case],Rval[i]*blInj[,2]/Rval[case])
  
  #Second part: bottom of chamber (aspect ratio = 10)
  tmpList[[2]] = cbind(Lval[i] + (blBC[,1]-Lval[case]),Rval[i]*blBC[,2]/Rval[case])
  
  #Third part: top of chamber (aspect ratio = cutX/cutY)
  tmpList[[3]] = cbind(Lval[i] + (blTC[,1]-Lval[case]), Rval[i] + (blTC[,2]-Rval[case]) )
  
#   numFinalGrid = nrow(grid1) + nrow(grid2) + nrow(grid3)
  
  CommonCoord[[i]] = tmpList
}

save(CommonCoord,file="../CommonCoord.RData")

#Check:
plot(CommonCoord[[case]][[1]][,1],CommonCoord[[case]][[1]][,2], pch=4, cex=0.3)
