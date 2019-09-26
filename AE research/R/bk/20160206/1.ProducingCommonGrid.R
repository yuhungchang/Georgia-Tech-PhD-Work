##########################################################
#         *** PRODUCING COMMON GRID ***                  #
#    MADE BY CHIH-LI and SIMON ON 10/5/2015              #               
#    NOTE: PLEASE READ THE MANUSCRIPT                    #
#    Reports & Presentations/Common Grid/Common Grid.pdf # 
##########################################################

# Command line for going to another version of R:
# module load R/3.2.2
# R
# versions in uranus-6: 2.13.1  2.15.0  2.15.2  2.15.2_shared_lib  3.0.1  3.1.0  3.1.2  3.2.2  Revo

#Loading libraries
library(colorRamps)
library(shape)
library(foreach)
library(doSNOW)
library(parallel)
library(bigalgebra)
library(bigpca)
library(calibrator)
library(rARPACK)

ptm <- proc.time()

### Step 0: input
setwd("RawData/")     # set working directory
des = read.csv("../DoE.csv")    # read in design points
num.case <- nrow(des)         # set number of sample cases

Lval = des[,2]/1000
Rval = des[,3]/1000
DLval = des[,6]/1000


###################################################
# Step 1: Find the densest grid among the 30 case #
###################################################

#Exporting unnormalized coordinates for all runs

Coordinate = vector("list",num.case)
Coordinate.index = vector("list",num.case)

for (i in 1:num.case){
  
  #Rename files
  week = ceiling(i/6)
  if (i == 14) week = 1
  
  Coordinate[[i]] <- read.table(paste0(i,"/Week",week,"_",i,"_",2,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,1:2]
  #Cut X-axis by Lval + 4 * Rval and Y-axis by 2 * Rval
  select.fg <- (Coordinate[[i]][,1] <= (Lval[i] + 4 * Rval[i])) &
    (Coordinate[[i]][,2] <= (2 * Rval[i]))
  Coordinate[[i]] <- Coordinate[[i]][select.fg, ]
  Coordinate.index[[i]] <- which(select.fg)
}
save(Coordinate, file = "../Coordinate.RData")
save(Coordinate.index, file = "../Coordinate.index.RData")

########## Checking Cut Grid ############# 
# for(i in 1:30){
#   png(paste0("../check/",i,".png"), width = 1200, height = 400)
#   plot(Coordinate[[i]])
#   dev.off()  
# }
##########################################

Grid.vt <- c()
for (i in 1:num.case){
  print(i)
  print(dim(Coordinate[[i]]))
  Grid.vt = c(Grid.vt,nrow(Coordinate[[i]]))
}
case <- which.max(Grid.vt)
print(case)

###################################################
# Step 2: Split the original grid into four parts #
###################################################

#case = 22 # Golden (densest) case

#Normalize - Divide into Four parts
CutCoordInd = vector("list",num.case)
for (i in 1:num.case){
  tmpList = vector("list",4)
  tmpList[[1]] = which((Coordinate[[i]][,1]<=DLval[i]) & (Coordinate[[i]][,2]<=Rval[i]))
  tmpList[[2]] = which((Coordinate[[i]][,1]>DLval[i]) & (Coordinate[[i]][,1]<=Lval[i]) & (Coordinate[[i]][,2]<=Rval[i]))
  tmpList[[3]] = which((Coordinate[[i]][,1]>Lval[i]) & (Coordinate[[i]][,2]<=Rval[i]))
  tmpList[[4]] = which((Coordinate[[i]][,1]>Lval[i]) & (Coordinate[[i]][,2]>Rval[i]))
  
  CutCoordInd[[i]] = tmpList
}
save(CutCoordInd,file="../CutCoordInd.RData")

#Golden (densest) injector
blInj1 = Coordinate[[case]][CutCoordInd[[case]][[1]],]
blInj2 = Coordinate[[case]][CutCoordInd[[case]][[2]],]
blBC = Coordinate[[case]][CutCoordInd[[case]][[3]],]
blTC = Coordinate[[case]][CutCoordInd[[case]][[4]],]

###################################################
# Step 3: Rescale each part to Golden case        #
###################################################

#Loading in factor settings


#Output common grid
CommonCoord = vector("list",num.case)
for (i in 1:num.case){
  tmpList = vector("list",num.case)
  
  #First part: part one of inside injector (aspect ratio = 10)
  tmpList[[1]] = cbind(DLval[i]*blInj1[,1]/DLval[case], Rval[i]*blInj1[,2]/Rval[case])
 
  #Second part: part two of inside injector (aspect ratio = 10)
  tmpList[[2]] = cbind(DLval[i] + (Lval[i] - DLval[i])* (blInj2[,1] - DLval[case])/(Lval[case] - DLval[case]), Rval[i]*blInj2[,2]/Rval[case])
  
  #Third part: bottom of chamber (aspect ratio = 10)
  tmpList[[3]] = cbind(Lval[i] + (blBC[,1]-Lval[case]) * Rval[i]/Rval[case], Rval[i]*blBC[,2]/Rval[case]) #edited on 1/13
  
  #Fourth part: top of chamber (aspect ratio = cutX/cutY)
  tmpList[[4]] = cbind(Lval[i] + (blTC[,1]-Lval[case]) * Rval[i]/Rval[case], Rval[i] + (blTC[,2]-Rval[case]) * Rval[i]/Rval[case]) #edited on 1/13
  
#   numFinalGrid = nrow(grid1) + nrow(grid2) + nrow(grid3)
  CommonCoord[[i]] = tmpList
}

save(CommonCoord,file="../CommonCoord.RData")

########### Check: Common Grid in Other cases##############
# load("../CommonCoord.RData")
# load("../Coordinate.RData")
# 
# for (i in 1:30){
#   print(i)
#   
#   #Rename files
#   #RealCoord = Coordinate[[i]][Coordinate[[i]][,1]< Lval[i] & Coordinate[[i]][,2]< Rval[i], ]
#   RealCoord = Coordinate[[i]]
#   png(paste0("../check/checkCommonGrid.", i , ".png"), width = 1200, height = 400)
#   par(mfrow = c(1,2))
#   plot(RealCoord[,1], RealCoord[,2], xlab = 'x (m)', ylab = 'y (m)', pch = 4) 
#   abline(v = des[i,6]/1000, col=2, lwd = 3)
#   segments(des[i,6]/1000,0,des[i,6]/1000,des[i,3]/1000, lwd = 4, col = 2)
#   abline(v = des[i,2]/1000, col=2, lwd = 3)
#   abline(h = des[i,3]/1000, col=2, lwd = 3)
#   plot(RealCoord[,1], RealCoord[,2], xlab = 'x (m)', ylab = 'y (m)', type = "n", main = "Rescaled") 
#   for(j in 1:4)  points(CommonCoord[[i]][[j]][,1], CommonCoord[[i]][[j]][,2], pch = 4)
#   abline(v = c((des[i,6] - des[i,5]/2)/1000,(des[i,6] + des[i,5]/2)/1000) , col=2)
#   dev.off()
# }  

###################################################
# Step 4: Fit common grid to rescaled grid and    #
#         then Rescale back to original scale     #
#         (interpolation)                         #
###################################################

#Set up clusters
cl <- makeCluster(detectCores())
registerDoSNOW(cl)


CommonCoord.MinDistIndex <- vector("list", num.case)
CommonCoord.Weight <- vector("list", num.case)
Number.NN <- 10 # number of nearest neighborhood

for (i in 1:num.case){
  cat("Case", i, "\n")
  
  RealCoord.i <- Coordinate[[i]]
  CommonCoord.i <- rbind(CommonCoord[[i]][[1]], CommonCoord[[i]][[2]],
                         CommonCoord[[i]][[3]], CommonCoord[[i]][[4]])
  
  ptm <- proc.time()
  coords <- t(RealCoord.i[,1:2])
  
  out.Matrix <- foreach(k = 1:nrow(CommonCoord.i), .combine = rbind) %dopar% {
    sort.out <- sort.int(colSums((coords - CommonCoord.i[k,])^2), 
                         decreasing = FALSE, index.return = TRUE)
    if(sort.out$x[1] < 10^(-16)){ #if exactly the same point, weight 1 o.w 0
      min.dist.ind <- sort.out$ix[1:Number.NN]
      output <- c(min.dist.ind, c(1, rep(0, Number.NN-1)))
    }else{
      min.dist <- sort.out$x[1:Number.NN]
      weight <- 1/min.dist
      weight <- weight/sum(weight)
      min.dist.ind <- sort.out$ix[1:Number.NN]
      output <- c(min.dist.ind, weight)
    }
    return(output)
  }
  CommonCoord.MinDistIndex[[i]] <- out.Matrix[,1:Number.NN]
  CommonCoord.Weight[[i]] <- out.Matrix[,(Number.NN + 1) : (2 * Number.NN)]
  
  rownames(CommonCoord.MinDistIndex[[i]]) <- NULL
  rownames(CommonCoord.Weight[[i]]) <- NULL

  print(proc.time() - ptm)
}

save(CommonCoord.MinDistIndex, file = "../CommonCoord.MinDistIndex.RData")
save(CommonCoord.Weight, file = "../CommonCoord.Weight.RData")

stopCluster(cl)

