##########################################################
#         *** Plot Validation Run                        #
#             and Quantify the difference ***            #
#    MADE BY CHIH-LI and SIMON ON 1/17/2015              #               
##########################################################

#Loading libraries
library(colorRamps)
library(shape)
library(foreach)
library(doSNOW)
library(parallel)
#install.packages("bigmemory", repos="http://cran.us.r-project.org")
library(bigmemory)

### Set work directory
setwd("/home/proj/jeffwu/isye-ae/Validation Data/")

### Read data
tt <- 0
down1.start.row <- 21
dowm1.nrow <- 129 * 225
down1.dat <- read.table(paste0("Mesh_plot_",tt,".dat"), header = FALSE, 
                   skip = down1.start.row, nrow = dowm1.nrow, sep = " ", colClasses=rep("numeric",14), 
                   comment.char = "")
down2.start.row <- down1.start.row + dowm1.nrow + 5
dowm2.nrow <- 65 * 225
down2.dat <- read.table(paste0("Mesh_plot_",tt,".dat"), header = FALSE, 
                   skip = down2.start.row, nrow = dowm2.nrow, sep = " ", colClasses=rep("numeric",14), 
                   comment.char = "")
up.start.row <- down2.start.row + dowm2.nrow + 5
up.nrow <- 193 * 129
up.dat <- read.table(paste0("Mesh_plot_",tt,".dat"), header = FALSE, 
                   skip = up.start.row, nrow = up.nrow, sep = " ", colClasses=rep("numeric",14), 
                   comment.char = "")
valid.dat <- rbind(down1.dat, down2.dat, up.dat)
Coordinate.newcase <- valid.dat[,2:3]

### Geometry parameter
newL = 0.022
newR = 0.00321546
newDL = 0.00341658

### Fit to common grid
#cut outside
select.fg <- (Coordinate.newcase[,1] <= (newL + 4 * newR)) &
  (Coordinate.newcase[,2] <= (2 * newR))
Coordinate.newcase <- Coordinate.newcase[select.fg, ]
Coordinate.newcase.index <- which(select.fg)

#Normalize - Divide into Four parts
tmpList = vector("list",4)
tmpList[[1]] = which((Coordinate.newcase[,1]<=newDL) & (Coordinate.newcase[,2]<=newR))
tmpList[[2]] = which((Coordinate.newcase[,1]>newDL) & (Coordinate.newcase[,1]<=newL) & (Coordinate.newcase[,2]<=newR))
tmpList[[3]] = which((Coordinate.newcase[,1]>newL) & (Coordinate.newcase[,2]<=newR))
tmpList[[4]] = which((Coordinate.newcase[,1]>newL) & (Coordinate.newcase[,2]>newR))

CutCoordInd.newcase = tmpList

#Golden (densest) injector
load("../Coordinate.RData")
load("../CutCoordInd.RData")
des = read.csv("../30mppts_by_week.csv")
Lval = des[,2]/1000
Rval = des[,3]/1000
DLval = des[,6]/1000
case = 22
blInj1 = Coordinate[[case]][CutCoordInd[[case]][[1]],]
blInj2 = Coordinate[[case]][CutCoordInd[[case]][[2]],]
blBC = Coordinate[[case]][CutCoordInd[[case]][[3]],]
blTC = Coordinate[[case]][CutCoordInd[[case]][[4]],]

#Output common grid
tmpList = vector("list",4)

#First part: part one of inside injector (aspect ratio = 10)
tmpList[[1]] = cbind(newDL*blInj1[,1]/DLval[case], newR*blInj1[,2]/Rval[case])

#Second part: part two of inside injector (aspect ratio = 10)
tmpList[[2]] = cbind(newDL + (newL - newDL)* (blInj2[,1] - DLval[case])/(Lval[case] - DLval[case]), newR*blInj2[,2]/Rval[case])

#Third part: bottom of chamber (aspect ratio = 10)
tmpList[[3]] = cbind(newL + (blBC[,1]-Lval[case]) * newR/Rval[case], newR*blBC[,2]/Rval[case]) #edited on 1/13

#Fourth part: top of chamber (aspect ratio = cutX/cutY)
tmpList[[4]] = cbind(newL + (blTC[,1]-Lval[case]) * newR/Rval[case], newR + (blTC[,2]-Rval[case]) * newR/Rval[case]) #edited on 1/13

#   numFinalGrid = nrow(grid1) + nrow(grid2) + nrow(grid3)
CommonCoord.newcase = tmpList
save(CommonCoord.newcase, file = "CommonCoord.newcase.RData")

#Set up clusters
cl <- makeCluster(detectCores())
registerDoSNOW(cl)

Number.NN <- 10 # number of nearest neighborhood

RealCoord.i <- Coordinate.newcase
CommonCoord.i <- rbind(CommonCoord.newcase[[1]], CommonCoord.newcase[[2]],
                       CommonCoord.newcase[[3]], CommonCoord.newcase[[4]])

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
CommonCoord.newcase.MinDistIndex <- out.Matrix[,1:Number.NN]
CommonCoord.newcase.Weight <- out.Matrix[,(Number.NN + 1) : (2 * Number.NN)]

rownames(CommonCoord.newcase.MinDistIndex) <- NULL
rownames(CommonCoord.newcase.Weight) <- NULL

stopCluster(cl)

#check interpolation result
check.fg <- TRUE
if(check.fg) check.tm.num <- 10

#Choose flow property
ind = c(4:9) 
cat("ind = ", ind, "...\n")

#Timesteps available
available_file.tmp <- 0:998
TT = length(available_file.tmp)
rescaleInd_by_1000 = 7

#Set up clusters
cl <- makeCluster(detectCores(),type="SOCK")
registerDoSNOW(cl)

#Compute interpolation
weight.mx <- CommonCoord.newcase.Weight
min.dist.mx <- CommonCoord.newcase.MinDistIndex

FT = as.big.matrix(foreach(tm = available_file.tmp, .combine=cbind, .packages = c("colorRamps", "shape")) %dopar% {
  #Read snapshot for case i at timestep tm
  down1.start.row <- 21
  dowm1.nrow <- 129 * 225
  down1.dat <- read.table(paste0("Mesh_plot_",tm,".dat"), header = FALSE, 
                          skip = down1.start.row, nrow = dowm1.nrow, sep = " ", colClasses=rep("numeric",14), 
                          comment.char = "")
  down2.start.row <- down1.start.row + dowm1.nrow + 5
  dowm2.nrow <- 65 * 225
  down2.dat <- read.table(paste0("Mesh_plot_",tm,".dat"), header = FALSE, 
                          skip = down2.start.row, nrow = dowm2.nrow, sep = " ", colClasses=rep("numeric",14), 
                          comment.char = "")
  up.start.row <- down2.start.row + dowm2.nrow + 5
  up.nrow <- 193 * 129
  up.dat <- read.table(paste0("Mesh_plot_",tm,".dat"), header = FALSE, 
                       skip = up.start.row, nrow = up.nrow, sep = " ", colClasses=rep("numeric",14), 
                       comment.char = "")
  snapT <- rbind(down1.dat, down2.dat, up.dat)
  snapT <- snapT[,-1]
  snapT <- snapT[,ind]
  
  select.fg <- Coordinate.newcase.index
  snapT <- snapT[select.fg,]
  snapT[,ind == rescaleInd_by_1000] <- snapT[,ind == rescaleInd_by_1000]/1000

  #interpolation: inverse distance weighting with nearest neighborhood k method
  #               k = 10
  interpolate.mx <- matrix(0, ncol = ncol(snapT), nrow = nrow(min.dist.mx))
  for(kk in 1:ncol(snapT)){
    interpolate.val <- rep(0, nrow(min.dist.mx))
    for(num in 1:ncol(min.dist.mx)){
      interpolate.val <- interpolate.val + snapT[min.dist.mx[,num], kk] * weight.mx[,num]
    }
    interpolate.mx[,kk] <- interpolate.val
  }
  if(check.fg){
    if(tm <= check.tm.num){
      numCol <- 1000
      colors <- matlab.like(numCol+1)
      RealCoord.i <- Coordinate.newcase
      
      for(kk in 1:length(ind)){
        png(paste0("check_plot/True-", ind[kk], "-", tm, ".png"), width = 600, height = 400)
        zcolor <- colors[(snapT[,kk] - min(snapT[,kk]))/diff(range(snapT[,kk]))*numCol + 1] 
        plot(RealCoord.i[,1],RealCoord.i[,2],col=zcolor,pch=15,cex=0.5,
             xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", tm, ", ind = ", ind[kk], ", validation"))
        colorlegend(col=colors,zlim=range(snapT[,kk]),left=T)
        dev.off()
        
        for(part in 1:4){
          row.index <- CutCoordInd.newcase[[part]]
          png(paste0("check_plot/True-", ind[kk], "-", tm, "-", part, ".png"), width = 600, height = 400)
          plot(RealCoord.i[row.index,1],RealCoord.i[row.index,2],col=zcolor[row.index],pch=15,cex=0.8,
               xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", tm, ", ind = ", ind[kk], ", validation, part = ", part))
          colorlegend(col=colors,zlim=range(snapT[,kk]),left=T)
          dev.off()
        }
      }
    }
  }
  
  
  return(interpolate.mx)
}, backingfile = paste0("FT_NewCaseAllInd.bck"), backingpath = "bigmemory_bk", descriptorfile = paste0("FT_NewCaseAllInd.dsc"))

for(kk in 1:length(ind)){
  column.index <- kk + 0:(length(available_file.tmp)-1) * length(ind)
  FT.tmp <- as.big.matrix(FT[, column.index], backingfile = "FT_newcase.bck", backingpath = paste0("bigmemory_bk/", ind[kk]), descriptorfile = "FT_newcase.dsc")
}

if(check.fg){
  numCol <- 1000
  colors <- matlab.like(numCol+1)
  
  for(kk in 1:length(ind)){
    FT.val <- attach.big.matrix(paste0("bigmemory_bk/", ind[kk], "/FT_newcase.dsc"))
    CommonCoord.i <- rbind(CommonCoord.newcase[[1]], CommonCoord.newcase[[2]],
                           CommonCoord.newcase[[3]], CommonCoord.newcase[[4]])
    for(tt in 1:check.tm.num){
      numCol <- 1000
      colors <- matlab.like(numCol+1)
      png(paste0("check_plot/FT-", ind[kk], "-", tt, ".png"), width = 600, height = 400)
      zcolor <- colors[(FT.val[,tt] - min(FT.val[,tt]))/diff(range(FT.val[,tt]))*numCol + 1] 
      plot(CommonCoord.i[,1], CommonCoord.i[,2], col = zcolor, pch=15, cex=0.5,
           xlab = 'x (m)', ylab = 'y (m)', main = paste0("FT at T = ", tt, ", ind = ", ind[kk], ", validation")) 
      colorlegend(col=colors,zlim=range(FT.val[,tt]),left=T)
      dev.off()
      
      for(part in 1:4){
        CommonCoord.i.j <- CommonCoord.newcase[[part]]
        if(part == 1){
          row.index <- 1:nrow(CommonCoord.i.j)
        }else{
          row.index <- (max(row.index) + 1) : (max(row.index) + nrow(CommonCoord.i.j))       
        }
        
        numCol <- 1000
        colors <- matlab.like(numCol+1)
        png(paste0("check_plot/FT-", ind[kk], "-", tt, "-", part, ".png"), width = 600, height = 400)
        plot(CommonCoord.i.j[,1], CommonCoord.i.j[,2], col = zcolor[row.index], pch=15, cex=0.8,
             xlab = 'x (m)', ylab = 'y (m)', 
             main = paste0("FT at T = ", tt, ", ind = ", ind[kk], ", validation, part = ", part)) 
        colorlegend(col=colors,zlim=range(FT.val[,tt]),left=T)
        dev.off()
      }
    }
  }
}

stopCluster(cl)