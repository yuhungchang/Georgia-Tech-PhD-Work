##########################################################
#         *** INTERPOLATION ***                          #
#    MADE BY CHIH-LI and SIMON ON 12/6/2015              #
#    DESCRIPTION: using inverse distance weighting       #
#                 with nearest neighborhood k method     #          
#    NOTE: PLEASE READ THE MANUSCRIPT                    #
#    Reports & Presentations/Common Grid/Common Grid.pdf # 
##########################################################

#Step 0: input
#setwd("/home/proj/jeffwu/isye-ae/RawData/") 
setwd("RawData")       # set working directory
num.case <- 30         # set number of sample cases
ind = c(4:9)           # select flow responses
cat("ind = ", ind, "...\n")
available_file.tmp <- 0:998     # Read in snapshots available
TT = length(available_file.tmp) # Total number of snapshots
rescaleInd_by_1000 = 7
#rescaleInd_by_million = 13

check.fg <- TRUE  #check FT
if(check.fg){
  check.tm.num <- 10 # the check number of snapshot
  load("../CommonCoord.RData")
  load("../Coordinate.RData")
  load("../CutCoordInd.RData")
  if(!file.exists("../check/check_FT/")) {
    dir.create("../check/")
    dir.create("../check/check_FT/")
  }
} 

# create folders for output
if(!dir.exists("../bigmemory_bk")) dir.create("../bigmemory_bk")
for(kk in 1:length(ind)){
  if(!dir.exists(paste0("../bigmemory_bk/", ind[kk]))) dir.create(paste0("../bigmemory_bk/", ind[kk]))
}


#Loading libraries
library(foreach)
library(doSNOW)
library(parallel)
library(bigmemory)
if(check.fg){
  library(colorRamps)
  library(shape)
}

#Load common coordinate indices
load("../CommonCoord.MinDistIndex.RData")
load("../CommonCoord.Weight.RData")
load("../Coordinate.index.RData")

#Set up clusters
cl <- makeCluster(detectCores(),type="SOCK")
registerDoSNOW(cl)

#Compute FT

# rerun 1:6, 21
for(i in 1:num.case){
  cat("Case", i, "\n")
  ptm <- proc.time()
  
  week = ceiling(i/6)
  if(i == 14) week = 1 
  
  weight.mx <- CommonCoord.Weight[[i]]
  min.dist.mx <- CommonCoord.MinDistIndex[[i]]
  
  FT = as.big.matrix(foreach(tm = available_file.tmp, .combine=cbind, .packages = c("colorRamps", "shape")) %dopar% {
    #Read snapshot for case i at timestep tm
    #snapT <- read.table(paste0(i,"/Week", week, "_", i, "_", tm ,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
    tm <- format(tm, width = 3)
    tm <- gsub(" ", 0, tm)
    snapT <- read.table(paste0(i,"/",tm,".dat",sep=""), header = FALSE, skip = 15, sep = '\t', colClasses=rep("numeric",14), comment.char = "")[,ind]
    select.fg <- Coordinate.index[[i]]
    snapT <- snapT[select.fg,]
    snapT[,ind == rescaleInd_by_1000] <- snapT[,ind == rescaleInd_by_1000]/1000
    #snapT[,ind == rescaleInd_by_million] <- snapT[,ind == rescaleInd_by_million]/1000000
    
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
        RealCoord.i <- Coordinate[[i]]
        
        for(kk in 1:length(ind)){
          png(paste0("../check/check_FT/True_", i, "-", ind[kk], "-", tm, ".png"), width = 600, height = 400)
          zcolor <- colors[(snapT[,kk] - min(snapT[,kk]))/diff(range(snapT[,kk]))*numCol + 1] 
          plot(RealCoord.i[,1],RealCoord.i[,2],col=zcolor,pch=15,cex=0.5,
               xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", tm, ", ind = ", ind[kk], ", case = ", i))
          colorlegend(col=colors,zlim=range(snapT[,kk]),left=T)
          dev.off()
          
          for(part in 1:4){
            row.index <- CutCoordInd[[i]][[part]]
            png(paste0("../check/check_FT/True_", i, "-", ind[kk], "-", tm, "-", part, ".png"), width = 600, height = 400)
            plot(RealCoord.i[row.index,1],RealCoord.i[row.index,2],col=zcolor[row.index],pch=15,cex=0.8,
                 xlab = 'x (m)', ylab = 'y (m)', main = paste0("True Snapshot at T = ", tm, ", ind = ", ind[kk], ", case = ", i, ", part = ", part))
            colorlegend(col=colors,zlim=range(snapT[,kk]),left=T)
            dev.off()
          }
        }
      }
    }
    
    
    return(interpolate.mx)
  }, backingfile = paste0("FT_tmp.bck"), backingpath = paste0("../bigmemory_bk"), descriptorfile = paste0("FT_tmp.dsc"))
  
  for(kk in 1:length(ind)){
    column.index <- kk + 0:(length(available_file.tmp)-1) * length(ind)
    FT.tmp <- as.big.matrix(FT[, column.index], backingfile = paste0("FT_",i,".bck"), backingpath = paste0("../bigmemory_bk/", ind[kk]), descriptorfile = paste0("FT_",i,".dsc"))
  }
  
  file.remove(c("../bigmemory_bk/FT_tmp.bck", "../bigmemory_bk/FT_tmp.dsc"))
  
  if(check.fg){
    numCol <- 1000
    colors <- matlab.like(numCol+1)
    
    for(kk in 1:length(ind)){
      FT.val <- attach.big.matrix(paste0("../bigmemory_bk/", ind[kk], "/FT_",i,".dsc"))
      CommonCoord.i <- rbind(CommonCoord[[i]][[1]], CommonCoord[[i]][[2]],
                             CommonCoord[[i]][[3]], CommonCoord[[i]][[4]])
      for(tt in 1:check.tm.num){
        numCol <- 1000
        colors <- matlab.like(numCol+1)
        png(paste0("../check/check_FT/FT_", i , "-", ind[kk], "-", tt, ".png"), width = 600, height = 400)
        zcolor <- colors[(FT.val[,tt] - min(FT.val[,tt]))/diff(range(FT.val[,tt]))*numCol + 1] 
        plot(CommonCoord.i[,1], CommonCoord.i[,2], col = zcolor, pch=15, cex=0.5,
             xlab = 'x (m)', ylab = 'y (m)', main = paste0("FT at T = ", tt, ", ind = ", ind[kk], ", case = ", i)) 
        colorlegend(col=colors,zlim=range(FT.val[,tt]),left=T)
        dev.off()
        
        for(part in 1:4){
          CommonCoord.i.j <- CommonCoord[[i]][[part]]
          if(part == 1){
            row.index <- 1:nrow(CommonCoord.i.j)
          }else{
            row.index <- (max(row.index) + 1) : (max(row.index) + nrow(CommonCoord.i.j))       
          }
          
          numCol <- 1000
          colors <- matlab.like(numCol+1)
          png(paste0("../check/check_FT/FT_", i , "-", ind[kk], "-", tt, "-", part, ".png"), width = 600, height = 400)
          plot(CommonCoord.i.j[,1], CommonCoord.i.j[,2], col = zcolor[row.index], pch=15, cex=0.8,
               xlab = 'x (m)', ylab = 'y (m)', 
               main = paste0("FT at T = ", tt, ", ind = ", ind[kk], ", case = ", i, ", part = ", part)) 
          colorlegend(col=colors,zlim=range(FT.val[,tt]),left=T)
          dev.off()
        }
      }
    }
  }
  print(proc.time() - ptm)
}

stopCluster(cl)

