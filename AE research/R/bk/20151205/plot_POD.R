#######################################################################
# Saving images - POD modes
#######################################################################
library(foreach)
library(doSNOW)
library(parallel)
library(colorRamps)
library(shape)

setwd("../../proj/jeffwu/isye-ae/RawData/")
colors <- matlab.like(numCol+1)

#Parameters
numTimes = 999
numModes = 500
ind = 9 #Flow temperature
load("../CommonCoord.RData")
#load(paste0("../rARPACK_output/POD_",numTimes,"t_",numModes,"m_",ind,"f.RData"))

numCol <- 1000 #numCol colors in legend
colors <- matlab.like(numCol+1)
#dir.create(paste0("../rARPACK_output/",numTimes,"t_",numModes,"m_",ind,"f_PODimages"))
#   Plot = foreach(j = 1:1, .packages = c("colorRamps", "shape")) %dopar% {
Range = c(-0.02,0.02)

#Use case 1 as basis case
Coordinate = rbind(CommonCoord[[1]][[1]],CommonCoord[[1]][[2]],CommonCoord[[1]][[3]],CommonCoord[[1]][[4]])

for(ind in 4:9){
  if (ind==4){
    stg = "Velocity (u)"
  }else if(ind==5){
    stg = "Velocity (v)"
  }else if(ind==6){
    stg = "Velocity (w)"
  }else if(ind==7){
    stg = "Pressure"
  }else if(ind==8){
    stg = "Temperature"
  }else if(ind==9){
    stg = "Density"
  }
  
  load(paste0("../bigmemory_bk/", ind, "/PODdat.RData"))
  na.index <- c(385,578)
  
  for (j in 1:10){  
    zcolor <- pmin(pmax((PsiTsvg[,j] - min(PsiTsvg[,j]))/(diff(range(PsiTsvg[,j])))*numCol,0),numCol) + 1
    
    print(paste0("POD mode ",j))      
    png(paste0("../POD_result/", ind, "/",numTimes,"t_",numModes,"m_",ind,"f_POD_",j,".png"), width = 600, height = 400)
    plot(Coordinate[-na.index,1],Coordinate[-na.index,2],
         col=colors[zcolor],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0(stg, " - Mode ", j)) 
    #colorlegend(col=colors,zlim=range(PsiTsvg[,j]),left=T, digit = 2)
    dev.off()
  }
}

############################################################################
# # Plots the POD modes
# ###########################################################################
# 
# library(colorRamps)
# library(shape)
# library(ggplot2)
# 
# #Loading coordinates and one POD mode
# load("output/Psi_100.RData")
# load("output/Coordinate.RData")
# qq <- 100
# 
# for (i in 1:10){
#   modeNum = i
#   mod <- PsiTsvg[,modeNum]
#   
#   #Plotting one POD mode
#   numCol <- 1000
#   zrange <- c(-0.032,0.032)
#   colors <- matlab.like(numCol+1)
#   zcolor <- colors[pmax((mod - zrange[1])/diff(range(zrange))*numCol,0) + 1] 
#   png(paste0('./output/POD figures/m',modeNum,'q',qq,'_rho.png'),width=480,height=360)
#   plot(Coordinate[,1],Coordinate[,2],col=zcolor,pch=15,cex=1.2,
#        xlab = 'x (m)', ylab = 'y (m)', main = paste0('POD mode ', modeNum, " (q=", qq, ")")) 
#   colorlegend(col=colors,zlim=zrange,digit=4,left=T)
#   dev.off()
# }

###########################################################################
# Plots the POD eigenvalues (energy)
###########################################################################
library(colorRamps)
library(shape)
library(ggplot2)  
library(grid)

#Set working directory
ind = 9
cut = 0.95
setwd(paste0("/home/proj/jeffwu/isye-ae/bigmemory_bk/", ind))
load("PODdat.RData")
eigval = pmax(eigval,rep(0,length(eigval)))#Temporary: will fix later

png(paste0("../../energy_",ind,".png"), width = 600, height = 600)
df <- data.frame(eig=1:length(eigval), energy = cumsum(eigval/sum(eigval)))
my_grob = grobTree(textGrob(paste0("10 modes - ",100*round(df$energy[10],2),"%\n",
                                   "50 modes - ",100*round(df$energy[50],2),"%\n",
                                   "100 modes - ",100*round(df$energy[100],2),"%\n",
                                   min(which(df$energy>cut))," modes - 95%"),
                            x=0.5,  y=0.6, hjust=0,gp=gpar(fontsize=20)))
ggplot(df,aes(x=eig,y=energy)) + geom_line(color="blue",size=2) + 
  xlab("Number of POD modes") + ylab("% energy captured") + ggtitle("POD energy spectrum - Temperature") + theme(plot.title=element_text(face="bold", size=20)) + 
  ylim(0,1) +
  geom_vline(xintercept=10,colour="red",size=1.2,show_guide=TRUE,linetype=2) + geom_vline(xintercept=50,colour="red",show_guide=TRUE,size=1.2,linetype=3) + geom_vline(xintercept=100,colour="red",show_guide=TRUE,size=1.2,linetype=4) +
  geom_vline(xintercept=min(which(df$energy>cut)),colour="red",show_guide=TRUE,size=1.2)+
  annotation_custom(my_grob)
dev.off()
