install.packages("MaxPro")
library(MaxPro)

#Generating points
numIter = 100 #repeat for numIter to find best points
n = 18 #number of points
p = 2 #dimensions

InitialDesign<-MaxProLHD(18,2)$Design #global search for LHD
DOX <- MaxPro(InitialDesign,iteration = 1000) #local search from opt LHD
D = DOX #design points

for (i in 2:numIter){
  print(i)
  InitialDesign <- MaxProLHD(18,2)$Design #global search for LHD
  DOX <- MaxPro(InitialDesign, iteration=1000) #local search from opt LHD 
  if (DOX$measure < D$measure){
    D = DOX
  } #design points
}

#plot
plot(1, type="n", xlab=expression(x[1]), ylab=expression(x[2]), xlim = c(0,1), ylim = c(0,1), cex.lab=1.25)
points(D$Design,col=4,pch=8)
title("2-factor design points")

#save
write.csv(D$Design,file="18mppts_v2.csv")