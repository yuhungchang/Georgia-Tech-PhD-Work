library(GPfit)
library(lhs)
source("MRGP_fit.R")
source("predict.MRGP.R")
load("../../proj/jeffwu/isye-ae/rARPACK_output/POD_200_500m.RData")
load("../../proj/jeffwu/isye-ae/CommonCoord.RData")
load("../../proj/jeffwu/isye-ae/rARPACK_output/Fmean_200.RData")
TimeStep = 201
numModes = 20
set.seed(5)
TestData = matrix(runif(5), nrow = 1)
X = read.csv("../../proj/jeffwu/isye-ae/desnorm.csv")

base = 2 #Base case
newCases = matrix(0, nrow = ncol(X), ncol = ncol(X))
for (i in 1:ncol(X)){
  newCases[i,] = as.matrix(X[base,] + 0.2 * diag(ncol(X))[i,])
}


ptm <- proc.time()
Y_hat = matrix(0, ncol = numModes, nrow = TimeStep)
Y_hat.ls = vector("list", nrow(newCases))
for(k in 1:nrow(newCases)) Y_hat.ls[[k]] = Y_hat
for(i in 1 : TimeStep){
  print(i)
  Y = t(BTsvg[1:numModes, i + TimeStep * (1:30 - 1)])
  ModelFit = MRGP_fit(X, Y)
  Pred = predict.MRGP(ModelFit, newCases)
  for(k in 1:nrow(newCases)) Y_hat.ls[[k]][i,] = Pred$Y_hat[k,]
}

print(proc.time() - ptm)

save(Y_hat.ls, file = "../../proj/jeffwu/isye-ae/MRGP_newCases_Y_hat.RData")

load("../../proj/jeffwu/isye-ae/MRGP_newCases_Y_hat.RData")
#Make the prediction

library(colorRamps)
library(shape)
load("../../proj/jeffwu/isye-ae/CommonCoord.RData")
desorig = read.csv("../../proj/jeffwu/isye-ae/30mppts_by_week.csv")[2:3]/1000
rangedes = rbind(apply(desorig,2,max),apply(desorig,2,min))

#Compute original coordinates for these cases
Coordinates = vector("list",5)
oneCase = 2

for (i in 1:5){
  Coordinates[[i]] = cbind(CommonCoord[[oneCase]][[1]][,1]/desorig[oneCase,1]*(newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                           CommonCoord[[oneCase]][[1]][,2]/desorig[oneCase,2]*(newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2]))
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind(CommonCoord[[oneCase]][[2]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                                                  CommonCoord[[oneCase]][[2]][,2]/desorig[oneCase,2]*(newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])))
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind(CommonCoord[[oneCase]][[3]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                                                  CommonCoord[[oneCase]][[3]][,2] - desorig[oneCase,2] + (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])))
}
dim(Coordinates[[1]])

for(i in 1:5){
  Y_hat = t(Y_hat.ls[[i]])
  newFlow = PsiTsvg[,1:numModes] %*% Y_hat + Fmean
  
  #Output image
  for (j in 1:100){
    Range = c(120, 309.8351) # Same range with true FTs
    # Range = c(-20,30)
    numCol = 100
    colors <- matlab.like(numCol+1)
    zcolor <- pmin(pmax((newFlow[,j] - Range[1])/(Range[2]-Range[1])*numCol,0),numCol) + 1
    png(paste0("m_pred/",i,"/newFlow_pred_", i, "_", j, ".png"), width = 600, height = 400)
    plot(Coordinates[[i]][,1],Coordinates[[i]][,2],col=colors[zcolor],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Flow prediction at T = ", j),
         xlim = c(0,0.10), ylim = c(0,0.012) )
    colorlegend(col=colors,zlim=Range,left=T)
    dev.off()
  }
}




turn.on = FALSE

if(turn.on){
  
  ptm <- proc.time()
  Y_hat = matrix(0, ncol = numModes, nrow = TimeStep)
  for(i in 1 : TimeStep){
    print(i)
    Y = t(BTsvg[1:numModes,1 + TimeStep * (1:30 - 1)])
    ModelFit = MRGP_fit(X, Y)
    Pred = predict.MRGP(ModelFit, TestData)
    Y_hat[i,] = Pred$Y_hat
  }
  
  print(proc.time() - ptm)
  
  save(Y_hat, file = "../../proj/jeffwu/isye-ae/MRGP_Y_hat.RData")
  
}

turn.on = FALSE

if(turn.on){
  
  library(foreach)
  library(doSNOW)
  library(parallel)
  
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  
  ptm <- proc.time()
  Y_hat <- foreach(i = 1 : TimeStep, .combine = rbind, .packages = "GPfit") %dopar% {
    
    Y = t(BTsvg[1:numModes,1 + TimeStep * (1:30 - 1)])
    ModelFit = MRGP_fit(X, Y)
    Pred = predict.MRGP(ModelFit, TestData)
    Pred$Y_hat
  }
  
  stopCluster(cl)
  
  proc.time() - ptm
  
  
  
  library(foreach)
  library(doSNOW)
  library(parallel)
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  
  
  Y_hat <- foreach(k = 1 : nrow(BTsvg), .combine = rbind, .packages = "GPfit") %dopar% {
    Y_hat_tmp = rep(0, length(TimeStep))
    for(i in 1L:TimeStep){
      Y = BTsvg[k,i + TimeStep * (1:30 - 1)]
      ModelFit = GP_fit(X, Y)
      Pred = predict.GP(ModelFit, TestData)
      Y_hat_tmp[i] = Pred$Y_hat
    }
    Y_hat_tmp
  }
  
  
  Y_hat <- foreach(i = 1 : TimeStep, .combine = cbind, .packages = "GPfit") %dopar% {
    Y = BTsvg[,i + TimeStep * (1:30 - 1)]
    Y_hat_tmp = rep(0, nrow(Y))
    for(k in 1L:nrow(Y)){
      ModelFit = GP_fit(X, Y[k,])
      Pred = predict.GP(ModelFit, TestData)
      Y_hat_tmp[k] = Pred$Y_hat
    }
  }
  
  stopCluster(cl)
  
  save(Y_hat, file = "../../proj/jeffwu/isye-ae/Test_Y_hat.RData")
  
}

