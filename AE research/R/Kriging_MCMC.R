##############################################################################
# MCMC kriging for flow
##############################################################################

library(randtoolbox)
library(doSNOW)

library(MCMCpack)
library(mvtnorm)
library(lattice)
library(gridExtra)
library(rgl)
library(parallel)
library(colorRamps)
library(shape)
library(foreach)

#Set working directory
setwd("/home/proj/jeffwu/isye-ae/")

##############################################################################
# Metropolis sampling

metrop <- function(log_post,xold,prop_sd){
  xnew <- rnorm(1,xold,prop_sd) #Proposal value
  lpold <- log_post(xold)
  lpnew <- log_post(xnew)
  acc <- exp(lpnew - lpold) #Acceptance probability
  if (runif(1) < acc){
    #Accept new proposal
    return (xnew);
  }
  else {
    #Reject new proposal
    return (xold);
  }
}

##############################################################################
# (Unnormalized) Log-posterior

logpost <- function(xx, ind){
  
  if ((xx <= 0)||(xx >= 1)){
    return (-Inf)
  }
  
  newrho = rho
  newrho[ind] = xx
  
  H = matrix(1,nrow=n,ncol=n)
  for (j in 1:d){
    H = H * (parRho[(i-1),j])^(4*dst[,,j])
  }
  
  cholH = chol(H)
  Hinv = chol2inv(cholH)
  
  detH = (prod(diag(cholH)))^2 #Determinant of H through Cholesky
  
  ret = -(m/2)*log(detH) - (1/2)*sum( (ymat%*%(Tinv%*%t(ymat)))*(Hinv) ) + (alpha-1)*sum(log(newrho)) + (beta-1)*sum(log(1-newrho))
  
  return (ret)
  
  #   sum( ( ymat%*%(Tinv%*%t(ymat)) ) * (Hinv) )
  #   t(y)%*%(kronecker(Hinv,Tinv)%*%y) #EQUAL
}

##############################################################################
# Prediction function
predfn = function(xx){
  # Emulator prediction
  vec <- apply(des,1,function(yy){(xx-yy)^2}) #distances
  hnew <- apply(vec,2,function(yy){prod(parRhom^(4*yy))}) #correlation
  return ( parMum + hnew%*%(Hinv%*%(ymat-mus)) )
}

# varfn = function(xx){ 
#   # Marginal emulator UQ
#   vec <- apply(x,1,function(yy){(xx-yy)^2}) #distances
#   hnew <- apply(vec,2,function(yy){prod(parRhom^(4*yy))}) #correlation
#   return ( c(1-hnew%*%(Hinv%*%hnew))*diag(parTm) )
# }

##############################################################################
# Bayesian emulator: MCMC fitting and prediction
##############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Simulation
numTimes = 200 #Number of desired timesteps
numModes = 200 #Number of desired POD modes
# available_file.tmp <- 3:(3+numTimes) #Files with the same format
# TT = length(available_file.tmp)
ind = 9

## Load in coordinates, mean, MCMC samples, POD data and design
load("CommonCoord.RData")
load(paste0("rARPACK_output/Fmean_",numTimes,"t_",ind,"f.RData"))
# load(paste0("rARPACK_output/MCMC_",numSamp,".RData"))
load(paste0("rARPACK_output/POD_",numTimes,"t_500m_", ind, "f.RData"))
TT = ncol(BTsvg)/30
des = read.csv("desnorm.csv");
desorig = read.csv("30mppts_by_week.csv")[2:6]/1000
rangedes = rbind(apply(desorig,2,max),apply(desorig,2,min))

#Assigning variables
n = 30 #Number of different control runs
m = numModes  #Number of observations per control setting
d = 5 #Number of control factors

#MCMC
burnin = 2000 #MCMC burn-in
numSamp = 5000 #MCMC samples
scal = 25 #larger scal -> stronger prior for IW
# nu = m+1
nu = scal*m
Psi = scal*diag(m)
alpha = 1
beta = 1
nois = 1e-8 #Noise for numerical stability
stepsizes = rep(0.02,d)

#Prediction
base = 2 #Base case 
newCases = matrix(0,nrow=d,ncol=d)
# for (i in 1:d){
#   newCases[i,] = as.matrix(des[base,] + 0.2*diag(d)[i,])
# }
newCases[1,] = as.matrix(des[base,] + 0.2*diag(d)[1,])
newCases[2,] = as.matrix(des[base,] + 0.6*diag(d)[2,])
newCases[3,] = as.matrix(des[base,] - 0.6846*diag(d)[3,])
newCases[4,] = as.matrix(des[base,] - 0.71567*diag(d)[4,])
newCases[5,] = as.matrix(des[base,] + 0.8*diag(d)[5,])

#Compute original coordinates for these cases
Coordinates = vector("list",d)
oneCase = 2

for (i in 1:d){

  # New dimensions
  newL = (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1])
  newR = (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])
  newdL = (newCases[i,5]*(rangedes[1,5]-rangedes[2,5]) + rangedes[2,5])
  
  # Injector start: scale by dL on x, scale by R on y
  Coordinates[[i]] = cbind(CommonCoord[[oneCase]][[1]][,1]/desorig[oneCase,5]*(newdL),
        CommonCoord[[oneCase]][[1]][,2]/desorig[oneCase,2]*(newR))
  
#   tmp = cbind(CommonCoord[[oneCase]][[1]][,1]/desorig[oneCase,5]*(newdL),
#               CommonCoord[[oneCase]][[1]][,2]/desorig[oneCase,2]*(newR))
#   apply(tmp,2,min)
#   apply(tmp,2,max)
  
  # Injector body: scale by (L-dL) on x, scale by R on y
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind((CommonCoord[[oneCase]][[2]][,1]-desorig[oneCase,5])/(desorig[oneCase,1]-desorig[oneCase,5])*(newL- newdL) + newdL,
                           CommonCoord[[oneCase]][[2]][,2]/desorig[oneCase,2]*(newR)) )
  
#   tmp = cbind((CommonCoord[[oneCase]][[2]][,1]-desorig[oneCase,5])/(desorig[oneCase,1]-desorig[oneCase,5])*(newL- newdL) + newdL,
#               CommonCoord[[oneCase]][[2]][,2]/desorig[oneCase,2]*(newR))
#   apply(tmp,2,min)
#   apply(tmp,2,max)
  
  # Bottom chamber: 
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind(CommonCoord[[oneCase]][[3]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                           CommonCoord[[oneCase]][[3]][,2]/desorig[oneCase,2]*(newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])))
  
#   tmp = cbind(CommonCoord[[oneCase]][[3]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
#               CommonCoord[[oneCase]][[3]][,2]/desorig[oneCase,2]*(newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2]))
#   apply(tmp,2,min)
#   apply(tmp,2,max)
  
  # Top chamber:
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind(CommonCoord[[oneCase]][[4]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                           CommonCoord[[oneCase]][[4]][,2] - desorig[oneCase,2] + (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])))
  
#   tmp = cbind(CommonCoord[[oneCase]][[4]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
#               CommonCoord[[oneCase]][[4]][,2] - desorig[oneCase,2] + (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2]))
#   apply(tmp,2,min)
#   apply(tmp,2,max)
  
}

# dim(Coordinates[[1]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#compute distance matrix
dst = array(0,c(n,n,d))
for (i in 1:d){
  dst[,,i] <- as.matrix(dist(des[,i]))^2
}

#Start cluster
cl <- makeCluster(detectCores())
registerDoSNOW(cl)


ptm <- proc.time()

#Run MCMC and prediction for each timestep
# tmp = foreach(curTime = 1:TT,.packages = c("colorRamps", "shape","mvtnorm","MCMCpack")) %dopar% {
for (curTime in 180:183){
  
  # ptm <- proc.time()
  
  #Vectorize response
  y = c(BTsvg[1:m,seq(from=curTime,by=TT,length.out=n)])
  # Matrix response (n x m)
  ymat = matrix(y,nrow=n,byrow=T)
    
  ## Priors:  - [\mu] \propto 1
  #           - [T] ~ IW(\nu,\Psi)
  #           - [rho] ~i.i.d. Beta(\alpha,\beta)
  
  ## MCMC:
  parT = array(diag(m),c(m,m,numSamp)) #T - Covariance matrix for response
  parMu = matrix(rep(colMeans(ymat),numSamp),ncol=m,byrow=T) #Mu - mean vector for each response
  parRho = matrix(0.5,nrow=numSamp,ncol=d) #Rho - correlation vector
  
  for (i in 2:numSamp){
    
    print(paste0(curTime,":",i))
    
    #H - Correlation matrix for control
    H = matrix(1,nrow=n,ncol=n)
    for (j in 1:d){
      H = H * (parRho[(i-1),j])^(4*dst[,,j])
    }
  #   print("Inverting H: ")
    tryCatch({
      cholH = chol(H+nois*diag(n))
      Hinv = chol2inv(cholH)
    }, error = function(e){})

  
    # Update mu
    tryCatch({
      parMu[i,] = rmvnorm( 1 , ( t(ymat)%*%rowSums(Hinv) ) / sum(Hinv), (parT[,,i-1]+nois*diag(m)) / sum(Hinv) )
    }, error=function(e){
      parMu[i,] = parMu[i-1,]
    })
    
    
    # Update T
    tmpMat = ymat - matrix(rep(parMu[i,],n),nrow=n,byrow=T)
    covMat = t(tmpMat)%*%(Hinv%*%tmpMat)
    covMat = (covMat+t(covMat))/2 #Symmetrify
    tryCatch({
      parT[,,i] = riwish( nu+n, covMat + Psi );
    }, error=function(e){
      parT[,,i] = parT[,,i-1]
    })
    parT[,,i] = ( parT[,,i]+t(parT[,,i]) )/2#Symmetrify
    
    # Metropolis steps for rho
  #   print("Inverting T: ")
    tryCatch({
      Tinv = chol2inv(chol(parT[,,i]+nois*diag(m)))
    }, error=function(e){
      
    })
    parRho[i,] = parRho[i-1,]
    for (j in 1:d){
      #Update H with current rho
      rho = parRho[i,]
      
      newlp <- function(xx){logpost(xx,j)}
      
      parRho[i,j] = metrop(newlp,parRho[i,j],stepsizes[j])
    }
    
  }
  
  # proc.time() - ptm
  #Timing: 1 sample = 1/17 seconds (200 modes)
  
#   ##Save MCMC samples
#   save(parT,parMu,parRho,scal,nu,Psi,alpha,beta,stepsizes,file=paste0("rARPACK_output/200t_200m_MCMC_",numSamp,".RData"))
  
  ##############################################################################
  # Bayesian joint emulator: Prediction and UQ - Plug-in for rho
  ##############################################################################
  
  ## Posterior means
  parMum = colMeans(parMu[(burnin+1):numSamp,]) #mu
  parTm = apply(parT[,,(burnin+1):numSamp],c(1,2),median) #T
  parRhom = colMeans(parRho[(burnin+1):numSamp,]) #rho
  
  # Recomputing H with plug-in rho
  H = matrix(1,nrow=n,ncol=n)
  for (j in 1:d){
    H = H * (parRhom[j])^(4*dst[,,j])
  }
  cholH = chol(H)
  Hinv = chol2inv(cholH)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Computing prediction and UQ

  #Make the prediction
  mus = matrix(rep(parMum,n),nrow=n,byrow=T)
  predval = apply(newCases,1,predfn) #Plug-in prediction of POD coefficients
  # predse = sqrt( apply(grdm,1,varfn) ); #Prediction st. err.
  newFlow = PsiTsvg[,1:m]%*%predval+Fmean
  
  #Output data
  save(newFlow,file=paste0("rARPACK_output/200t_200m_pred/",ind,"f/newFlow_pred_",curTime,"_dat.Rdata"))
  # load(paste0("rARPACK_output/",numModes,"m_pred (",ind,")/newFlow_pred_",curTime,"_dat.Rdata"))

#   load("rARPACK_output/FT_200t_4f.RData")
  # Min:max of each response
  # ind = 4: -101:69
  # ind = 5: -46:56
  # ind = 6: -65:12
  # ind = 7: -45:56
  # ind = 8: 120:310
  # ind = 9: 125:1009
  
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
  
    
  #Output image
  for (j in 1:d){
    
    #Save image
    numCol = 100
    colors <- matlab.like(numCol+1)
    zcolor <- pmin(pmax((newFlow[,j] - Range[1])/(Range[2]-Range[1])*numCol,0),numCol) + 1
    png(paste0("rARPACK_output/200t_200m_pred/",ind,"f/",j,"/newFlow_pred_",curTime,".png"), width = 600, height = 400)
    plot(Coordinates[[j]][,1],Coordinates[[j]][,2],col=colors[zcolor],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Flow prediction at T = ", curTime),
         xlim = c(0,0.10), ylim = c(0,0.013) )
    colorlegend(col=colors,zlim=Range,left=T)
    dev.off()
    
    #Save txt
    svg = cbind(Coordinates[[j]][,1],Coordinates[[j]][,2],newFlow[,j])
    colnames(svg) = c("x","y",paste0("flow",ind))
    write.csv(svg,row.names=F,eol="\r\n",
              file=paste0("rARPACK_output/200t_200m_pred/",ind,"f/",j,"/newFlow_pred",j,"_",ind,"f_",curTime,".dat"))
  }
  
  # proc.time() - ptm
  
}

proc.time() - ptm

