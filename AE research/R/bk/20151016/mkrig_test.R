##############################################################################
# Multiple kriging: a quick test of MLE and Bayesian emulators
##############################################################################

library(randtoolbox)
library(MCMCpack)
library(mvtnorm)
library(lattice)
library(gridExtra)
library(rgl)

##############################################################################
# 2-D Test function: Franke (scaled to [0,1]^2)
##############################################################################
frankesh <- function(xx,shift)
{
  #Scales function to [0,1]^2
  if ((xx[1]<0)|(xx[1]>1)|(xx[2]<0)|(xx[2]>1)){
    return (0)
  }
  else{
    
    lower1 = shift
    upper1 = 1+shift
    lower2 = shift
    upper2 = 1+shift
    
    x1 <- xx[1]*(upper1-lower1) + lower1
    x2 <- xx[2]*(upper2-lower2) + lower2
    
    term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
    term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
    term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
    term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)
    
    y <- term1 + term2 + term3 + term4
    return(y)
  }
}

# #Testing plot
# shift = 0.05
# numPred = 40 #number of prediction in each dim
# xvect = seq(from=0,to=1,length.out=numPred);
# yvect = seq(from=0,to=1,length.out=numPred);
# grd = expand.grid(x=xvect, y=yvect);
# grd$truefn = apply(grd,1,frankesh,shift);
# z = matrix(grd$truefn,nrow=numPred)
# filled.contour(xvect,yvect,z)

n = 30 #Number of different control runs
m = 3  #Number of observations per control setting
d = 2  #Number of control factors
x = sobol(n,2) #Control settings
dst1 = as.matrix(dist(x[,1]))  #Distance matrix for dim. 1
dst2 = as.matrix(dist(x[,2]))  #Distance matrix for dim. 1

y1 = apply(x,1,frankesh,0) #First response: no shift
y2 = apply(x,1,frankesh,0.05) #Second response: shifted 0.05
y3 = -apply(x,1,frankesh,-0.05) #Third response: shifted -0.05 and negated
y = c(rbind(y1,y2,y3)) #Vectorized response: (y11, ..., y1m, ..., yn1, ..., ynm)

##############################################################################
# Metropolis sampling
##############################################################################
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

# #Testing metrop
# numSamp = 1000
# samp = rep(0,numSamp)
# fn <- function(xx){log(dnorm(xx))}
# prop_sd <- 0.5
# for (i in 2:numSamp){
#   samp[i] <- metrop(fn,samp[i-1],prop_sd)
# }
# hist(samp,breaks=50)

##############################################################################
# (Unnormalized) Log-posterior
##############################################################################

logpost <- function(xx, ind){
  
  if ((xx <= 0)||(xx >= 1)){
    return (-Inf)
  }
  
  newrho = rho
  newrho[ind] = xx
  
  H = (newrho[1])^(4*dst1) * (newrho[2])^(4*dst2)
  cholH = chol(H)
  Hinv = chol2inv(cholH)
  
  detH = (prod(diag(cholH)))^2 #Determinant of H through Cholesky
  
  ret = -(m/2)*log(detH) - (1/2)*sum( (ymat%*%(Tinv%*%t(ymat)))*(Hinv) ) + (alpha-1)*sum(log(newrho)) + (beta-1)*sum(log(1-newrho))
    
  return (ret)
    
#   sum( ( ymat%*%(Tinv%*%t(ymat)) ) * (Hinv) )
#   t(y)%*%(kronecker(Hinv,Tinv)%*%y) #EQUAL
}

##############################################################################
# Bayesian joint emulator: MCMC fitting
##############################################################################
n = 30 #Number of different control runs
m = 3  #Number of observations per control setting
d = 2  #Number of control factors
x = mMdes = as.matrix(read.csv(sprintf("./test_designs/mM2-%d.csv",n))[,2:3]) #Control settings
# plot(1, type="n", xlab=expression(x[1]), ylab=expression(x[2]), xlim = c(0,1), ylim = c(0,1), cex.lab=1.25)
# points(x,pch=16,cex=1.25,col="black",lwd=2)
dst1 = as.matrix(dist(x[,1]))^2  #Squared-distance matrix for dim. 1
dst2 = as.matrix(dist(x[,2]))^2  #Squared-distance matrix for dim. 2

fn1 <- function(xx){frankesh(xx,0)} #First response: no shift
fn2 <- function(xx){frankesh(xx,0.05)} #Second response: shifted 0.05
fn3 <- function(xx){-frankesh(xx,-0.05)} #Third response: shifted -0.05 and negated
y1 = apply(x,1,fn1) 
y2 = apply(x,1,fn2)
y3 = apply(x,1,fn3)
y = c(rbind(y1,y2,y3)) #Vectorized response: (y11, ..., y1m, ..., yn1, ..., ynm)
ymat = cbind(y1,y2,y3) #Matrix response (n x m)
  
## Priors:  - [\mu] \propto 1
#           - [T] ~ IW(\nu,\Psi)
#           - [rho] ~i.i.d. Beta(\alpha,\beta)

nu = m+1
Psi = diag(m)
alpha = 1
beta = 1
stepsizes = rep(0.02,d)

## MCMC:
numSamp = 10000
parT = array(diag(m),c(m,m,numSamp)) #T - Covariance matrix for response
parMu = matrix(0,nrow=numSamp,ncol=m) #Mu - mean vector for each response
parRho = matrix(0.5,nrow=numSamp,ncol=d) #Rho - correlation vector

for (i in 2:numSamp){
  
  #H - Correlation matrix for control
  H = (parRho[(i-1),1])^(4*dst1) * (parRho[(i-1),2])^(4*dst2)
  cholH = chol(H)
  Hinv = chol2inv(cholH)
  
  # Update mu
  parMu[i,] = rmvnorm( 1 , ( t(ymat)%*%rowSums(Hinv) ) / sum(Hinv), 
                        parT[,,i-1] / sum(Hinv) )
  
  # Update T
  tmpMat = ymat - matrix(rep(parMu[i,],n),nrow=n,byrow=T)
  parT[,,i] = riwish( nu+n, t(tmpMat)%*%(Hinv%*%tmpMat) + Psi );
  
  # Metropolis steps for rho
  Tinv = chol2inv(chol(parT[,,i]))
  parRho[i,] = parRho[i-1,]
  for (j in 1:d){
    #Update H with current rho
    rho = parRho[i,]
    
    newlp <- function(xx){logpost(xx,j)}
    
    parRho[i,j] = metrop(newlp,parRho[i,j],stepsizes[j])
  }
  
}

## Posterior means
burnin = 5000
parMum = colMeans(parMu[(burnin+1):numSamp,]) #mu
parTm = apply(parT[,,(burnin+1):numSamp],c(1,2),median) #T
parRhom = colMeans(parRho[(burnin+1):numSamp,]) #rho

##############################################################################
# Bayesian joint emulator: Prediction and UQ - Fully Bayesian
##############################################################################


##############################################################################
# Bayesian joint emulator: Prediction and UQ - Plug-in for rho
##############################################################################
# Plot parameters
ind = "1" #1 - first response, etc.
sebd = 0.5 #bound for plotting standard errors

# Recomputing H with plug-in rho
H = (parRhom[1])^(4*dst1) * (parRhom[2])^(4*dst2)
cholH = chol(H)
Hinv = chol2inv(cholH)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computing prediction and UQ
numPred = 40 #number of prediction in each dim
xvect = seq(from=0,to=1,length.out=numPred);
yvect = seq(from=0,to=1,length.out=numPred);
grdm = expand.grid(x1=xvect, x2=yvect);

predfn = function(xx){
  # Emulator prediction
  vec <- apply(x,1,function(yy){(xx-yy)^2}) #distances
  hnew <- apply(vec,2,function(yy){prod(parRhom^(4*yy))}) #correlation
  return ( parMum + hnew%*%(Hinv%*%(ymat-mus)) )
}

varfn = function(xx){ 
  # Marginal emulator UQ
    vec <- apply(x,1,function(yy){(xx-yy)^2}) #distances
    hnew <- apply(vec,2,function(yy){prod(parRhom^(4*yy))}) #correlation
    return ( c(1-hnew%*%(Hinv%*%hnew))*diag(parTm) )
}

trueval = rbind(apply(grdm,1,fn1),apply(grdm,1,fn2),apply(grdm,1,fn3))
mus = matrix(rep(parMum,n),nrow=n,byrow=T)
predval = apply(grdm,1,predfn)
predse = sqrt( apply(grdm,1,varfn) ); #Prediction st. err.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prediction s.e. for each of the three response
grd = cbind(grdm,
            y1=trueval[1,],y2=trueval[2,],y3=trueval[3,],
            yfit1=predval[1,],yfit2=predval[2,],yfit3=predval[3,],
            se1=predse[1,],se2=predse[2,],se3=predse[3,])

#True response
switch(ind,"1"={yvec=grd$y1},"2"={yvec=grd$y2},"3"={yvec=grd$y3})
plottr=wireframe(yvec~x1*x2,data = grd, 
                   scales = list(arrows = FALSE), drape = T,
                   xlab = expression(x[1]), ylab = expression(x[2]), zlab = "",
                   colorkey=FALSE,
                   screen = list(z=120,x=-60,y=0), main=list("True",cex=1.8))

#Fitted response
switch(ind,"1"={fitvec=grd$yfit1},"2"={fitvec=grd$yfit2},"3"={fitvec=grd$yfit3})
plotfit=wireframe(fitvec~grd$x1*grd$x2, 
                  scales = list(arrows = FALSE), drape = T,
                  xlab = expression(x[1]), ylab = expression(x[2]), zlab = "",
                  colorkey=FALSE,
                  screen = list(z=120,x=-60,y=0), main=list("Emulator",cex=1.8))

#Pred. S.E. plots (2 s.e. bound)
zlim = c(0,sebd)
switch(ind,"1"={sevec=grd$se1},"2"={sevec=grd$se2},"3"={sevec=grd$se3})
plotse=wireframe(sevec~x1*x2,data = grd, 
                scales = list(arrows = FALSE), drape=T,
                xlab = expression(x[1]), ylab = expression(x[2]), zlab = "",
                colorkey=FALSE, 
                zlim = zlim, screen = list(z=120,x=-60,y=0), main=list("Pred S.E.",cex=1.8))

#Plot the plots
grid.arrange(plottr,plotfit,plotse,
             nrow=1,ncol=3,byrow = F)

##############################################################################
# Bayesian independent emulator: MCMC
##############################################################################

##############################################################################
# Comparison of joint CI's between independent and correlated emulator
##############################################################################
varfnsc = function(xx){ 
  # Scaled UQ
  vec <- apply(x,1,function(yy){(xx-yy)^2}) #distances
  hnew <- apply(vec,2,function(yy){prod(parRhom^(4*yy))}) #correlation
  return ( 1-hnew%*%(Hinv%*%hnew) )
}

predsesc = sqrt( apply(grdm,1,varfnsc) );
normerr = (trueval-predval)/predsesc
normerr_mod = normerr[,which(colSums(abs(normerr)>3)==0)]

open3d()
plot3d(t(normerr_mod), size=3, box=FALSE,
       xlab="y1",ylab="y2",zlab="y3")
plot3d( ellipse3d(parTm, centre=c(0,0,0),level=0.9999), col="green", alpha=0.2, add = TRUE)
plot3d( ellipse3d(det(parTm)^(1/3)*diag(3), centre=c(0,0,0),level=0.9999), col="blue", alpha=0.2, add = TRUE)


# # Plot the true error
# open3d()
# bg3d("white")
# persp3d(xvect, yvect, grd$y1-grd$yfit1, aspect=c(1, 1, 1), col = "red",
#         xlab = expression(x[1]), ylab = expression(x[2]), zlab = "Pred. error",
#         zlim = c(-sebd,sebd))
# # Plot 95% confidence intervals
# open3d()
# bg3d("white")
# persp3d(xvect, yvect, abs(grd$y2-grd$yfit2)-3*grd$se2, aspect=c(1, 1, 1), col = "red",
#         xlab = expression(x[1]), ylab = expression(x[2]), zlab = "Pred. error",
#         zlim = c(-sebd,sebd))
# #persp3d(xvect, yvect, 3*grd$se1, col = "blue", alpha = 0.5)
# # persp3d(xvect, yvect, 3*grd$se2, col = "blue", alpha = 0.5,add=T)
# # Plot evaluated points
# points3d(x[,1],x[,2],0,col="black")
# 
# # Rotate for better view
# play3d(spin3d(axis = c(0, 0, 1), rpm = 5), duration = 5)
# play3d(spin3d(axis = c(1, 0, 0), rpm = 5), duration = 5)
# play3d(spin3d(axis = c(0, 1, 0), rpm = 5), duration = 5)
# 
# plotse=wireframe(2*se3~x1*x2,data = grd, 
#                   scales = list(arrows = FALSE), 
#                   zlab = "Pred. S.E.", xlab = expression(x[1]), ylab = expression(x[2]),
#                   colorkey=FALSE, col.regions="transparent",
#                   zlim = zlim, screen = list(z=120,x=-60,y=0))
# ploterr=wireframe(abs(y3-yfit3)~x1*x2,data = grd, 
#                   scales = list(arrows = FALSE), drape = T,
#                   xlab = expression(x[1]), ylab = expression(x[2]),
#                   main=list(expression(y[3]),cex=1.8),colorkey=FALSE,
#                   zlab = "", zlim = zlim, screen = list(z=120,x=-60,y=0))
# print(ploterr,more=T)
# print(plotse,more=F)
# 
