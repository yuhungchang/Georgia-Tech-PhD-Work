library(GPfit)
library(randtoolbox)
library(nloptr)
library(lattice)
library(BB)

############################################################################
#Import data
############################################################################
#set working directory to folder containing this file
setwd("C:/Users/SimonLenovo/Dropbox/Research Projects/Aerospace injector/R")
dat = read.csv("../Data/export.csv");
des = read.csv("../Data/desnorm.csv");
factNames = names(des);

############################################################################
#Initial Kriging
############################################################################
nRuns = 30;

attach(dat)
curdat = dat[time==50,c(3,4)];
attach(curdat);
attach(des);
hfit = GP_fit(des[1:nRuns,],h); #fitted kriging model for h
anglefit = GP_fit(des[1:nRuns,],angle); # ... for angle

summary.GP(hfit)
summary.GP(anglefit)

############################################################################
# Refitting kriging removing inert factors
############################################################################

# #Factors 1 and 2 inert for h
# #Factor 1 inert for angle
# hfit = GP_fit(des[1:nRuns,3:5],h);
# anglefit = GP_fit(des[1:nRuns,2:5],angle); # ... for angle
#
# summary.GP(hfit)
# summary.GP(anglefit)

#########################################################
#Feasibility (inverse problem)
#########################################################

#########################################################
#R^{-1} preprocessing for speed in hfit
M = 1
X_h = hfit$X
Y_h = hfit$Y
n = nrow(X_h)
d = ncol(X_h)
corr_h = hfit$correlation_param
power = corr_h$power
nu_h = corr_h$nu
beta_h = hfit$beta
sig2_h = hfit$sig2
delta_h = hfit$delta
R_h = corr_matrix(X_h, beta_h, corr_h)
One = rep(1, n)
LO = diag(n)
Sig_h = R_h + delta_h * LO
L_h = chol(Sig_h)
if (delta_h == 0) {
  Sig_invOne_h = solve(L_h, solve(t(L_h), One))
  Sig_invY_h = solve(L_h, solve(t(L_h), Y_h))
  Sig_invLp_h = solve(L_h, solve(t(L_h), LO))
}  else {
  s_One_h = One
  s_Yi_h = Y_h
  s_Li_h = LO
  Sig_invOne_h = matrix(0, ncol = 1, nrow = n)
  Sig_invY_h = matrix(0, ncol = 1, nrow = n)
  Sig_invLp_h = matrix(0, ncol = n, nrow = n)
  for (it in 1:M) {
    s_Onei_h = solve(L_h, solve(t(L_h), delta_h * s_Onei_h))
    Sig_invOne_h = Sig_invOne_h + s_Onei_h/delta_h
    s_Yi_h = solve(L_h, solve(t(L_h), delta_h * s_Yi_h))
    Sig_invY_h = Sig_invY_h + s_Yi_h/delta_h
    s_Li_h = solve(L_h, solve(t(L_h), delta_h * s_Li_h))
    Sig_invLp_h = Sig_invLp_h + s_Li/delta_h
  }
}
#########################################################
#R^{-1} preprocessing for speed in anglefit
M = 1
X_a = anglefit$X
Y_a = anglefit$Y
n = nrow(X_a)
d = ncol(X_a)
corr_a = anglefit$correlation_param
power = corr_a$power
nu_a = corr_a$nu
beta_a = anglefit$beta
sig2_a = anglefit$sig2
delta_a = anglefit$delta
R_a = corr_matrix(X_a, beta_a, corr_a)
One = rep(1, n)
LO = diag(n)
Sig_a = R_a + delta_a * LO
L_a = chol(Sig_a)
if (delta_a == 0) {
  Sig_invOne_a = solve(L_a, solve(t(L_a), One))
  Sig_invY_a = solve(L_a, solve(t(L_a), Y_a))
  Sig_invLp_a = solve(L_a, solve(t(L_a), LO))
}  else {
  s_One_a = One
  s_Yi_a = Y_a
  s_Li_a = LO
  Sig_invOne_a = matrix(0, ncol = 1, nrow = n)
  Sig_invY_a = matrix(0, ncol = 1, nrow = n)
  Sig_invLp_a = matrix(0, ncol = n, nrow = n)
  for (it in 1:M) {
    s_Onei_a = solve(L_a, solve(t(L_a), delta_a * s_Onei_a))
    Sig_invOne_a = Sig_invOne_a + s_Onei_a/delta_a
    s_Yi_a = solve(L_a, solve(t(L_a), delta_a * s_Yi_a))
    Sig_invY_a = Sig_invY_a + s_Yi_a/delta_a
    s_Li_a = solve(L_a, solve(t(L_a), delta_a * s_Li_a))
    Sig_invLp_a = Sig_invLp_a + s_Li/delta_a
  }
}

#########################################################
#Prediction functions (faster than predict.GP for optimization)
hpred = function(xnew){
    #htarget - desired thickness
    xn_h = t(matrix(xnew[hfact]))
    r = exp(-(abs(X_h - as.matrix(rep(1, n)) %*% (xn_h))^power) %*%
              (10^beta_h))
    yhat = (((1 - t(r) %*% Sig_invOne_h)/(t(One) %*% Sig_invOne_h)) %*%
              t(One) + t(r)) %*% Sig_invY_h
    return(yhat-htarget)
}

anglepred = function(xnew){
    #angletarget - desired angle
    xn_a = t(matrix(xnew[anglefact]))
    r = exp(-(abs(X_a - as.matrix(rep(1, n)) %*% (xn_a))^power) %*%
              (10^beta_a))
    yhat = (((1 - t(r) %*% Sig_invOne_a)/(t(One) %*% Sig_invOne_a)) %*%
              t(One) + t(r)) %*% Sig_invY_a
    return(yhat-angletarget)
}

obj = function(xnew){
  return(c(hpred(xnew),anglepred(xnew)))
}

#########################################################
#Check our predictions with predict.GP

# hfact = c(3,4,5) #sig. factors for h
# anglefact = c(2,3,4,5) # ... for angle
hfact = 1:5
anglefact = 1:5

xnew = rep(0.5,5) #desired prediction point

htarget = 0.5
angletarget = 40

hpred(xnew)
predict.GP(hfit,t(matrix(xnew[hfact])))$Y_hat-htarget
#same

anglepred(xnew)
predict.GP(anglefit,t(matrix(xnew[anglefact])))$Y_hat-angletarget
#also same

#one initial point
htarget = 0.75
angletarget = 60
onePt = multiStart(rep(0.5,5),obj,upper=rep(1,5),lower=rep(0,5))

onePt$par #the inverse point
apply(onePt$par,1,obj) #check == 0

#space-filling grid of initial points
htarget = 0.75
angletarget = 60
numPts = 25
multPts = multiStart(sobol(numPts,5),obj,upper=rep(1,5),lower=rep(0,5))

#get rid of anything outside of [0,1]^5
idx = rowSums((multPts$par[multPts$converged,]<0)+(multPts$par[multPts$converged,]>1))
multPts = multPts$par[multPts$converged,][-which(idx>0),]

apply(multPts,1,obj) #check == 0
scatter3d(multPts[,1],multPts[,2],multPts[,3])
