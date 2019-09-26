library(Rcpp)
library(bigmemory)
library(BH)
library(calibrator)
library(bigpca)
library(rARPACK)
library(parallel)
library(colorRamps)
library(shape)
library(nloptr)

#RCpp code
sourceCpp(code='
          
          #include <iostream>
          #include <omp.h>
          
          //  Objective function
          //
          // [[Rcpp::depends(RcppArmadillo)]]
          #include <RcppArmadillo.h>
          #include <iostream>
          
          using namespace Rcpp;
          
          // [[Rcpp::plugins(openmp)]]
          // [[Rcpp::export]]
          double fastGPdev(NumericVector& corr, NumericMatrix& xMat, NumericVector& Rcpp_yVec, double pw) {
          //xMat - Design matrix
          //yVec - Observed response
          //pw   - Power for exponential correlation
          
          int n = xMat.nrow();
          int p = xMat.ncol();
          arma::vec yVec(Rcpp_yVec.begin(),n,false);
          arma::mat R(n,n); //Correlation matrix
          R.fill(1.0); //Fill with ones
          
          //Update R
          for (int i = 0; i < n; i++){
          for (int j = (i+1); j < n; j++){
          for (int k = 0; k < p; k++){
          R(i,j) = R(i,j) * exp( - pow(10.0,corr(k)) * pow(xMat(i,k)-xMat(j,k), 2) );
          //std::cout << - pow(10.0,corr(k)) << std::endl;
          //std::cout << tmpdbl << std::endl;
          //std::cout << pow(xMat(i,k)-xMat(j,k), pw) ) << std::endl;
          //std::cout << "(i,j): " << i << ", " << j << "  " << R(i,j) << std::endl;
          }
          R(j,i) = R(i,j); //Symmetrize
          }
          }

          //Nugget
          for (int i = 0; i < n; i++){
          R(i,i) += R(i,i) + 1e-6;
          }
          
          arma::mat Rinv = inv_sympd(R);
          
          double mu = (1/accu(Rinv))*accu(Rinv*yVec);
          arma::vec onevec(n);
          onevec.fill(1.0);
          arma::vec dif = yVec-mu*onevec;
          double sigma2 = dot(dif, ( Rinv * dif ) );
          double ret = (double)n*log(sigma2) + log(det(R));
          return(ret);
          
          }')

sourceCpp(code='
          
          #include <iostream>
          #include <omp.h>
          
          //  Objective function
          //
          // [[Rcpp::depends(RcppArmadillo)]]
          #include <RcppArmadillo.h>
          #include <iostream>
          
          using namespace Rcpp;
          
          // [[Rcpp::plugins(openmp)]]
          // [[Rcpp::export]]
          List getParam(NumericVector& corr, NumericMatrix& xMat, NumericVector& Rcpp_yVec, double pw) {
          //xMat - Design matrix
          //yVec - Observed response
          //pw   - Power for exponential correlation
          
          int n = xMat.nrow();
          int p = xMat.ncol();
          arma::vec yVec(Rcpp_yVec.begin(),n,false);
          arma::mat R(n,n); //Correlation matrix
          R.fill(1.0); //Fill with ones
          
          //Update R
          for (int i = 0; i < n; i++){
          for (int j = (i+1); j < n; j++){
          for (int k = 0; k < p; k++){
          R(i,j) = R(i,j) * exp( - pow(10.0,corr(k)) * pow(xMat(i,k)-xMat(j,k), 2) );
          //std::cout << - pow(10.0,corr(k)) << std::endl;
          //std::cout << tmpdbl << std::endl;
          //std::cout << pow(xMat(i,k)-xMat(j,k), pw) ) << std::endl;
          //std::cout << "(i,j): " << i << ", " << j << "  " << R(i,j) << std::endl;
          }
          R(j,i) = R(i,j); //Symmetrize
          }
          }
          
          //Nugget
          for (int i = 0; i < n; i++){
            R(i,i) += R(i,i) + 1e-6;
          }

          arma::mat Rinv = inv_sympd(R);
          
          double mu = (1/accu(Rinv))*accu(Rinv*yVec);
          arma::vec onevec(n);
          onevec.fill(1.0);
          arma::vec dif = yVec-mu*onevec;
          double sigma2 = dot(dif, ( Rinv * dif ) );
          //double ret = (double)n*log(sigma2) + log(det(R));
          return (List::create(Named("corr") = corr, Named("mu") = mu, Named("sigma2") = sigma2));
          
          }')


#RCpp code
sourceCpp(code='
          
          #include <iostream>
          #include <omp.h>
          
          //  Objective function
          //
          // [[Rcpp::depends(RcppArmadillo)]]
          #include <RcppArmadillo.h>
          #include <iostream>
          
          using namespace Rcpp;
          
          // [[Rcpp::plugins(openmp)]]
          // [[Rcpp::export]]
          List fastGPpred(NumericMatrix& xNewMat, NumericVector& corr, NumericMatrix& xMat, NumericVector& Rcpp_yVec, double pw) {
          //xNewMat - New settings to predict
          //xMat - Design matrix
          //yVec - Observed response
          //pw   - Power for exponential correlation
          
          int newn = xNewMat.nrow();
          int n = xMat.nrow();
          int p = xMat.ncol();
          arma::vec yVec(Rcpp_yVec.begin(),n,false);
          arma::mat R(n,n); //Correlation matrix
          R.fill(1.0); //Fill with ones
          arma::mat Rnew(newn,n); //New Correlation matrix
          Rnew.fill(1.0); //Fill with ones
          
          //Update R
          for (int i = 0; i < n; i++){
          for (int j = (i+1); j < n; j++){
          for (int k = 0; k < p; k++){
          R(i,j) = R(i,j) * exp( - pow(10.0,corr(k)) * pow(xMat(i,k)-xMat(j,k), 2) );
          }
          R(j,i) = R(i,j); //Symmetrize
          }
          }
          
          //Nugget
          for (int i = 0; i < n; i++){
          R(i,i) += R(i,i) + 1e-6;
          }
          
          //Update Rnew
          for (int i = 0; i < newn; i++){
          for (int j = 0; j < n; j++){
          for (int k = 0; k < p; k++){
          Rnew(i,j) = Rnew(i,j) * exp( - pow(10.0,corr(k)) * pow(xNewMat(i,k)-xMat(j,k), 2) );
          }
          }
          }     
          
          arma::mat Rinv = inv_sympd(R);
          
          double mu = (1/accu(Rinv))*accu(Rinv*yVec);
          arma::vec onevec(n);
          arma::vec onenewvec(newn);
          onevec.fill(1.0);
          arma::vec dif = yVec-mu*onevec;
          double sigma2 = dot(dif, ( Rinv * dif ) );
          
          arma::vec ret = Rnew * (Rinv * dif);
          ret = ret + mu*onenewvec;
          
          //           std::cout << Rnew.n_rows << ", " << Rnew.n_cols << std::endl;
          //           std::cout << Rinv.n_rows << ", " << Rinv.n_cols << std::endl;
          arma::mat tmpmat = (Rnew * (Rinv * Rnew.t()));
          
          arma::vec var(newn);
          for (int i = 0; i < newn; i++){
          var(i) = sigma2 * ( 1 - tmpmat(i,i) );
          }
          
          return (List::create(Named("pred") = ret, Named("var") = var));
          
          }')


# #Test for fastGPdev
# fastGPdev(rep(0.1,5),des,curBT[1,],2)
# GP_dev(rep(0.1,5),des,curBT[1,],corr = list(type="exponential",power=2))
# 
# #Optimize
# GPfitsol <- GP_fit(des,curBT[1,])
# fastGPdev(GPfitsol$beta,des,curBT[1,],2)
# mySol <- lbfgs(x0=rep(0.5,5), fn=fastGPdev, xMat=des, Rcpp_yVec=curBT[1,], pw=2)
# fastGPdev(mySol$par,des,curBT[1,],2)
# getParam(mySol$par,des,curBT[1,],2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting the kriging model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set working directory
ind = 5
setwd(paste0("/home/proj/jeffwu/isye-ae/bigmemory_bk/", ind))

#Load common coordinate indices and data
load("../../CommonCoord.RData")
load("FTmean.vector.RData")
load("PODdat.RData")
na.index <- c(385,578)
des = as.matrix(read.csv("../../desnorm.csv"));

#Get the POD modes for one timestep
BT <- t(BTsvg)
# tme <- 1 #Current timestep
# GP_fit(des,curBT[1,],control=c(2,2,1),maxit=10,corr=list(type="exponential",power=2))

#Now do this for all timesteps
numModes = 500;

ptm <- proc.time();
for (t in 501:999){
# for (t in 67:(ncol(BT)/30)){
  
  t1 <- Sys.time()
  
  print(paste0("Fitting - ind: ",ind," time: ", t))
  idx <- seq(from=t,to=ncol(BT),by=ncol(BT)/30)
  curBT <- BT[,idx]
  params = vector("list",numModes)
  for (i in 1:numModes){
    sv <- lbfgs(rep(0.5,5), fn=fastGPdev, control = list(xtol_rel = 1e-1), xMat=des, Rcpp_yVec=curBT[i,], pw=2)
    params[[i]] = getParam(sv$par,des,curBT[i,],2)
  }
  
  proc.time() - ptm; #1.771 sec for each
  
  save(params,file=paste0("../../kriging_indep/",ind,"/params_resp",ind,"_time",t,".RData"))
  
  t2 <- Sys.time()
  print(t2-t1)
  
}

proc.time() - ptm; #1.771 sec for each

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prediction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Set working directory
ind = 5
setwd(paste0("/home/proj/jeffwu/isye-ae/bigmemory_bk/", ind))
base = 2 #Base case 
d = 5 #Number of control factors
#Compute new design matrix
des = as.matrix(read.csv("../../desnorm.csv"));
newCases = matrix(0,nrow=d,ncol=d)
newCases[1,] = as.matrix(des[base,] + 0.2*diag(d)[1,])
newCases[2,] = as.matrix(des[base,] + 0.2*diag(d)[2,])
newCases[3,] = as.matrix(des[base,] + 0.2*diag(d)[3,])
newCases[4,] = as.matrix(des[base,] + 0.2*diag(d)[4,])
newCases[5,] = as.matrix(des[base,] + 0.2*diag(d)[5,])
desorig = read.csv("../../30mppts_by_week.csv")[2:6]/1000
rangedes = rbind(apply(desorig,2,max),apply(desorig,2,min))

#Compute new coordinates
load("../../CommonCoord.RData")
load("FTmean.vector.RData")
load("PODdat.RData")
Coordinates = vector("list",d) #One for each new case
oneCase = 2
for (i in 1:d){
  # New dimensions
  newL = (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1])
  newR = (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])
  newdL = (newCases[i,5]*(rangedes[1,5]-rangedes[2,5]) + rangedes[2,5])
  # Injector start: scale by dL on x, scale by R on y
  Coordinates[[i]] = cbind(CommonCoord[[oneCase]][[1]][,1]/desorig[oneCase,5]*(newdL),
                           CommonCoord[[oneCase]][[1]][,2]/desorig[oneCase,2]*(newR))
  # Injector body: scale by (L-dL) on x, scale by R on y
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind((CommonCoord[[oneCase]][[2]][,1]-desorig[oneCase,5])/(desorig[oneCase,1]-desorig[oneCase,5])*(newL- newdL) + newdL,
                                                  CommonCoord[[oneCase]][[2]][,2]/desorig[oneCase,2]*(newR)) )
  # Bottom chamber: 
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind(CommonCoord[[oneCase]][[3]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                                                  CommonCoord[[oneCase]][[3]][,2]/desorig[oneCase,2]*(newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])))
  # Top chamber:
  Coordinates[[i]] = rbind(Coordinates[[i]],cbind(CommonCoord[[oneCase]][[4]][,1] - desorig[oneCase,1] + (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1]),
                                                  CommonCoord[[oneCase]][[4]][,2] - desorig[oneCase,2] + (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])))
}

#set range
numCol <- 1000 #numCol colors in legend
colors <- matlab.like(numCol+1)
if (ind==4){
  Range = c(-101, 69)
  uqRange = c(0,20)
}else if(ind==5){
  Range = c(-46, 56)
  uqRange = c(0,7)
}else if(ind==6){
  Range = c(-65, 12)
  uqRange = c(0,40)
}else if(ind==7){
  Range = c(-500,1180)
  uqRange = c(0,2000)
}else if(ind==8){
  Range = c(120, 310)
  uqRange = c(0,50)
}else if(ind==9){
  Range = c(125, 1010)
  uqRange = c(0,20)
}else if(ind==13){
  Range = c(-400000,1000000)
  uqRange = c(0,20)
}
na.index <- c(385,578)

#Do predictions
numModes = 500
BT <- t(BTsvg)
# for (t in 1:500){
for (t in 130:500){
  
  #Compute POD coefficient prediction
  coefpred <- matrix(0,nrow=numModes,ncol=d)
  varpred <- matrix(0,nrow=numModes,ncol=d)
  print(paste0("Predicting - ind: ",ind," time: ", t))
  idx <- seq(from=t,to=ncol(BT),by=ncol(BT)/30)
  curBT <- BT[,idx]
  load(paste0("../../kriging_indep/",ind,"/params_resp",ind,"_time",t,".RData"))
  
  for (i in 1:numModes){
    tmpobj <- fastGPpred(newCases,params[[i]]$corr,des,curBT[i,],2)
    coefpred[i,] <- tmpobj$pred
    varpred[i,] <- tmpobj$var/30
  }
  save(coefpred,file=paste0("../../kriging_indep/",ind,"/coef_resp",ind,"_time",t,".RData"))
  
  #Reconstruct snapshot prediction
  predsnap <- PsiTsvg%*%coefpred
  sdsnap <- 1.96*sqrt((PsiTsvg^2)%*%varpred) #95% confidence
  predsnap <- FTmean.vector+predsnap
  save(predsnap,sdsnap,file=paste0("../../kriging_indep/",ind,"/pred_resp",ind,"_time",t,".RData"))
  
  #Plot image
  for (j in 1:5){
    
    #Prediction
#     zcolor <- pmin(pmax((predsnap[,j] - Range[1])/(Range[2]-Range[1])*numCol,0),numCol) + 1
#     png(paste0("../../kriging_indep/",ind,"/flow",j,"_resp",ind,"_time",t,".png"), width = 600, height = 400)
#     plot(Coordinates[[j]][-na.index,1],Coordinates[[j]][-na.index,2],col=colors[zcolor],pch=15,cex=1.2,
#          xlab = 'x (m)', ylab = 'y (m)', main = paste0("Flow prediction at T = ", t),
#          xlim = c(0,0.10), ylim = c(0,0.013) )
#     colorlegend(col=colors,zlim=Range,left=T)
#     dev.off()
    
    #UQ
    zcolor <- pmin(pmax((sdsnap[,j] - uqRange[1])/(uqRange[2]-uqRange[1])*numCol,0),numCol) + 1
    png(paste0("../../kriging_indep/",ind,"/UQ/UQ_flow",j,"_resp",ind,"_time",t,".png"), width = 600, height = 400)
    plot(Coordinates[[j]][-na.index,1],Coordinates[[j]][-na.index,2],col=colors[zcolor],pch=15,cex=1.2,
         xlab = 'x (m)', ylab = 'y (m)', main = paste0("Flow prediction at T = ", t),
         xlim = c(0,0.10), ylim = c(0,0.013) )
    colorlegend(col=colors,zlim=uqRange,left=T)
    dev.off()
  }
}

