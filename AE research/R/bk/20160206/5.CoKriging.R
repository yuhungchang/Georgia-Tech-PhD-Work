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
library(glasso)
library(huge)
library(SLHD)
library(MBESS)

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
          List computeS(NumericVector& corr, NumericMatrix& xMat, NumericMatrix& Rcpp_yMat, double pw) {
          //computeS returns covariance matrix S to use in glasso
          //corr - initial correlation
          //xMat - Design matrix (n x p)
          //Rcpp_yMat - Observed response matrix (n x (numModes.numResponse))
          //pw   - Power for exponential correlation
          
          int n = xMat.nrow();
          int p = xMat.ncol();
          int numCoef = Rcpp_yMat.ncol();
          arma::mat yMat(Rcpp_yMat.begin(),n,numCoef,false);
          arma::mat R(n,n); //Correlation matrix
          R.fill(1.0); //Fill with ones
          
          //Update R
          for (int i = 0; i < n; i++){
            for (int j = (i+1); j < n; j++){
              for (int k = 0; k < p; k++){
                R(i,j) = R(i,j) * exp( - pow(10.0,corr(k)) * pow(xMat(i,k)-xMat(j,k), 2.0) );
              }
            R(j,i) = R(i,j); //Symmetrize
            }
          }

          ////Nugget
          //for (int i = 0; i < n; i++){
          //  R(i,i) += 1e-6;
          //}
          
          arma::mat Rinv = inv_sympd(R);
                  
          //Compute difference matrix
          arma::mat mumat = (1/accu(Rinv))*(Rinv*yMat);
          arma::vec mu = sum(mumat.t(), 1);
          arma::mat diff(n,numCoef); //B - 1 otimes mu
          for (int i=0;i<n;i++){
            for (int j=0;j<numCoef;j++){
              diff(i,j) = mu(j);
            }
          }
          diff = yMat - diff;
          
          //Compute R^(1/2)%*%diff and call .glasso from Fortran
          arma::vec eigval;
          arma::mat eigvec;
          eig_sym(eigval, eigvec, Rinv);

          arma::mat Rsqrt = eigvec*( diagmat(sqrt(eigval))*eigvec.t() );
          arma::mat Rsd = (Rsqrt * diff);   
          
          return (List::create(Named(\"mu\") = mu, Named(\"Rinv\") = Rinv, Named(\"R\") = R, Named(\"Rsd\") = Rsd));

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
          using namespace std;
          
          // [[Rcpp::plugins(openmp)]]
          // [[Rcpp::export]]
          double computeDev(NumericVector& mu, NumericMatrix& Rcpp_Rinv, NumericVector& corr, NumericMatrix& Rcpp_yMat, NumericMatrix& Rcpp_Tinv, bool indep_ind) {
          //computeDev returns the deviance of the GP model
          //corr - initial correlation
          //mu - Computed mu
          //Rinv - Computed inverse of R
          //Rcpp_yMat - Observed response matrix (n x (numModes.numResponse))
          //Tinv - Computed inverse of T through glasso
          
          int n = Rcpp_Rinv.nrow();
          int numCoef = Rcpp_yMat.ncol();
          arma::rowvec mua(mu.begin(),numCoef,false);
          arma::mat yMat(Rcpp_yMat.begin(),n,numCoef,false);
          arma::mat Rinv(Rcpp_Rinv.begin(),n,n,false);
          arma::mat Tinv(Rcpp_Tinv.begin(),numCoef,numCoef,false);
          double ret = 0.0;

          //Precompute T^{-1}(B-mu), if computing T
          if (!indep_ind){
            arma::mat mumat = repmat(mua,n,1);
            arma::mat diff = (yMat - mumat);
            arma::mat TB = Tinv * diff.t();
            
            arma::rowvec tmpi(numCoef); 
            arma::vec tmpj(numCoef); 
            for (int i=0;i<n;i++){
              for (int j=0;j<n;j++){
                tmpi = diff.row(i);
                tmpj = TB.col(j);
                ret += Rinv(i,j) * dot( tmpi , tmpj);
              }
            }
          }

          //Log-det is easy if Tinv is diagonal, otherwise do log_det
          double val = 0.0;
          double sign = 1.0;
          if (!indep_ind){
            log_det(val,sign,Tinv);
          }
          else{
            for (int k=0; k<numCoef; k++){
              val += log(Tinv(k,k));
            }
          }
          ret += n * ( log(n)*numCoef - sign*val ) - numCoef*log(det(Rinv));
          //ret += (- n * sign*val ) + numCoef*log(det(Rinv));

          //std::cout << \"final: \" << ret << std::endl;
 
          return(ret);
          }')
  
GPdev <- function(corr, xMat, yMat, pw=2){
#   corr <- allPar[which.min(allVal),]
  
  #See sourceCpp function for parameter details
  n <- nrow(xMat)
    
  #Update parameters
  params <- computeS(corr, xMat, yMat, pw);   #params[1] - mu, params[2] - Rinv, params[3] - data matrix for T

  #Tune Tinv
  glhuge = huge(x = params$Rsd, method='glasso') #Screening for speed-up
  glsel = huge.select(glhuge)
  optind = which.min(glsel$ebic.score[-1])+1 #Find the best model which gives correlations
  vars = apply(params$Rsd,2,function(xx){mean(xx^2)})*30/29  #Adjust by variance
  Tinv = sweep(glsel$icov[[optind]],MARGIN=2,vars*diag(glsel$icov[[optind]]),'/')
  Tinv = as.matrix(Tinv)
  
  #Compute deviance
  ret = computeDev(c(params$mu), params$Rinv, corr, yMat, Tinv, 0) #0 indicates not indep
  print(ret)
  return(ret)
}

GPdev.ind <- function(corr, xMat, yMat, pw=2){
  #   corr <- allPar[which.min(allVal),]
  
  #See sourceCpp function for parameter details
  n <- nrow(xMat)
  
  #Update parameters
  params <- computeS(corr, xMat, yMat, pw);   #params[1] - mu, params[2] - Rinv, params[3] - data matrix for T
  
  #Testing as independent kriging
  vars = apply(params$Rsd,2,function(xx){mean(xx^2)}) #Adjust by variance
  if (length(vars)==1){
    Tinv = matrix(1/vars, ncol=1, nrow=1)
  }else{
    Tinv = diag(1/vars)
  }
  
  #Compute deviance
#   dim(yMat)
  ret = computeDev(c(params$mu), params$Rinv, corr, yMat, Tinv, 1) #1 indicates indep
#   print(ret)
  return(ret)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Computing cut-off for POD modes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set working directory
filepath <- paste0(getwd(),"/")
#filepath <- ""              #set filepath if directory other than execution directory
#setwd(filepath)

indvec = 4:9
tvec = 1:999 #Timesteps to tune

#Load design
des = as.matrix(read.csv('desnorm.csv'));

#Set cutoff for POD modes at 95% energy
numModes = rep(-1,length(indvec)) #Number of modes to choose from
BTlist = vector('list',length(indvec))
cut = 0.95

for (i in 1:length(indvec)){
  #Load data
  ind = indvec[i] #Current flow response
  setwd(paste0(filepath,"bigmemory_bk/", ind))
  load('PODdat.RData')
  BTlist[[i]] = BTsvg
  
  #Compute cut-off
  cs = cumsum(eigval)/sum(eigval)
  cs = pmax(cs-cut,0)
  cs[cs==0] = Inf
#   numModes[i] = which.min(cs)
  numModes[i] = ncol(PsiTsvg)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tuning correlation parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(filepath)
pars <- matrix(0,nrow=length(tvec),ncol=ncol(des))

#Tuning correlation parameters
for (t in 1:length(tvec)){
  print(paste0('Tuning corr. param. for t = ', t, '...'))
  
  yMat = matrix(0,nrow=30,ncol=0) #Matrix of POD coefficients (n x (numModes.numResponse))
  slc = seq(from=tvec[t],by=999,length.out=30)
  
  for (i in 1:length(indvec)){
    #Find indices for each case
    curNumMode = numModes[i]
    yMat = cbind( yMat, BTlist[[i]][slc,1:curNumMode] ) 
  }
  
#   load(paste0('yMat_',t,'.RData'))
 
  #If t == 1, do a global search for correlation parameters
  if (t == 1){
    lower = -5
    upper = 5
    numIni = 100
    Mm = (maximinSLHD(1,numIni,5)$StandDesign)*(upper-lower)+lower
    svList = vector('list',numIni)
    allVal = rep(-1,numIni)
    allPar = matrix(-1,nrow=numIni,ncol=5)
    for (i in 1:numIni){
      print(i)
      svList[[i]] <- lbfgs(Mm[i,], fn=GPdev.ind, control = list(xtol_rel = 1e-2), 
                           lower = rep(lower,5), upper = rep(upper,5), xMat=des, yMat=yMat, pw=2)
      allVal[i] = svList[[i]]$value
      allPar[i,] = svList[[i]]$par
    }
    pars[t,] = allPar[which.min(allVal),]
  }
  else{
    pars[t,] <- lbfgs(pars[(t-1),], fn=GPdev.ind, control = list(xtol_rel = 1e-2), 
          lower = rep(lower,5), upper = rep(upper,5), xMat=des, yMat=yMat, pw=2)$par
  }
}

save(pars,file='kriging_T/corrparams.RData')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tuning correlation matrix T (parallelized)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tvec = 1:999 #Timesteps to tune

setwd(filepath)
load('kriging_T/corrparams.RData')

for (t in 1:length(tvec)){
  
  print(paste0('Tuning T for t = ', t, '...'))

  yMat = matrix(0,nrow=30,ncol=0) #Matrix of POD coefficients (n x (numModes.numResponse))
  slc = seq(from=tvec[t],by=999,length.out=30)
  
  for (i in 1:length(indvec)){
    #Find indices for each case
    curNumMode = numModes[i]
    yMat = cbind( yMat, BTlist[[i]][slc,1:curNumMode] ) 
  }
  
  # #  Actual tuning through profile likelihood
  #   ptm <- proc.time()
  #   corr <- lbfgs(allPar[which.min(allVal),], fn=GPdev, control = list(xtol_rel = 1e-1), 
  #                        lower = rep(-10,5), upper = rep(2,5), xMat=des, yMat=yMat, pw=2)$par
  #   proc.time() - ptm
  
  # Plug-in (GLASSO)
  corr <- pars[t,]
  pm <- computeS(corr, des, yMat, 2);   #params[1] - mu, params[2] - Rinv, params[3] - data matrix for T
  glhuge = huge(x = pm$Rsd, method = 'glasso', lambda.min.ratio=0.1)
  glsel = huge.select(glhuge)
  optind = which.min(glsel$ebic.score[-1])+1 #Find the best model which gives correlations
  vars = apply(pm$Rsd,2,function(xx){mean(xx^2)}) #Adjust by variance
  invMat = cov2cor(as.matrix(glsel$icov[[optind]]))
  Tinv = cor2cov( (invMat + t(invMat))/2, 1/sqrt(vars*diag(glsel$icov[[optind]])) )
  Tcov = solve(Tinv)
  
#   # Plug-in (CT)
#   corr <- pars[t,]
#   pm <- computeS(corr, des, yMat, 2);   #params[1] - mu, params[2] - Rinv, params[3] - data matrix for T
#   glhuge = huge(x = pm$Rsd, method='ct', lambda.min.ratio = 0.5)
#   glsel = huge.select(glhuge)
#   Tinv = glsel$path[[glsel$opt.index]]
#   optind = which.min(glsel$ebic.score[-1])+1 #Find the best model which gives correlations
#   vars = apply(pm$Rsd,2,function(xx){mean(xx^2)})*29/30 #Adjust by variance
#   Tinv = as.matrix(sweep(glsel$icov[[optind]],MARGIN=2,vars*diag(glsel$icov[[optind]]),'/'))
#   Tcov = as.matrix(sweep(glsel$cov[[optind]],MARGIN=2,vars/diag(glsel$cov[[optind]]),'*'))
  
  #Find which correlations are significant
  testmat = (Tinv)
  idx = which(testmat!=0,arr.ind=T)  
  offdiag = apply(idx,1,function(ind){
    if (ind[1]==ind[2]){
      return (0)
    }
    else{
      return (1)
    }
  })
  idxnz = idx[which(offdiag==1),]
#     for (k in 1:nrow(idxnz)){
#       print(paste0(idxnz[k,1],', ',idxnz[k,2],' - Covariance: ',testmat[idxnz[k,1],idxnz[k,2]]))
#     }
  
  #Save to cluster
  params = list(mu=pm$mu, Ticov=Tinv, Tcov=Tcov, corr=corr, idxnz=idxnz, yMat=yMat, numModes=numModes)
  save(params, file=paste0('kriging_T/allparams_t',tvec[t],'.RData'))
}
