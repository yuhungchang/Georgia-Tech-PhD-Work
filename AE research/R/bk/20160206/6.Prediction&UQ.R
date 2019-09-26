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

####################################################################################
#RCpp code
####################################################################################
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
          List GPpred(NumericMatrix& xNewMat, NumericMatrix& xMat, NumericMatrix& Rcpp_yMat, NumericVector& corr, NumericVector& mu, NumericMatrix& Tcov, double pw) {
          //xNewMat - New settings to predict
          //xMat - Design matrix
          //yMat - POD coefficients from data
          //corr - Vector of correlation parameters
          //mu   - Mean of POD coef.
          //Tcov - Covariance matrix 
          //pw   - Power for exponential correlation
          
          int newn = xNewMat.nrow();
          int n = xMat.nrow();
          int p = xMat.ncol();
          int numCoef = Rcpp_yMat.ncol();
          arma::mat yMat(Rcpp_yMat.begin(),n,numCoef,false);
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
          
          ////Nugget
          //for (int i = 0; i < n; i++){
          //R(i,i) += 1e-6;
          //}
          
          //Update Rnew
          for (int i = 0; i < newn; i++){
            for (int j = 0; j < n; j++){
              for (int k = 0; k < p; k++){
                Rnew(i,j) = Rnew(i,j) * exp( - pow(10.0,corr(k)) * pow(xNewMat(i,k)-xMat(j,k), 2) );
              }
            }
          }     
          
          arma::mat Rinv = inv_sympd(R);

          //Compute prediction
          arma::mat pred(newn,numCoef);
          arma::vec onevec(n);
          arma::vec onenewvec(newn);
          onevec.fill(1.0);
          onenewvec.fill(1.0);
          arma::vec dif;
          arma::mat Rnewinv = Rnew * Rinv;
          for (int i=0; i<numCoef; i++){//For each POD coefficient...
            dif = yMat.col(i)-mu(i)*onevec;
            pred.col(i) = mu(i)*onenewvec + Rnewinv*dif;
          }

          //Compute UQ
          arma::mat var(newn,numCoef);
          arma::mat tmpmat = (Rnew * (Rinv * Rnew.t()));
          for (int i=0; i<numCoef; i++){
            for (int j=0; j<newn; j++){
              var(j,i) = Tcov(i,i) * ( 1 - tmpmat(j,j) );              
            }
          }
                  
          return (List::create(Named(\"Rnew\") = Rnew, Named(\"Rinv\") = Rinv, Named(\"pred\") = pred, Named(\"var\") = var));
          
          }')

####################################################################################
#Loading files and coordinates
####################################################################################

#Set working directory
filepath <- paste0(getwd(),"/")
#filepath <-"" # if filepath not execution directory
setwd(filepath)

base = 1 #Base case 
d = 5 #Number of control factors
indvec = 4:9
tvec = 1:999 #Times to predict

#Compute new design matrix
des = as.matrix(read.csv('desnorm.csv'));
desorig = read.csv('DoE.csv')[2:6]/1000
rangedes = rbind(apply(desorig,2,max),apply(desorig,2,min))
numPred = 1;
newCases = matrix(0,nrow=numPred,ncol=d)
Lp = desorig[base,1]*0.1/(rangedes[1,1]-rangedes[2,1])
thetap = desorig[base,3]*0.1/(rangedes[1,3]-rangedes[2,3])
deltap = desorig[base,4]*0.1/(rangedes[1,4]-rangedes[2,4])
newCases[1,] = as.matrix(des[base,] + Lp*diag(d)[1,] + thetap*diag(d)[3,] + deltap*diag(d)[4,])

#Compute new coordinates
# load('CommonCoord.RData')
# Coordinates = vector('list',numPred) #One for each new case
# oneCase = 1
# for (i in 1:numPred){
#   # New dimensions
#   newL = (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1])
#   newR = (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])
#   newdL = (newCases[i,5]*(rangedes[1,5]-rangedes[2,5]) + rangedes[2,5])
#   # Injector start: scale by dL on x, scale by R on y
#   Coordinates[[i]] = cbind(CommonCoord[[oneCase]][[1]][,1]/desorig[oneCase,5]*(newdL),
#                            CommonCoord[[oneCase]][[1]][,2]/desorig[oneCase,2]*(newR))
#   # Injector body: scale by (L-dL) on x, scale by R on y
#   Coordinates[[i]] = rbind(Coordinates[[i]],cbind((CommonCoord[[oneCase]][[2]][,1]-desorig[oneCase,5])/(desorig[oneCase,1]-desorig[oneCase,5])*(newL- newdL) + newdL,
#                                                   CommonCoord[[oneCase]][[2]][,2]/desorig[oneCase,2]*(newR)) )
#   # Bottom chamber: 
#   Coordinates[[i]] = rbind(Coordinates[[i]],cbind( ( CommonCoord[[oneCase]][[3]][,1] - desorig[oneCase,1] ) * (newR / desorig[oneCase,2]) + newL,
#                                                   CommonCoord[[oneCase]][[3]][,2] * (newR / desorig[oneCase,2]) ) )
#   # Top chamber:
#   Coordinates[[i]] = rbind(Coordinates[[i]],cbind( ( CommonCoord[[oneCase]][[4]][,1] - desorig[oneCase,1] ) * (newR / desorig[oneCase,2]) + newL ,
#                                                   ( CommonCoord[[oneCase]][[4]][,2] - desorig[oneCase,2] ) * (newR / desorig[oneCase,2]) + newR) )
# }

load('CommonCoord.RData')
oneCase = 22
blInj1 = CommonCoord[[oneCase]][[1]]
blInj2 = CommonCoord[[oneCase]][[2]]
blBC = CommonCoord[[oneCase]][[3]]
blTC = CommonCoord[[oneCase]][[4]]
Lval = desorig[,1]
Rval = desorig[,2]
DLval = desorig[,5]

Coordinates = vector('list', numPred) #One for each new case
for (i in 1:numPred){
  
  newL = (newCases[i,1]*(rangedes[1,1]-rangedes[2,1]) + rangedes[2,1])
  newR = (newCases[i,2]*(rangedes[1,2]-rangedes[2,2]) + rangedes[2,2])
  newDL = (newCases[i,5]*(rangedes[1,5]-rangedes[2,5]) + rangedes[2,5])
  
  tmpList = vector('list',4)
  
  #First part: part one of inside injector (aspect ratio = 10)
  tmpList[[1]] = cbind(newDL*blInj1[,1]/DLval[oneCase], newR*blInj1[,2]/Rval[oneCase])
  
  #Second part: part two of inside injector (aspect ratio = 10)
  tmpList[[2]] = cbind(newDL + (newL - newDL)* (blInj2[,1] - DLval[oneCase])/(Lval[oneCase] - DLval[oneCase]), newR*blInj2[,2]/Rval[oneCase])
  
  #Third part: bottom of chamber (aspect ratio = 10)
  tmpList[[3]] = cbind(newL + (blBC[,1]-Lval[oneCase]) * newR/Rval[oneCase], newR*blBC[,2]/Rval[oneCase])
  
  #Fourth part: top of chamber (aspect ratio = cutX/cutY)
  tmpList[[4]] = cbind(newL + (blTC[,1]-Lval[oneCase]) * newR/Rval[oneCase], newR + (blTC[,2]-Rval[oneCase])* newR/Rval[oneCase])
  
  #   numFinalGrid = nrow(grid1) + nrow(grid2) + nrow(grid3)
  Coordinates[[i]] = rbind(tmpList[[1]], tmpList[[2]], tmpList[[3]], tmpList[[4]])
}

#Set colors of plot
numCol <- 1000 #numCol colors in legend
colors <- matlab.like(numCol+1)

#Loading number of modes and POD data
numModes = rep(-1,length(indvec)) #Number of modes to choose from
BTlist = vector('list',length(indvec))
PsiTlist = vector('list',length(indvec))
FTmeanlist = vector('list',length(indvec))
cut = 0.95

for (i in 1:length(indvec)){
  #Load data
  ind = indvec[i] #Current flow response
  setwd(paste0(filepath,"bigmemory_bk/", ind))
  load('PODdat.RData')
  load('FTmean.vector.RData')
  
  BTlist[[i]] = BTsvg
  PsiTlist[[i]] = PsiTsvg
  FTmeanlist[[i]] = FTmean.vector
  
  #Compute cut-off
  cs = cumsum(eigval)/sum(eigval)
  cs = pmax(cs-cut,0)
  cs[cs==0] = Inf
  #numModes[i] = which.min(cs)
  numModes[i] = ncol(PsiTsvg)
}
numCoef = sum(numModes)


####################################################################################
#Computing prediction and UQ
####################################################################################

# Set working directory
setwd(filepath)

na.index <- c(192, 385, 578, 771, 964)
numPred = nrow(newCases)

for (tt in 1:length(tvec)){
  
  print(paste0('Predict/UQ for t = ', tt, '...'))
  
  #Load fitted parameters
  load(paste0('kriging_T/allparams_t',tvec[tt],'.RData'))
  
  #Compute prediction and UQ at new settings
#   newCases <- matrix(des[1,],nrow=1)
  predcoef <- GPpred(newCases, des, params$yMat, params$corr, params$mu, as.matrix(params$Tcov), 2)  
  coefp <- predcoef$pred
  varp <- predcoef$var
#   coefp-params$yMat[1,]
  
  #Reconstruct prediction for each flow
  cModes = c(0,cumsum(numModes))
  predsnap = vector('list',length(indvec))
  sdsnap = vector('list',length(indvec))
  for (i in 1:length(indvec)){
    predsnap[[i]] <- FTmeanlist[[i]] + PsiTlist[[i]][,1:numModes[i]]%*%(coefp[,(cModes[i]+1):cModes[i+1]])
    sdsnap[[i]] <- 1.96 * sqrt( (PsiTlist[[i]][,1:numModes[i]]^2)%*%(varp[,(cModes[i]+1):cModes[i+1]]) / 30)
  }
  save(predcoef,predsnap,sdsnap,file=paste0('kriging_T/pred/pred_',tvec[tt],'.RData'))
  
  #Plot prediction and UQ for all flows
  for (i in 1:length(indvec)){
    
    #Find range
    if (indvec[i]==4){
      Range = c(-50, 50)
      sdRange = c(0, 5)
    }else if(indvec[i]==5){
      Range = c(-30, 30)
      sdRange = c(0, 2)
    }else if(indvec[i]==6){
      Range = c(-30, 30)
      sdRange = c(0, 3)
    }else if(indvec[i]==7){
      Range = c(-400,400)
      sdRange = c(0, 60)
    }else if(indvec[i]==8){
      Range = c(100, 350)
      sdRange = c(0, 20)      
    }else if(indvec[i]==9){
      Range = c(100, 1200)
      sdRange = c(0, 100)
    }    
    
    for (j in 1:numPred){
      
      #Compute color and output prediction
      zcolor <- colors[ pmin(pmax( (predsnap[[i]][,j] - Range[1])/diff(Range)*numCol, 0) + 1,numCol) ]
      png(paste0('kriging_T/image/pred-c',j,'_ind',indvec[i],'_t',tvec[tt],'.png'), width = 600, height = 400)
      plot(Coordinates[[j]][-na.index,1],Coordinates[[j]][-na.index,2],col= zcolor,pch=15,cex=0.5,
           xlab = 'x (m)', ylab = 'y (m)', main = paste0('Prediction at T = ', tvec[tt],' ms'),
           xlim = c(0,0.04), ylim = c(0,0.008) ) 
      colorlegend(col=colors,zlim=Range,left=T)
      dev.off()
      
      #Compute color and output UQ
      zcolor <- colors[ pmin(pmax( (sdsnap[[i]][,j] - sdRange[1])/diff(sdRange)*numCol, 0) + 1,numCol) ]
      png(paste0('kriging_T/image/uq-c',j,'_ind',indvec[i],'_t',tvec[tt],'.png'), width = 600, height = 400)
      plot(Coordinates[[j]][-na.index,1],Coordinates[[j]][-na.index,2],col= zcolor,pch=15,cex=0.5,
           xlab = 'x (m)', ylab = 'y (m)', main = paste0('UQ at T = ', tvec[tt],' ms'),
           xlim = c(0,0.04), ylim = c(0,0.008) ) 
      colorlegend(col=colors,zlim=sdRange,left=T)
      dev.off()
    }
    
  }
  
  #Export as .dat file:
  #Load in prediction and UQ
  predMat <- matrix(0,nrow=length(predsnap[[1]]),ncol=2+length(indvec))
  uqMat <- matrix(0,nrow=length(predsnap[[1]]),ncol=2+length(indvec))
  predMat[,1:2] = Coordinates[[1]][-na.index,]
  uqMat[,1:2] = Coordinates[[1]][-na.index,]
  for (i in 1:length(indvec)){
    predMat[,i+2] = predsnap[[i]]
    uqMat[,i+2] = sdsnap[[i]]
  }
  colnames(predMat) = c('x','y','u','v','w','P','T','rho')
  colnames(uqMat) = c('x','y','u','v','w','P','T','rho')
  
  #Write predictions and UQ
  write.table(predMat, file=paste0('kriging_T/pred/pred_',tvec[tt],'.dat'), row.names = F)
  write.table(uqMat, file=paste0('kriging_T/pred/uq_',tvec[tt],'.dat'), row.names = F)

}
