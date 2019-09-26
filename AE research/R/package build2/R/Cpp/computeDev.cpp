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
}