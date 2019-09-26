#include <iostream>
#include <omp.h>

//  Objective function
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericMatrix computeR(NumericVector& corr, NumericMatrix& xMat, NumericMatrix& Rcpp_yMat, double pw) {
  //computeR returns covariance matrix S to use in glasso
  
  int n = xMat.nrow();
  int p = xMat.ncol();
  int numCoef = Rcpp_yMat.ncol();
  NumericMatrix R(n,n); //Correlation matrix
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
  
  return (R);
}