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
  
  return (List::create(Named("Rnew") = Rnew, Named("Rinv") = Rinv, Named("pred") = pred, Named("var") = var));
  
}