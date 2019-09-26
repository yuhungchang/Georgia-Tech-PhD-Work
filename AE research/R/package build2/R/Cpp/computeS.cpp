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
  
  return (List::create(Named("mu") = mu, Named("Rinv") = Rinv, Named("R") = R, Named("Rsd") = Rsd));
}