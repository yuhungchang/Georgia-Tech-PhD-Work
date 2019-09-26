#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <numeric>
#include <iostream>
#include <omp.h>

//  function for bigInnerProd
//
// [[Rcpp::plugins(openmp)]]
template <typename T>
NumericMatrix bigInnerProd(XPtr<BigMatrix> pMat, MatrixAccessor<T> mat) {
  //sym - flag for whether matrix is symmetric
  NumericMatrix ret(pMat->ncol(),pMat->ncol());
#pragma omp parallel for
  for (size_t i=0; i < pMat->ncol(); ++i){
    for (size_t j=i; j < pMat->ncol(); ++j){
      //      cout << "(i,j): " << i << ", " << j << endl;
      //      cout << mat[i] << endl;
      //      cout << mat[j] << endl;
      //      cout << std::inner_product(mat[i], mat[i]+pMat->nrow(), mat[j], 0.0) << endl;
      ret(i,j) = std::inner_product(mat[i], mat[i]+pMat->nrow(), mat[j], 0.0);
    }
  }
  return (ret);
}

// Dispatch function for bigInnerProd
//

// [[Rcpp::export]]
NumericVector bigInnerProd(SEXP pBigMat) {
  // First we have to tell Rcpp what class to use for big.matrix objects.
  // This object stores the attributes of the big.matrix object passed to it
  // by R.
  XPtr<BigMatrix> xpMat(pBigMat);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return bigInnerProd(xpMat, MatrixAccessor<char>(*xpMat));
  case 2:
    return bigInnerProd(xpMat, MatrixAccessor<short>(*xpMat));
  case 4:
    return bigInnerProd(xpMat, MatrixAccessor<int>(*xpMat));
  case 8:
    return bigInnerProd(xpMat, MatrixAccessor<double>(*xpMat));
  default:
    // This case should never be encountered unless the implementation of
    // big.matrix changes, but is necessary to implement shut up compiler
    // warnings.
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}