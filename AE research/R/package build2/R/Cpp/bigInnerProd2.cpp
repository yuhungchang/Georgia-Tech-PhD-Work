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
NumericMatrix bigInnerProd2(XPtr<BigMatrix> pMat, XPtr<BigMatrix> pMat2, MatrixAccessor<T> mat, MatrixAccessor<T> mat2) {
  NumericMatrix ret(pMat->ncol(),pMat2->ncol());
#pragma omp parallel for
  for (size_t i=0; i < pMat->ncol(); ++i){
    for (size_t j=0; j < pMat2->ncol(); ++j){
      //      cout << "(i,j): " << i << ", " << j << endl;
      //      cout << mat[i] << endl;
      //      cout << mat[j] << endl;
      //      cout << std::inner_product(mat[i], mat[i]+pMat->nrow(), mat[j], 0.0) << endl;
      ret(i,j) = std::inner_product(mat[i], mat[i]+pMat->nrow(), mat2[j], 0.0);
    }
  }
  return (ret);
}

// Dispatch function for bigInnerProd
//

// [[Rcpp::export]]
NumericMatrix bigInnerProd2(SEXP pBigMat, SEXP pBigMat2) {
  // First we have to tell Rcpp what class to use for big.matrix objects.
  // This object stores the attributes of the big.matrix object passed to it
  // by R.
  XPtr<BigMatrix> xpMat(pBigMat);
  XPtr<BigMatrix> xpMat2(pBigMat2);
  
  switch(xpMat->matrix_type()) {
  case 1:
    return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<char>(*xpMat), MatrixAccessor<char>(*xpMat2));
  case 2:
    return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<short>(*xpMat), MatrixAccessor<short>(*xpMat2));
  case 4:
    return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<int>(*xpMat), MatrixAccessor<int>(*xpMat2));
  case 8:
    return bigInnerProd2(xpMat, xpMat2, MatrixAccessor<double>(*xpMat), MatrixAccessor<double>(*xpMat2));
  default:
    // This case should never be encountered unless the implementation of
    // big.matrix changes, but is necessary to implement shut up compiler
    // warnings.
    throw Rcpp::exception("unknown type detected for big.matrix object!");
  }
}