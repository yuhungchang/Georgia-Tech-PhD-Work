dev_wrapper <- function(beta,X,Y,nug_thres,...){
  if (is.matrix(X) == FALSE){
    x = as.matrix(X)
  }
  d = ncol(X);
  beta = rep(beta,d);
  dev_val = MRGP_deviance(beta,X,Y,nug_thres,...);
  return(dev_val)
}