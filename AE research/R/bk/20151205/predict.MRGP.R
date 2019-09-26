predict.MRGP <- function (object, xnew = object$X, M = 1, ...) 
{
  if (is.matrix(xnew) == FALSE) {
    xnew = as.matrix(xnew)
  }
  if (M <= 0) {
    M = 1
    warning("M was assigned non-positive, changed to 1. \n")
  }
  X = object$X
  Y.mx = object$Y
  Y_hat.mx = matrix(0, nrow = nrow(xnew), ncol = ncol(Y.mx))
  for(dd in 1:ncol(Y.mx)){
    Y = Y.mx[,dd]
    n = nrow(X)
    d = ncol(X)
    corr = object$correlation_param
    power = corr$power
    nu = corr$nu
    if (d != ncol(xnew)) {
      stop("The training and prediction data sets are of \n\tdifferent dimensions. \n")
    }
    beta = object$beta
    sig2 = object$sig2
    delta = object$delta
    dim(beta) = c(d, 1)
    H = corr_matrix(X, beta, corr)
    One = rep(1, n)
    LO = diag(n)
    Sig = H + delta * LO
    L = chol(Sig)
    if (delta == 0) {
      Sig_invOne = solve(L, solve(t(L), One))
      Sig_invY = solve(L, solve(t(L), Y))
      Sig_invLp = solve(L, solve(t(L), LO))
    }
    else {
      s_Onei = One
      s_Yi = Y
      s_Li = LO
      Sig_invOne = matrix(0, ncol = 1, nrow = n)
      Sig_invY = matrix(0, ncol = 1, nrow = n)
      Sig_invLp = matrix(0, ncol = n, nrow = n)
      for (it in 1:M) {
        s_Onei = solve(L, solve(t(L), delta * s_Onei))
        Sig_invOne = Sig_invOne + s_Onei/delta
        s_Yi = solve(L, solve(t(L), delta * s_Yi))
        Sig_invY = Sig_invY + s_Yi/delta
        s_Li = solve(L, solve(t(L), delta * s_Li))
        Sig_invLp = Sig_invLp + s_Li/delta
      }
    }
    nnew = nrow(xnew)
    Y_hat = rep(0, nnew)
    MSE = rep(0, nnew)
    if (corr$type == "exponential") {
      for (kk in 1:nnew) {
        xn = matrix(xnew[kk, ], nrow = 1)
        r = exp(-(abs(X - as.matrix(rep(1, n)) %*% (xn))^power) %*% 
                  (10^beta))
        yhat = (((1 - t(r) %*% Sig_invOne)/(t(One) %*% Sig_invOne)) %*% 
                  t(One) + t(r)) %*% Sig_invY
        Y_hat[kk] = yhat
      }
    }
    if (corr$type == "matern") {
      for (kk in 1:nnew) {
        xn = matrix(xnew[kk, ], nrow = 1)
        temp = 10^beta
        temp = matrix(temp, ncol = d, nrow = (length(X)/d), 
                      byrow = TRUE)
        temp = 2 * sqrt(nu) * abs(X - as.matrix(rep(1, n)) %*% 
                                    (xn)) * (temp)
        ID = which(temp == 0)
        rd = (1/(gamma(nu) * 2^(nu - 1))) * (temp^nu) * besselK(temp, 
                                                                nu)
        rd[ID] = 1
        r = matrix(apply(rd, 1, prod), ncol = 1)
        yhat = (((1 - t(r) %*% Sig_invOne)/(t(One) %*% Sig_invOne)) %*% 
                  t(One) + t(r)) %*% Sig_invY
        Y_hat[kk] = yhat
      }
    }
    Y_hat.mx[,dd] = Y_hat
  }
  
  prediction = NULL
  names = c()
  for (i in 1:d) {
    names[i] = paste("xnew.", i, sep = "")
  }
  names[(d + 1) : (d + ncol(Y_hat.mx))] = paste0("Y_hat.", 1:ncol(Y_hat.mx))
  full_pred = cbind(xnew, Y_hat.mx)
  colnames(full_pred) = names
  prediction$Y_hat = Y_hat.mx
  #prediction$MSE = MSE
  prediction$complete_data = full_pred
  return(prediction)
}