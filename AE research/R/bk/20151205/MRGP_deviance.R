MRGP_deviance <- function (beta, X, Y, nug_thres = 20, corr = list(type = "exponential", 
                                                  power = 1.95)) # Y is an n by r matrix
{
  if (is.matrix(X) == FALSE) {
    X = as.matrix(X)
  }
  n = nrow(X)
  d = ncol(X)
  r = ncol(Y)
  if (n != nrow(Y)) {
    stop("The dimensions of X and Y do not match. \n")
  }
  if (d != length(beta)) {
    stop("The dimensions of beta and X do not match \n")
  }
  if (n != nrow(Y)) {
    stop("The dimensions of X and Y do not match \n")
  }
  if (nug_thres < 10) {
    warning("nug_thres outside of the normal range of [10, 25]")
  }
  if (nug_thres > 25) {
    warning("nug_thres outside of the normal range of [10, 25]")
  }
  One = rep(1, n)
  dim(beta) = c(1, d)
  H = corr_matrix(X, beta, corr)
  
  temp = eigen(H, symmetric = TRUE, only.values = TRUE)
  eig_val = temp$values
  condnum = kappa(H, triangular = TRUE, exact = TRUE)
  max_eigval = eig_val[1]
  delta = max(c(0, abs(max_eigval) * (condnum - exp(nug_thres))/(condnum * 
                                                                   (exp(nug_thres) - 1))))
  LO = diag(n)
  Sig = H + delta * LO
  L = chol(Sig)
  Sig_invOne = solve(L, solve(t(L), One))
  Sig_invY = solve(L, solve(t(L), Y))
  Sig_invLp = solve(L, solve(t(L), LO))
  
  mu_hat = solve(t(One) %*% Sig_invOne, t(One) %*% Sig_invY)
  mu_hat_matrix <- matrix(rep(mu_hat, each = n), ncol = r)
  Sig_invb = solve(L, solve(t(L), (Y - mu_hat_matrix)))
  temp2 = eigen(L, only.values = TRUE)
  eig_valL = temp2$values
  T_hat = t(Y - mu_hat_matrix) %*% Sig_invb / n
  part1 = 2 * r * sum(log(abs(eig_valL)))
  part2 = n * log(det(T_hat))
  devval = part1 + part2
  #devval = devval[1, 1]
  if (is.finite(devval) == FALSE) {
  #  stop("Infinite values of the Deviance Function, \n\t\tunable to find optimum parameters \n")
  }
  return(devval)
}