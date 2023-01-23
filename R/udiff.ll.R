udiff.ll <- function(theta, Y, X, Z) {
  if (nrow(Y) != nrow(X)) {stop("wrong structure!")}
  if (nrow(Y) != nrow(Z)) {stop("wrong structure!")}
  if (nrow(X) != nrow(Z)) {stop("wrong structure!")}
  #browser()
  if (length(theta) != ncol(Y)*(ncol(Z) + 1) + ncol(Y)*ncol(X) + ncol(Z))
  {stop("the length of theta is not correct!")}
  end1 = ncol(Y)*(ncol(Z)+1)
  end2 = ncol(Y)*ncol(X)
  end3 = ncol(Z)
  #browser()
  W = cbind(matrix(1,nrow(X),1), Z)
  theta.y = theta[1:end1]
  theta.y = matrix(c(theta.y), nrow=ncol(Z)+1, ncol=ncol(Y))
  psi.y = theta[(end1+1):(end1+end2)]
  psi.y = matrix(c(psi.y), nrow=ncol(X), ncol=ncol(Y))
  phi = theta[(end1+end2+1):(end1+end2+end3)]
  phi = matrix(c(phi), nrow=ncol(Z), ncol=1)
  expZ = exp(Z %*% phi) # [nrow * 1]
  sump = W %*% theta.y  + (X %*% psi.y * c(expZ)) # add expZ for each col
  ll = rowSums(Y * sump) - log(1+rowSums(exp(sump)))
  #print(sum(log(-ll)))
  #sum(log(-ll))
  #print(sum(ll))
  -sum(ll)
}


