udiff <- function (formula, data, subset, weights, na.action,
                contrasts = NULL, offset, ...)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mf[] <- lapply(mf, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
  })
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  w <- as.vector(model.weights(mf))
  if(!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if(!is.null(offset)) {
    if(length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  
  if (is.empty.model(mt)) {
    stop("Model cannot be empty!")
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
 #   y <- as.factor(y)
    x <- data.frame(x[,2:ncol(x)])
    if (ncol(x) != 2) {stop("The independent variables should be two factors")}
#    x[,1] <- as.factor(x[,1])
#    x[,2] <- as.factor(x[,2])
    
    z <- if(is.null(w)) udiff.fit(x, y, offset = offset, ...)
    else udiff.wfit(x, y, w, offset = offset, ...)
  }
  class(z) <- c(if(is.matrix(y)) "udiff")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z
}

udiff.fit <- function (x, y, offset = NULL, tol = 1e-07,
                      ...)
{
  if (is.null(n <- nrow(x))) stop("'x' must be a matrix")
  if(n == 0L) stop("0 (non-NA) cases")
#  if (!is.factor(y)) { stop("The dependent variable should be a factor.")}
  if (ncol(x) != 2) { stop("The independent variables should be two factors.")}
#  if (!is.factor(x[,1])) {stop("The independent variables should be two factors.")}
#  if (!is.factor(x[,2])) {stop("The independent variables should be two factors.")}
  if (length(unique(y)) == 1) { stop("The dependent variables have only one level.")}
  if (length(unique(x[,1])) == 1) { stop("The independent variables have only one level.")}
  if (length(unique(x[,2])) == 1) { stop("The independent variables have only one level.")}
  vY = vdummy(y)[,2:length(unique(y))]
  vX = vdummy(x[,1])[,2:length(unique(x[,1]))]
  vZ = vdummy(x[,2])[,2:length(unique(x[,2]))]
  nelem = ncol(vY)*(ncol(vZ) + 1) + ncol(vY)*ncol(vX) + ncol(vZ)
  theta = runif(rep(1, nelem),-1, 1)
  chkDots(...)
  ## Generate the start value for theta
  ll <- partial(udiff.ll,
                Y = vY, X = vX, Z = vZ)
  gr <- gen_gr(ll)
  z <- optim(theta, ll, gr, method = "BFGS", hessian = T)
  # z <- optim(theta, ll)
  coef <- z$par
  df <-  nrow(vY) - length(coef) - 1
  std <- sqrt(diag(solve(z$hessian)))
  tval <- coef / std
  pval <- 2 * pt(-abs(tval), df)
  coefmat <- t(rbind(coef, std, tval, pval))
  names.x <- unique(x[,1])[2:length(unique(x[,1]))]
  names.z <- unique(x[,2])[2:length(unique(x[,2]))]
  names.y <- unique(y)[2:length(unique(y))]
  names_coef <- c() 
  for (i in names.y) {
    names_coef <- c(names_coef, paste("theta", i, "Cons", sep="."))
    for (j in names.z) {
      ntheta <- paste("theta", i, j, sep=".")
      names_coef <- c(names_coef, ntheta)
    }
  }
  
  for (i in names.y) {
    for (j in names.x) {
      ntheta <- paste("Psi", i, j, sep=".")
      names_coef <- c(names_coef, ntheta)
    }
  }
  
  for (j in names.z) {
    ntheta <- paste("Phi", j, sep=".")
    names_coef <- c(names_coef, ntheta)
  }
  browser()
  rownames(coefmat) <- names_coef
  obs <- n
  AIC <- 2 * df - 2 * z$value
  BIC <- df * log(obs) - 2 * z$value
  # pivot <- z$pivot
  ## careful here: the rank might be 0
  # r1 <- seq_len(z$rank)
  # dn <- colnames(x); if(is.null(dn)) dn <- paste0("x", 1L:p)
  # nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
  # r2 <- if(z$rank < p) (z$rank+1L):p else integer()\
 
  # r1 <- y - z$residuals ; if(!is.null(offset)) r1 <- r1 + offset
  ## avoid unnecessary copy
  list(coefficients=coefmat, obs=obs, AIC=AIC, BIC=BIC)
}

