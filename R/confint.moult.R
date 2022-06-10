confint.moult <- function(object, parm, level = 0.95, ..., B = 1000, add.plot = TRUE) {
  
  cf <- object$coefficients
  ses <- object$standard.errors
  pnames <- names(object$coefficients)
  ff <- object$formula
  
  B <- B
  n <- object$n
  np <- length(cf)
  
  dur.est <- numeric(B)
  start.est <- numeric(B)
  sd.est <- numeric(B)
  half.est <- numeric(B)
  end.est <- numeric(B)
  
  ests <- matrix(NA, nrow = B, ncol = np)
  
  for (i in 1:B) {
    samp.ind <- sample(1:n, replace = TRUE)
    boot.samp <- object$X[samp.ind,]
    boot.out <- moult(ff, data = boot.samp)
    est <- boot.out$coefficients
    ests[i, ] <- unlist(est)
    half.est[i] <- est$mean + 0.5 * est$duration
    end.est[i] <- est$mean + est$duration
  }
  
  all <- cbind(half.est, end.est)
  all.df <- cbind(ests, all)
  colnames(all.df) <- c(pnames, "half.est", "end.est")
 
  if (add.plot == TRUE) {
    pairs(all.df, ... = ...)
  }
  
  ## Covariance Matrix 
  varcov <- var(all.df, na.rm = TRUE)
  
  ### Standard Errors
  ses <- sqrt(diag(varcov))
  
  ## Bootstrap Percentile Intervals
  alpha <- 1 - level
  ci.boot <- apply(all.df, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2))
  
  return(list(bootstrap.distribution = all.df, bootstrap.percentile.ci = ci.boot, 
              bootstrap.vcov = varcov, 
              bootstrap.SE = ses))
}



