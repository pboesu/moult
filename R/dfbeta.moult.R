dfbeta.moult <- function(model, ...) {
  
  c.original <- unlist(coef(model))
  c.original <- matrix(c.original, nrow = 1)
  
  n <- model$n
  np <- length(c.original)
  
  cfs <- matrix(NA, nrow = n, ncol = np)
  cfs.se <- matrix(NA, nrow = n, ncol = np)
  
  for (i in 1:n) {
    datlessi <- data.frame(model$X[-i,])
    
    mm <- moult(model$formula, data = datlessi)
    cfs[i, ] <- unlist(mm$coef)
    cfs.se[i, ] <- unlist(mm$standard.errors)
  }
  
  dfb <- sweep(cfs, 2, c.original) * -1
  dfb.st <- dfb / cfs.se
  
  plot.ts(dfb.st, plot.type = "single", col = 1:np, ylab = "dfbetas",  
          xlab = "observation number", ... = ... )
  abline(h = c(0, 2/sqrt(n), -2/sqrt(n)), col = "grey")
  legend("topleft", col = 1:np, legend = names(model$coefficients), bty = "n", lty = 1)

  return(list(dfbeta = dfb, dfbetas = dfb.st))
}

