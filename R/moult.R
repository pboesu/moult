moult <- function(formula, data = NULL, start = NULL, type = 2, method = "BFGS")

{ f.dat <- data
  
  FF <- Formula(formula)

  switch(length(FF)[2], FF <- update(FF, . ~ . -1 | 1 | 1 | 1),
                        FF <- update(FF, . ~ . -1 | . | 1 | 1),
                        FF <- update(FF, . ~ . -1 | . | . | 1),
                        FF <- update(FF, . ~ . -1 | . | . | .))

  call <- match.call()

  if (missing(data)) 
    f.dat <- model.frame(FF)
  
  f.dat <- na.omit(f.dat)
  
  mf <- model.frame(FF, data = f.dat)

  mm <- model.matrix(FF, rhs = 3, data = f.dat)
  md <- model.matrix(FF, rhs = 2, data = f.dat)

  msd <- model.matrix(FF, rhs = 4, data = f.dat)

  if (dim(msd)[2] > 1) {
    FF <- update(FF, . ~ . | . | . | . -1)
    msd <- model.matrix(FF, rhs = 4, data = f.dat)
  } 
  
  Day <- model.matrix(FF, rhs = 1, data = f.dat)
  PFMG <- as.vector(model.part(FF, lhs = 1, data = f.dat))
  M0 <- Day[PFMG == 0]
  MInd <- PFMG[PFMG > 0 & PFMG < 1]
  MTime <- Day[PFMG > 0 & PFMG < 1]
  M1 <- Day[PFMG == 1]
  
  mf0 <- matrix(mm[PFMG == 0,], ncol = ncol(mm))                 # mean start day 
  mfInd <- matrix(mm[PFMG > 0 & PFMG < 1,], ncol = ncol(mm))
  mf1 <- matrix(mm[PFMG == 1,], ncol = ncol(mm))
 
  df0 <- matrix(md[PFMG == 0,], ncol = ncol(md))                     # duration
  dfInd <- matrix(md[PFMG > 0 & PFMG < 1,], ncol = ncol(md))
  df1 <- matrix(md[PFMG == 1,], ncol = ncol(md))
 
  sd0 <- matrix(msd[PFMG == 0,], ncol = ncol(msd))                   # SD(start)
  sdInd <- matrix(msd[PFMG > 0 & PFMG < 1,], ncol = ncol(msd))
  sd1 <- matrix(msd[PFMG == 1,], ncol = ncol(msd))

	## linear regression (only data of birds actively moulting):
  Flin1 <- lm(MTime ~ MInd)                           
  
  p1 <- ncol(dfInd)
  p2 <- ncol(mfInd)
  p3 <- ncol(msd)
 
  no.params <- p1 + p2 + p3  

  colnames(md)[1] <- "intercept.1"
  colnames(mm)[1] <- "intercept"

  paramnames <- c(colnames(md), colnames(mm), colnames(msd))


# ----- 1. Construct likelihood function -----
  
  switch(type,     # --- select likelihood function according to data type (1 to 5)

    { # ----- data of type 1 -----

      LogLikM <- function(p)
      { P <- sum(log(1 - pnorm(M0, mean = mf0 %*% p[(p1 + 1):(p1 + p2)], 
                               sd = sd0 %*% exp(p[(p1+p2+1):no.params]) )))
   
        Q <- sum (log(pnorm(MTime, mean = mfInd %*% p[(p1+1):(p1+p2)], 
                            sd = sdInd %*% exp(p[(p1+p2+1):no.params])) - 
                      pnorm(MTime - (dfInd %*% p[1:p1]), mean = mfInd %*% p[(p1 + 1):(p1 + p2)], 
                            sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params])) ))
   	
        R <- sum(log(pnorm(M1 - (df1 %*% p[1:p1]), mean = mf1 %*% p[(p1 + 1) : (p1 + p2)], 
                           sd = sd1 %*% exp(p[(p1+p2+1):no.params]))))
  
        loglik <- P + Q + R 
        return(-loglik)
      }

    },
   
    { # ----- data of type 2 -----
      LogLikM <- function(p)
      { P.p <- log ( 1 - pnorm(M0, mean=mf0%*%p[(p1+1):(p1+p2)], sd= sd0 %*% exp(p[(p1+p2+1):no.params])))
        P <- sum(P.p[P.p > -Inf])
   
        q <- sum ( log( (dfInd %*% p[1:p1]) * dnorm (MTime - (MInd * (dfInd %*% p[1:p1] )), 
                                                     mean = mfInd %*%p [(p1 + 1):(p1 + p2)],
                                                     sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params]) )))
   	
        R <- sum ( log ( pnorm (M1 - (df1 %*% p[1:p1]), mean = mf1 %*% p[(p1 + 1):(p1 + p2)], 
                                sd = sd1 %*% exp(p[(p1 + p2 + 1):no.params])) ))

        loglik <- P + q + R
        return(-loglik)
      }
    },

    { # ----- data of type 3 -----
        
      LogLikM <- function(p)
      { qq <- sum(log((dfInd %*% p[1:p1]) * dnorm(MTime - (MInd * (dfInd %*% p[1:p1])),
                                                  mean = mfInd%*%p[(p1 + 1):(p1 + p2)], 
                                                  sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params]))))
	
        Q <- sum(log(pnorm(MTime, mean = mfInd %*% p[(p1 + 1):(p1 + p2)],
                           sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params])) - 
                 pnorm(MTime - (dfInd %*% p[1:p1]), mean = mfInd %*% p[(p1 + 1):(p1 + p2)],
                       sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params]))))
   
        loglik <- qq - Q                      # eq.6
        return(-loglik)
      }
    },

    { # ----- data of type 4 -----

      LogLikM <- function(p)
      { qq <- sum(log((dfInd %*% p[1:p1]) * dnorm(MTime - (MInd * (dfInd %*% p[1:p1])),
                                                  mean = mfInd %*% p[(p1 + 1):(p1 + p2)], 
                                                  sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params]))))
	
        P1 <- sum ( log (pnorm(MTime, mean = mfInd %*% p[(p1 + 1):(p1 + p2)],
                               sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params]))))
    
        R <- sum ( log ( pnorm (M1 - (df1 %*% p[1:p1]), 
                                mean = mf1 %*% p[(p1 + 1):(p1 + p2)], 
                                sd = sd1 %*% exp(p[(p1 + p2 + 1):no.params])) ))

        P2 <- sum ( log (pnorm(M1, mean = mf1 %*% p[(p1 + 1):(p1 + p2)], 
                               sd = sd1 %*% exp(p[(p1 + p2 + 1):no.params]))))

        loglik <- qq - P1 + R - P2
        return(-loglik)
      }
    },

    { # ----- data of type 5 -----

      LogLikM <- function(p)
      { P <- sum ( log ( 1 - pnorm(M0, mean = mf0 %*% p[(p1 + 1):(p1 + p2)], 
                                   sd = sd0 %*% exp(p[(p1 + p2 + 1):no.params]))))
    
        R1 <- sum ( log ( 1 - pnorm (M0 - (df0 %*% p[1:p1]), 
                                     mean = mf0 %*% p[(p1 + 1):(p1 + p2)], 
                                     sd = sd0 %*% exp(p[(p1 + p2 + 1):no.params])) ))

        qq <- sum( log ((dfInd %*% p[1:p1]) * dnorm(MTime - (MInd * (dfInd %*% p[1:p1])),
                                                    mean = mfInd %*% p[(p1 + 1):(p1 + p2)],
                                                    sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params]))))
    
        R2 <- sum ( log (1 - pnorm (MTime - (dfInd %*% p[1:p1]), 
                                    mean = mfInd %*% p[(p1 + 1):(p1 + p2)], 
                                    sd = sdInd %*% exp(p[(p1 + p2 + 1):no.params])) ))
   

        loglik <- P - R1 + qq - R2
        return(-loglik)

      }
    },
    
    print("moult: not a valid data type") ) # --- end of switch

   # ----- 2. Run optimisation -----

  if (missing(start))
    InitP <- c(coef(Flin1)[2], rep(0, times = p1 - 1), coef(Flin1)[1], rep(0, times = p2 - 1), 
               rep(log(sd(MTime, na.rm = TRUE)), times = p3))
  else 
    InitP <- start

  fit <- optim(par = InitP, LogLikM, hessian = TRUE, control = list(maxit = 100000), method = method)

   # ----- 3. Return output -----

  vcnames <- c(paste("d.", colnames(md), sep = ""), 
               paste("s.", colnames(mm), sep = ""), 
               paste("sd.", colnames(msd), sep = ""))
  names(fit$par) <- paramnames
  colnames(fit$hessian) <- vcnames
  rownames(fit$hessian) <- vcnames

  out <- list("estimates" = fit$par, "likelihood" = -fit$value, "convergence" = fit$convergence, 
              "message" = fit$message, "hessian" = fit$hessian)
  
  coefd <- fit$par[1:p1]
  coefm <- fit$par[(p1 + 1):(p1 + p2)]
  coefsd <- fit$par[(p1 + p2 + 1):no.params]
  km <- p2
  kd <- p1
  nobs <- length(PFMG)

  dur.est <- md %*% coefd
  mean.est <- mm %*% coefm

  moult.est <- ifelse (Day < mean.est, 0, ifelse (Day > mean.est + dur.est, 1,
                      - mean.est / dur.est + 1 / dur.est * Day))

  residuals <- PFMG - moult.est

  vcov <- solve(as.matrix(fit$hessian))
  dvc <- diag(vcov)
 
  ses <- sqrt(dvc)
  ses[(p1 + p2 + 1):no.params] <- exp(coefsd) * sqrt(ses[(p1 + p2 + 1):no.params])
  
  names(ses) <- paramnames
   
  formula.duration <- formula(FF, rhs = 2, lhs = -1)
  formula.mean <- formula(FF, rhs = 3, lhs = -1)
  formula.sd <- formula(FF, rhs = 4, lhs = -1)
  
  out <- list(coefficients = list(duration = coefd, mean = coefm, sd = exp(coefsd)),
              loglik = - fit$value,
              vcov = vcov,
              standard.errors = ses,
              type = type,
              residuals = residuals,
              fitted.values = moult.est,
              n = nobs,
              df.residual = nobs - no.params,
              terms = list(full = FF, duration = formula.duration, mean = formula.mean, sd = formula.sd), 
              call = call,
              formula = FF, 
              optim = fit,             
              converged = fit$convergence < 1,
              convergence.value = fit$convergence)
  
  class(out) <- "moult"
  return(out)
}



