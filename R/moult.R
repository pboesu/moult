moult <- function(formula, data = NULL, start = NULL, type = 2, method = "BFGS",
                  fixed = NULL, fixed.val = NULL, prec = 0.02) 
    
{ f.dat <- data.frame(data)
  
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
  
  Day <- model.frame(FF, rhs = 1, data = f.dat) 
  PFMG <- model.part(FF, lhs = 1, data = f.dat)
  
  if (max(PFMG) > 1) warning("Some moult index values > 1")
  
  M0 <- Day[PFMG == 0, 2]
  MInd <- PFMG[PFMG > 0 & PFMG < 1]
  MTime <- Day[PFMG > 0 & PFMG < 1, 2]
  M1 <- Day[PFMG == 1, 2]
  
  mf0 <- matrix(mm[PFMG == 0,], ncol = ncol(mm))                 # mean start day 
  mfInd <- matrix(mm[PFMG > 0 & PFMG < 1,], ncol = ncol(mm))
  mf1 <- matrix(mm[PFMG == 1,], ncol = ncol(mm))
 
  df0 <- matrix(md[PFMG == 0,], ncol = ncol(md))                     # duration
  dfInd <- matrix(md[PFMG > 0 & PFMG < 1,], ncol = ncol(md))
  df1 <- matrix(md[PFMG == 1,], ncol = ncol(md))
 
  sd0 <- matrix(msd[PFMG == 0,], ncol = ncol(msd))                   # SD(start)
  sdInd <- matrix(msd[PFMG > 0 & PFMG < 1,], ncol = ncol(msd))
  sd1 <- matrix(msd[PFMG == 1,], ncol = ncol(msd))

  MIndmin <- pmax(0, MInd - prec)
  MIndmax <- pmin(MInd + prec, 1)
  
  n1 <- dim(mf0)[1]
  n2 <- dim(mfInd)[1]
  n3 <- dim(mf1)[1]
  
	## linear regression (only data of birds actively moulting):
  Flin1 <- lm(MTime ~ MInd)                           
  
  p1 <- ncol(dfInd)
  p2 <- ncol(mfInd)
  p3 <- ncol(msd)
 
  no.params <- p1 + p2 + p3
   
  colnames(md)[1] <- "duration"
  colnames(mm)[1] <- "mean-start-day"
  colnames(msd)[1] <- "SD-start-day"
  
  paramnames <- c(colnames(md), colnames(mm), colnames(msd))

  lik.P <- function(pdur, pstart, psd) {
    P.p <- log ( 1 - pnorm(M0, mean = mf0 %*% pstart,
                           sd = sd0 %*% psd))
    return(P.p[P.p > -Inf])
  }
  
  lik.Q <- function(pdur, pstart, psd) {
    Q <- log(pnorm(MTime, mean = mfInd %*% pstart, 
                   sd = sdInd %*% psd) - 
               pnorm(MTime - (dfInd %*% pdur), mean = mfInd %*% pstart, 
                     sd = sdInd %*% psd))
  return(Q)
}

lik.R <- function(pdur, pstart, psd) {
  
  R <- log(pnorm(M1 - (df1 %*% pdur), mean = mf1 %*% pstart, 
                 sd = sd1 %*% psd))
  return(R)
}

lik.q <- function(pdur, pstart, psd) {
  q <- log(pnorm(MTime - MIndmin * (dfInd %*% pdur),
                 mean = mfInd %*% pstart, 
                 sd = sdInd %*% psd) -
             pnorm(MTime - MIndmax * (dfInd %*% pdur),
                   mean = mfInd %*% pstart, 
                   sd = sdInd %*% psd))
  return(q)
}

lik.P1 <- function(pdur, pstart, psd) {
  P1 <- log (pnorm(MTime, mean = mfInd %*% pstart,
                   sd = sdInd %*% psd))
  return(P1)
}

lik.P2 <- function(pdur, pstart, psd) {
  P2 <- log (pnorm(M1, mean = mf1 %*% pstart, 
                   sd = sd1 %*% psd))
  return(P2)
}

lik.R1 <- function(pdur, pstart, psd) {
  R1 <- log ( 1 - pnorm (M0 - (df0 %*% pdur), 
                         mean = mf0 %*% pstart, 
                         sd = sd0 %*% psd))
  return(R1)
}

lik.R2 <- function(pdur, pstart, psd) {
  R2 <- log (1 - pnorm (MTime - (dfInd %*% pdur), 
                        mean = mfInd %*% pstart, 
                        sd = sdInd %*% psd))
  return(R2)
}


  # ----- 1. Construct likelihood function -----
  
  LogLikM <- function(p, .fixed = fixed, .fixed.val = fixed.val) {
    
    if (is.null(.fixed)) {
      params <- p
    } else {
      if(length(p) + length(.fixed.val) != length(.fixed))
        stop("number of parameters not right, fixed")
      params <- numeric(length(.fixed))
      params[.fixed] <- .fixed.val
      params[!.fixed] <- p
      params
    }
    
    pdur <- p[1:p1]
    pstart <- p[(p1 + 1):(p1 + p2)]
    psd <- exp(p[(p1 + p2 + 1):no.params])
 
    switch(type,
             { P <- sum(lik.P(pdur, pstart, psd))
               Q <- sum(lik.Q(pdur, pstart, psd))
               R <- sum(lik.R(pdur, pstart, psd))
               loglik <- P + Q + R  # type 1
               return(-loglik) 
             },
             
             { P <- sum(lik.P(pdur, pstart, psd))
               q <- sum(lik.q(pdur, pstart, psd))
               R <- sum(lik.R(pdur, pstart, psd))
               loglik <- P + q + R   # type 2
               return(-loglik) 
             },
             
             { q <- sum(lik.q(pdur, pstart, psd))
               Q <- sum(lik.Q(pdur, pstart, psd))
          
               if (!is.na(Q) & !is.na(q)) {
                 if (Q < q) {
                   stop("Q < qq, type 3")
                   cat(list(Q = Q, q = q, params = params))
                 }
               }
               loglik <- q - Q                      # type 3
               return(-loglik)
             },
             
             { P1 <- sum(lik.P1(pdur, pstart, psd))
               P2 <- sum(lik.P2(pdur, pstart, psd))
               q <- sum(lik.q(pdur, pstart, psd))
               R <- sum(lik.R(pdur, pstart, psd))
               
               loglik <- q - P1 + R - P2           # type 4
               return(-loglik)
             },
             
             { R1 <- sum(lik.R1(pdur, pstart, psd))
               R2 <- sum(lik.R2(pdur, pstart, psd))
               q <- sum(lik.q(pdur, pstart, psd))
               P <- sum(lik.P(pdur, pstart, psd))
             
               loglik <- P - R1 + q - R2                     # type 5
               return(-loglik)
             },
             
             print("moult: not a valid data type") ) # --- end of switch)
}
      nobs <- switch(type,
                   sum(length(M0), length(M1), length(MInd)), 
                   sum(length(M0), length(M1), length(MInd)),
                   length(MInd),
                   sum(length(M1), length(MInd)),
                   sum(length(M0), length(MInd)))

    ## ===============================    
    ## ----- 2. Run optimisation -----
    ## ===============================
    
    if (missing(start))
        InitP <- c(max(coef(Flin1)[2], 1), rep(0, times = p1 - 1),
                   coef(Flin1)[1], rep(0, times = p2 - 1), 
                   rep(log(sd(MTime, na.rm = TRUE)), times = p3))
    else 
        InitP <- start
    
    if(!is.null(fixed)) InitP <- InitP[!fixed]
          
    fit <- optim(par = InitP, LogLikM, hessian = TRUE, control = list(maxit = 100000),
                 method = method)

    ## ============================
    ## ----- 3. Return output -----
    ## ============================
    
    vcnames <- c(paste("d.", colnames(md), sep = ""), 
                 paste("s.", colnames(mm), sep = ""), 
                 paste("sd.", colnames(msd), sep = ""))
    
    p.est <- numeric(length(no.params))
    if (!is.null(fixed)) {
        p.est[!fixed] <- fit$par
        p.est[fixed] <- fixed.val
    } else {
        p.est <- fit$par
    }
    
    names(p.est) <- paramnames
    
    coefd <- p.est[1:p1]
    coefm <- p.est[(p1 + 1):(p1 + p2)]
    coefsd <- p.est[(p1 + p2 + 1):no.params]
    km <- p2
    kd <- p1

    dur.est <- md %*% coefd
    mean.est <- mm %*% coefm

    moult.est <- ifelse (Day[, 2] < mean.est, 0,
                         ifelse (Day[, 2] > mean.est + dur.est, 1,
                                 - mean.est / dur.est + 1 / dur.est * Day[, 2]))
    residuals <- PFMG - moult.est

    
    H.inv <- solve(fit$hessian)

    vcov <- matrix(0, nrow = no.params, ncol = no.params)

    if (is.null(fixed)) {
        vcov <- H.inv
    } else {
        vcov[!fixed, !fixed] <- H.inv
    }
    
    dvc <- diag(vcov)
    
    colnames(vcov) <- vcnames
    rownames(vcov) <- vcnames
    
    ses <- sqrt(dvc)
    ses[(p1 + p2 + 1):no.params] <- exp(coefsd) * (ses[(p1 + p2 + 1):no.params])
  
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
                terms = list(full = FF, duration = formula.duration, mean = formula.mean,
                             sd = formula.sd), 
                call = call,
                X = mf,
                y = PFMG,
                Day = Day[, 2],
                formula = FF, 
                optim = fit,             
                converged = fit$convergence < 1,
                convergence.value = fit$convergence)

    class(out) <- "moult"
    return(out)
}



