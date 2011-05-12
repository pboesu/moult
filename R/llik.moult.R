llik.moult <-
function(p, moult.scores, days, type = 2)
{ Day <- days
  PFMG <- moult.scores
 
  M.all <- cbind(Day,PFMG) 
  M.all <- na.omit(M.all)
 
  Day <- M.all[,1]
  PFMG <- M.all[,2]

  M0 <- Day[PFMG == 0]
  MInd <- PFMG[PFMG > 0 & PFMG < 1]
  MTime <- Day[PFMG > 0 & PFMG < 1]
  M1 <- Day[PFMG == 1]

  switch(type, 
      
      { P <- sum(log(1 - pnorm(M0,mean=p[2],sd=exp(p[3]))))			# ---- type 1
        Q <- sum (log(pnorm(MTime, mean = p[2], sd = exp(p[3])) - 
                      pnorm(MTime - p[1], mean = p[2], sd = exp(p[3])) ))
        R <- sum(log(pnorm(M1 - p[1],mean=p[2],sd=exp(p[3]))))
        loglik <- P + Q + R 
      },

      { P <- sum ( log ( 1 - pnorm(M0, mean=p[2], sd=exp(p[3]))))			# ---- type 2
        q <- sum ( log( p[1] * dnorm (MTime - (MInd * p[1]), mean = p[2], sd = exp(p[3]) )))
        R <- sum ( log ( pnorm (M1 - p[1], mean = p[2], sd = exp(p[3])) ))
        loglik <- P + q + R 
      },

      { qq <- sum(log(p[1] * dnorm(MTime - (MInd * p[1]),mean=p[2],sd=exp(p[3]))))			# ---- type 3
        Q <- sum(log(pnorm(MTime,mean=p[2],sd=exp(p[3])) - 
                     pnorm(MTime - p[1],mean=p[2],sd=exp(p[3]))))
        loglik <- qq - Q                      # eq.6
      },

      { qq <- sum(log(p[1] * dnorm(MTime - (MInd * p[1]),mean=p[2],sd=exp(p[3]))))			# ---- type 4
        P1 <- sum ( log (pnorm(MTime, mean=p[2], sd=exp(p[3]))))
        R <- sum ( log ( pnorm (M1 - p[1], mean = p[2], sd = exp(p[3])) ))
        P2 <- sum ( log (pnorm(M1, mean=p[2], sd=exp(p[3]))))
        loglik <- qq - P1 + R - P2
      },

      { P <- sum ( log (1 - pnorm(M0, mean=p[2], sd=exp(p[3]))))			# ---- type 5
        R1 <- sum ( log ( 1 - pnorm (M0 - p[1], mean = p[2], sd = exp(p[3])) ))
        qq <- sum( log (p[1] * dnorm(MTime - (MInd * p[1]),mean=p[2],sd=exp(p[3]))))
        R2 <- sum ( log (1 - pnorm (MTime - p[1], mean = p[2], sd = exp(p[3])) ))
        loglik <- P - R1 + qq - R2
      }, 

      print("llik.moult: invalid type"))
   
  return(loglik)

}

