
predict.moult <- function(object, newdata = NULL, what = "distrib", intervals = 0.1, ...)
{ 		# --- predicted proportions not in moult, in moult, and completed

   if (missing(newdata) | what == "response")       # will just predict moult score (fitted values)
     
    { return(object$fitted.values) }

   else         # if newdata is given predict probability that in certain moult stage
    
    { p <- c(object$coefficients,recursive = TRUE)
      p <- as.vector(p)

      type <- object$type

      p1 <- length(object$coefficients$duration)
      p2 <- length(object$coefficients$mean)
      p3 <- length(object$coefficients$sd)
      no.params <- p1 + p2 + p3

      p.dur <- p[1:p1]
      p.mean <- p[(p1+1):(p1+p2)]
      p.sd <- p[(p1+p2+1):no.params]

#      ff <- object$formula
#      ff <- update(ff, ~. - moult.scores, evaluate = FALSE)
#      print(ff)

#      mf <- model.frame(ff,data=newdata)

#      print(head(mf))

      day <- newdata[,1]
      mm <- model.matrix(object$mean.formula, newdata)
      md <- model.matrix(object$duration.formula, newdata)

      if(!is.null(object$sd.formula))
       { msd <- model.matrix(object$sd.formula,newdata) }
      else
       { msd <- model.matrix(~ 1, newdata) }

      # -----      

      no.moult <- 1 - pnorm(day, mean = mm %*% p.mean, sd = msd %*% p.sd)
      
      in.moult <- pnorm(day, mean = mm %*% p.mean, sd = msd %*% p.sd) - 
                   pnorm(day - p[1], mean = mm %*% p.mean, sd = msd %*% p.sd ) 
      compl.moult <- pnorm(day - p[1], mean = mm %*% p.mean, sd = msd %*% p.sd )

      
      NIC <- cbind(no.moult, in.moult, compl.moult)
      colnames(NIC) <- c("not started", "in moult", "completed")

      ints <- seq(0, 1, by=intervals)
           
      IM <- matrix(0,ncol=length(ints)-1,nrow=length(day))
      rownames(IM) <- day
      colnames(IM) <- paste(ints[1:length(ints)-1],"-",ints[2:length(ints)])
      for (i in 2:length(ints))
      { IM[,i-1] <- pnorm(day - md %*% p.dur * ints[i-1], mean = mm %*% p.mean, sd = msd %*% p.sd) - 
                    pnorm(day - md %*% p.dur * ints[i], mean = mm %*% p.mean, sd = msd %*% p.sd) 
      }

      switch(type, 
       { M <- round(cbind(NIC[,1],IM, NIC[,3]),digits=3) },  # --- type 1
       { M <- round(cbind(NIC[,1],IM, NIC[,3]),digits=3) }, # --- type 2
       { M <- round(cbind(NIC[,1],IM / in.moult , NIC[,3]), digits=3) },  # --- type 3
       { M <- round(cbind(NIC[,1],IM / ( 1 - no.moult), compl.moult/(1 - no.moult)), digits=3) },  # --- type 4
       { M <- round(cbind(no.moult/(1 - compl.moult), IM / ( 1 - compl.moult), compl.moult), digits=3) },     # --- type 5                       
       print("predict.moult: not a valid type") )                   # ---- end of switch

      colnames(M) <- c("0",colnames(IM),"1")
      out <- list("M" = M)

   }      

   return(out)

}



