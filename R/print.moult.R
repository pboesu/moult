print.moult <- function(x, ...)      # see print.lm
 { cat("Coefficients",  "\n")

   print(coef(x)) 

   cat("Converged = ", x$converged,"\n")

 }