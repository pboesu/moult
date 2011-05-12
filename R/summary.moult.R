summary.moult <- function(object, ...)      
 { cat("Coefficients",  "\n")
   print(coef(object)) 

   cat("Converged = ", object$converged,"\n")

   cat("Standard Errors", "\n")
   print(object$standard.errors)
 }

