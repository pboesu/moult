\name{sanderlings}
\alias{sanderlings}
\docType{data}
\title{Sanderling Moult Data
}
\description{This data set gives moult indices for 164 Sanderlings trapped on 11 days.
}
\usage{data(sanderlings)}
\format{
  A data frame with 164 observations on the following 2 variables.
  \describe{
    \item{\code{Day}}{a numeric vector of day bird was measured, 1 = 1 July}
    \item{\code{MIndex}}{a numeric vector of moult indices, 0 = bird has not started moult, 1 = bird has completed moult }
  }
}
\details{This data set gives moult indices for 164 Sanderlings trapped on 11 days in the southwestern Cape, South Africa, between October 1978 and April 1979. Day 1 = 1 July). Moult indices are a transformation of moult scores so that moult index increases linearly with time. See Underhill and Zucchini (1988) for details.
}
\source{Underhill and Zucchini (1998)}

\references{Underhill, L. G. and Zucchini, W. (1998). A model for avian primary moult. Ibis 130: 358-372.

}
\examples{
data(sanderlings)

days <- sanderlings$Day                   
MIndex <- sanderlings$MIndex                 # moult indeces

m1 <- moult(MIndex,days,type=1)              # fit model of type 1 to data
m1
m2 <- moult(MIndex,days)                     # model of type 2 (default)
m2
m3 <- moult(MIndex,days,type=3)              # model of type 3
m3

uza <- - m2$coefficients$mean / m2$coefficients$duration    # intercept
uzb <- 1 / m2$coefficients$duration                         # slope


                			# ----- extract how many birds observed on each of the days
nn <- as.numeric(tapply(sanderlings$Day, sanderlings$Day, length))        
nn
					# ----- extract days of observations
day <- unique(sanderlings$Day)                                            
day

p0 <- predict(m2)
p0

p1 <- predict(m2, newdata=data.frame(day))        ## probabilities of moult stages 
p1$M*nn						## Table 6 in Underhill and Zucchini 1988

days2 <- seq(70,310,by=10)
p2 <- predict(m2, newdata = data.frame(days2))
p2$M*100                                  ## Table 7 in Underhill and Zucchini 1988

p3 <- predict(m3,newdata=data.frame(day))      
p3

# ----- Comparison with regression models -----

MInd <- sanderlings$MIndex[sanderlings$MIndex > 0 & sanderlings$MIndex < 1]
MTime <- sanderlings$Day[sanderlings$MIndex > 0 & sanderlings$MIndex < 1]

lm1 <- lm(MTime~MInd)                           
lm1.int <- coef(lm1)[1]
lm1.slope <- coef(lm1)[2]

lm2 <- lm(MInd~MTime)

plot(MTime, MInd, pch = 19, cex=0.7)
	    # regression of Index on Time
abline(lm2, col = "blue", lwd = 2) 
	    # regression of Time on Index: gives better estimates 
            # for mean start day and duration of moult
abline(-lm1.int/lm1.slope, 1/lm1.slope, col = "orange", lwd =2) 

abline(uza, uzb, col = "red", lty = 2, lwd = 2)
}
\keyword{datasets}
