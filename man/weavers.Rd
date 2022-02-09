\name{weavers}
\alias{weavers}
\docType{data}
\title{Weaver Moult Data}
\description{Weaver moult data from the Western Cape, South Africa}

\usage{data(weavers)}

\format{
  A data frame with 7543 observations on the following 4 variables.
  \describe{
    \item{\code{RDate}}{a character vector with dates on which
      individuals were caught, format: yyyy-mm-dd.}
    \item{\code{Sex}}{a numeric vector, 0 = unknown, 1 = male, 2 =
      female, 3 = possibly male, 4 = possibly female.}
    \item{\code{Year}}{year in which individual was caught.}
    \item{\code{Moult}}{a character vector with moult scores for 
  individual primary feathers, either nine or ten, starting with
  innermost primary feather.  0: old feather, 5: new feather, 1 to 4,
  feathers at various stages of growth in between.}
  }
}

\source{
  Oschadleus, D. (2005). Patterns of primary moult in weavers. PhD Thesis. University of Cape Town.
  
  SAFRING, URL: \url{https://safring.birdmap.africa/}
}

\examples{
data(weavers)
head(weavers)

if (is.numeric(weavers$Moult)) {
scores <- format(weavers$Moult, scientific = FALSE, trim = TRUE)
} else {
scores <- weavers$Moult
}

mscores <- substr(scores, 1, 9)
feather.mass <- c(10.4, 10.8, 11.5, 12.8, 14.4, 15.6, 16.3, 15.7, 15.7)

## convert moult scores to proportion of feather mass grown

weavers$pfmg <- ms2pfmg(mscores, feather.mass)
weavers$day <- date2days(weavers$RDate, dateformat = "yyyy-mm-dd", startmonth = 8)
ssex <- ifelse(weavers$Sex == 1 | weavers$Sex == 3, "male", 
  ifelse(weavers$Sex == 2 | weavers$Sex == 4, "female", NA))
weavers$ssex <- as.factor(ssex)

## model with duration and mean start date of moult depending on sex

mmf <- moult(pfmg ~ day | ssex | ssex, data = weavers, type = 3)
summary(mmf)

## predict duration and start of moult (then both) for males and females

ssex <- c("male", "female") 
day <- 150

(p1 <- predict.moult(mmf, newdata = data.frame(day, ssex), predict.type = "duration"))
(p2 <- predict.moult(mmf, newdata = data.frame(day, ssex), predict.type = "start"))
(p3 <- predict.moult(mmf, newdata = data.frame(day, ssex), predict.type = "both"))

}

\keyword{datasets}
