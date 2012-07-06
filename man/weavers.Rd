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
  
  SAFRING, URL: \url{http://safring.adu.org.za/}
}

\examples{
data(weavers)
head(weavers)
}

\keyword{datasets}
