\name{ms2pfmg}
\alias{ms2pfmg}
\title{Convert moult scores to proportion of feather mass grown}
\description{Convert moult scores obtained for single feathers into 
             overall proportion of feather mass grown.} 

\usage{ms2pfmg(ms, fm)}

\arguments{
  \item{ms}{vector of moult scores.  Ideally, each moult score is a
    character string of one individual's individual feather moult
    scores, each between 0 and 5, e.g. "555444000" if only nine
    primaries are of interest.  } 
  \item{fm}{vector of relative feather mass, corresponding to each
    feather in ms.}
}

\value{\code{ms2pfmg} returns a single value between 0 and 1: proportion of
  total feather mass grown.}

\details{\code{ms} will usually be a vector of 9 or 10 primary feather scores,
         but single feathers can be given, in which case \code{fm}
         is ignored.  The method used here assumes that a moult score of 1 for 
         any feather corresponds to 1/8th of the feather grown, 2 
         corresponds to 3/8th = 0.375, etc..  The proportion of feather mass 
         grown is then a weighted sum over all feathers, with weights 
         equal to the relative weight (compared to the total weight) of each 
         feather (Underhill and Zucchini 1988).}

\references{Underhill, L. G. and Zucchini, W. (1988) A model for avian primary 
            moult. Ibis \bold{130}, 358--372.}
            
\author{Bo T. Bonnevie}

\seealso{\code{\link{predict.moult}}, \code{\link{moult}}}

\examples{
## relative primary feather mass of the 10 primary feathers 
## (as proportion of total feather mass) for Sanderlings
fm.sand <- c(0.0385, 0.0458, 0.0544, 0.0680, 0.0827, 0.1019,
             0.1199, 0.1417, 0.1604, 0.1867)   

ms2pfmg(3, 0.2)              # single feather
ms2pfmg(3, 1)
ms2pfmg("5555500000", fm.sand)

}