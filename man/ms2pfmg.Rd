\name{ms2pfmg}
\alias{ms2pfmg}
\title{Convert moult scores to proportion feather mass grown}
\description{Convert moult scores obtained for single feathers into 
             overall proportion of feather mass grown.} 

\usage{ms2pfmg(ms, fm)}

\arguments{
  \item{ms}{vector or single value of moult scores (0 to 5)}
  \item{fm}{vector of relative feather mass, corresponding to each 
            feather in ms}}

\value{Returns a single value between 0 and 1: proportion of 
       total feather mass grown}

\details{ms will usually be a vector of 9 or 10 primary feather scores,
         but single feathers can be given, in which case fm
         will be ignored.}
\references{L. Underhill, W. Zucchini, 1988. A model for avian 
            primary moult. Ibis 130: 358-372.}
\author{Bo Bonnevie}
\note{ms2pfmg assumes that for each feather a moult score of 0 corresponds to 0\%
      of feather grown, 1: 0.125 of feather grown, 2: 0.375, 3: 0.625,
      4: 0.875, 5: 1.0}

\examples{

   	      # relative primary feather mass of the 10 primary feathers 
            # (as proportion of total feather mass) for Sanderlings
fm.sand <- c(0.0385,0.0458,0.0544,0.0680,0.0827,0.1019,0.1199,
             0.1417,0.1604,0.1867)   

ms2pfmg(3,0.2)              # single feather
ms2pfmg(3,1)
ms2pfmg("5555500000",fm.sand)

}