\name{moult}
\alias{moult}
\alias{print.moult}
\alias{summary.moult}
\alias{logLik.moult}
\title{Models to analyse data of moult in birds}
\description{Estimate duration and mean start date of moult from moult score data by maximum likelihood. Covariates to model duration and start of moult are possible. }

\usage{moult(moult.scores, days, formula = NULL, sd = NULL, start = NULL, type = 2)}
    
\arguments{
  \item{moult.scores}{vector of moult scores}
  \item{days}{vector of days corresponding to moult.scores}
  \item{formula}{symbolic description of the model for mean start date | duration}
  \item{sd}{if provided must be one grouping factor / vector allowing
  standard deviation for start date to be estimated separately for each
  group}
  \item{start}{starting values for parameters}
  \item{type}{integer (one of 1,2,3,4,5) referring to type of moult data
    and consequently model to be fitted (see details).}
  
}

\details{formula is specified in 2 parts, the first corresponds to a
  formula for mean start date, the second part to a formula for duration
  of moult, see examples. If only one part is supplied, this will be
  taken to be the model for the mean start date. If no formula is
  supplied, a model without covariates will be fitted (i.e. only
  duration, mean start date and standard deviation of start date will be
estimated). 

  
`type' refers to the type of moult data available (see Underhill and
Zucchini (1998) and Underhill, Zucchini and Summers (1990)). 

type 1: Sample is representative of entire population (not yet moulting, in moult, and
birds which have completed moult). For type 1 data, any value between 0 and 1 
(> 0 and < 1) can be used as the moult index for birds in active moult. The value 
used does not matter, only the fact that they are in moult. 

type 2 (default): Sample is representative of entire population (not yet moulting, in moult, and
birds which have completed moult). Moult scores are required.

type 3: Sample is representative only of birds in moult. Individuals
with moult scores 0 or 1 are ignored. 

type 4: Sample is representative only of birds in moult and those that
have completed moult. Individuals with moult scores 0 are ignored. 

type 5: Sample is representative only of birds that have not started
moult or that are in moult. Individuals with moult scores 1 are ignored. 

\code{moult} uses the R function \code{optim} to minimise the negative log-likelihood. 

The standard deviation parameters are estimated on the log-scale. Thus
the corresponding elements in the covariance matrix are also on the
log-scale. Standard errors for the standard deviation parameters are on
the scale of the standard deviation (in days), estimated using the delta
method.

}

\value{
    \item{coefficients}{parameter estimates split up into duration, mean
      start date and standard deviation of start date}
    \item{loglik}{log-likelihood at parameter estimates}
    \item{vcov}{variance covariance matrix for paramter estimates}
    \item{standard.errors}{vector of standard errors for parameter
      estimates, obtained from diagonal elements of vcov, see details}
    \item{type}{type of data assumed, see details above}
    \item{residuals}{vector of residuals: observed - fitted moult index}
    \item{fitted.values}{a vector of fitted values (moult scores)}
    \item{n}{number of observations}
    \item{df.residual}{residual degrees of freedom for fitted model}
    \item{terms}{duration formula, mean formula, full formula}
    \item{mean.formula}{model formula for mean start date}
    \item{duration.formula}{model formula for duration of moult}
    \item{formula}{model formula for mean start and duration of moult}
    \item{sd.formula}{grouping variable used to estimate standard
      deviations in mean start dates, different for each group}
    \item{optim}{object returned by call to optim, which minimises negative log-likelihood}
    \item{converged}{logical value indicating whether algorithm has converged or not}
    \item{convergence.value}{value for convergence returned by optim, see \code{optim} for details}
  }

\references{Underhill, L. G. and Zucchini, W. (1998). A model for avian primary moult. Ibis 130: 358-372.

Underhill, L. G. and Zucchini, W. and Summers, R. W. (1990). A model for avian primary moult-data types based on migration strategies and an example using the Redshank \emph{Tringa totanus}. Ibis 132: 118-123. }

\author{Birgit.Erni <birgit.erni@uct.ac.za>}

\examples{

data(sanderlings)

m2 <- moult(sanderlings$MIndex,sanderlings$Day) 
m2$coefficients
}

