\name{confint.moult}
\alias{confint.moult}
\title{Conficence Intervals for Moult Parameters}

\description{
Calculate normal theory and bootstrap percentile intervals for moult parameters. Also returns bootstrap covariance matrix and standard error estimates for parameters. }
\usage{
\method{confint}{moult}(object, parm, level = 0.95, ..., B = 1000, add.plot = TRUE)
}

\arguments{
  \item{object}{model output returned by call to moult.}
  \item{parm}{currently ignored (should be a vector of numbers or names, specifying for which parameters to calculate confidence intervals).}
  \item{level}{confidence level.}
  \item{\dots}{additional arguments for plot.}
  \item{B}{integer: number of bootstrap samples to take.}
  \item{add.plot}{logical: should scatterplot of bootstrapped parameters be added?}
}

\details{End date is calculated as \emph{mean start date + duration}. Half-date is calculated as \emph{mean start date + 0.5 * duration}.
}
\value{Returns a list with elements:
\item{bootstrap.distribution}{\eqn{B \times p} matrix of bootstrapped parameters.}
\item{bootstrap.percentile.ci}{Table of bootstrap percentile intervals for all parameters. In addition, confidence intervals for half-date and end date (see details) are given.}
\item{bootstrap.vcov}{Bootstrap estimates for variances and covariances between all parameters.}
\item{bootstrap.SE}{Bootstrap estimates of parameter standard errors.}
}
\author{Birgit Erni \email{birgit.erni@uct.ac.za}}

\examples{
data(sanderlings)
m2 <- moult(MIndex ~ Day, data = sanderlings)
confint(m2, B = 100) # increase B for better precision
}