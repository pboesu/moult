\name{dfbeta.moult}
\alias{dfbeta.moult}

\title{Influence Statistics for Moult Parameters}

\description{Calculates dfbeta (change in coefficients) and dfbetas (scaled by standard error) for moult parameters.}

\usage{
\method{dfbeta}{moult}(model, ...)
}

\arguments{
  \item{model}{a model object returned by \code{\link{moult}}.}
  \item{\dots}{further arguments.}
}
\details{
Both dfbeta (absolute change in coefficients) and dfbetas (change in coefficients scaled by standard error of coefficient) are returned. 
\deqn{dfbetas_i = \frac{\hat{b} - \hat{b}_i}{SE(\hat{b}_i}),}

where the \eqn{\hat{b}_i} estimate is obtained with observation \eqn{i} removed. 

In the optional plot of dfbetas, cutoff lines at \eqn{\pm 2 / \sqrt{n}} are added. These are the limits used in linear regression models, but there is no reason that the same limits are valid for moult models. Therefore these cutoff lines can only be used as very rough guidelines.
}

\value{
\item{dfbeta}{\eqn{n \times p} matrix with absolute change in coefficients when observation \eqn{i} is deleted. }
\item{dfbetas}{\eqn{n \times p} matrix with scaled change in coefficients when observation \eqn{i} is deleted.}
}

\seealso{\code{\link{dfbeta.lm}}}

\references{Fox, J.D. (2020). \emph{Regression Diagnostics: an Introduction.} 2nd edition. SAGE Publications.}


\author{Birgit Erni \email{birgit.erni@uct.ac.za}}


\examples{
data(sanderlings)
m2 <- moult(MIndex ~ Day, data = sanderlings) 
dfbeta(m2)
}

