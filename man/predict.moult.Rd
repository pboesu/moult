\name{predict.moult}
\alias{predict.moult}
\alias{coef.moult}
\title{Predict method for moult Models}

\description{Will predict either the proportion of birds in a certain moult 
             stage (as defined in intervals) on a specified day, or the average 
             moult score on a specified day}

\usage{
  \method{predict}{moult}(object, newdata, predict.type = "prob",
  intervals = 0.1, ...)}

\arguments{
  \item{object}{moult model objects}
  \item{newdata}{optional dataframe with explanatory variables for which to 
                 make predictions.  Note the first column must contain the days 
                 (as used when fitting) for which to make predictions.}
  \item{predict.type}{specifies form of predictions, see details.}
  \item{intervals}{length of moult categories for
    which probabilities/proportions should be calculated.  The default (= 0.1) 
    will calculate the probability that a bird will fall in this moult category
    on the specified day for each of the following categories: 0, 0 to < 0.1, 
    0.1 to < 0.2, 0.2 to < 0.3, ..., 0.9 to < 1.0, 1.}
  \item{...}{other arguments passed to the predict method}
}

\details{\code{predict.type} has the following options: 

  \describe{
    \item{"response"}{the average moult index (proportion of
      feather mass grown) for each of the days specified.}
    \item{"prob"}{default, the proportion of birds in
      each of the moult categories as defined by intervals
      is predicted. }
    \item{"start"}{predicts the mean start date of moult
      for the covariate combination defined in \code{newdata},
      with standard errors. }
  }
}  
\value{
  \item{M}{if \code{predict.type = TRUE}, or \code{newdata} is missing, the
    expected moult scores at the specified days are returned,
    else a matrix of predicted probabilities for being in each of the
    moult categories defined by intervals.}
}

\references{Erni, B., Bonnevie, B. T., Oschadleus, H.-D., Altwegg,
        R. and Underhill, L. G. (2013) moult: An R package to analyze
        moult in birds. Journal of Statistical Software, \bold{52}(8), 1--23.
        \url{http://www.jstatsoft.org/v52/i08/}

	Underhill, L. G. and Zucchini, W. (1988) A model for avian primary 
        moult. \emph{Ibis} \bold{130}, 358--372.

	Underhill, L. G. and Zucchini, W. and Summers, R. W. (1990) A
	model for avian primary moult-data types based on migration
	strategies and an example using the Redshank \emph{Tringa
	totanus}. \emph{Ibis} \bold{132}, 118--123. }

\author{Birgit Erni <birgit.erni@uct.ac.za>}

\seealso{\code{\link{moult}}}
\examples{

data(sanderlings)

m2 <- moult(MIndex ~ Day, data = sanderlings)                     
m3 <- moult(MIndex ~ Day, data = sanderlings, type = 3)           

## number of birds observed on each of the days
(nn <- as.numeric(table(sanderlings$Day)))

## unique observation days 
(day <- unique(sanderlings$Day))

## Table 6 in Underhill and Zucchini 1988
## expected number of birds in each moult stage on each day
p1 <- predict(m2, newdata = data.frame(day))        
p1$M * nn						

## Table 7 in Underhill and Zucchini 1988
## expected number of birds in each moult stage for days 70 to 310
days2 <- seq(70, 310, by = 10)
p2 <- predict(m2, newdata = data.frame(days2))
p2$M * 100        
                    
(p3 <- predict(m3, newdata = data.frame(day)))

}
