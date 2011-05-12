\name{predict.moult}
\alias{predict.moult}
\title{Predictions from Moult Models}

\description{Will make predictions on the number of birds in a certain moult stage for a given set of parameters}

\usage{
  \method{predict}{moult}(object, newdata, what = "distrib", intervals = 0.1, ...)}

\arguments{
  \item{object}{moult model objects}
  \item{newdata}{optional dataframe with explanatory variables for which to make predictions. Note the first column must contain the days (as used when fitting) for which to make predictions.}
  \item{what}{optional argument for form of predictions. The options are "distrib" (default) for predicted probabilities in categories defined by intervals, "response" for predicted moult index (proportion of feather mass grown) on specified day}
  \item{intervals}{optional argument for moult categories lengths for
    which probabilities should be calculated, default = 0.1 will
    calculate probabilites in for moult categories 0.1 to 0.2, 0.2 to
   0.3 etc.}
 \item{...}{other arguments passed to the predict method}
}

\value{
  \item{M}{If what = "response", or newdata is missing, the fitted values (predicted moult scores at the specified days), else a matrix of predicted probabilities for being in each of the moult categories defined by intervals.}
}

\references{Underhill, L. G. and Zucchini, W. (1998). A model for avian primary moult. Ibis 130: 358-372.

Underhill, L. G. and Zucchini, W. and Summers, R. W. (1990). A model for avian primary moult-data types based on migration strategies and an example using the Redshank Tringa totanus. Ibis 132: 118-123. }

\author{Birgit Erni <birgit.erni@uct.ac.za>}

\seealso{\code{moult}}
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

					# ----- extract how many birds observed on each of the days
nn <- as.numeric(tapply(sanderlings$Day, sanderlings$Day, length))        
nn
					# ----- extract days of observations
day <- unique(sanderlings$Day)                                            
day

p1 <- predict(m2, newdata=data.frame(day))        ## probabilities of moult stages 
p1$M*nn						## see Table 6 in Underhill and Zucchini 1988

days2 <- seq(70,310,by=10)
p2 <- predict(m2, newdata = data.frame(days2))
p2$M*100                                  ## see Table 7 in Underhill and Zucchini 1988


p3 <- predict(m3,newdata=data.frame(day))      
p3

}