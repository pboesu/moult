\name{day2date}
\alias{day2date}
\title{Convert a days-since-count into a date of format c("dd", "month")} 
\description{Convert a count (number of days since beginning of 
             startmonth) into a date of format c("dd", "month")} 
\usage{
day2date(dayno, startmonth)
}
\arguments{
  \item{dayno}{days-since count}
  \item{startmonth}{integer between 1 and 12 indicating month from which to start 
                    counting}}
\details{Counting starts with 1 on the first day of startmonth.}
\value{Returns a date of format c("dd", "month")} 
\author{Bo Bonnevie}
\seealso{\code{date2days}}

\examples{

day2date(31,1)
day2date(31,2)
day2date(365,1)         

}