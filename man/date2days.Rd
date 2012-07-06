\name{date2days}
\alias{date2days}
\title{Convert date to day count}
\description{Convert date mm/dd/yyyy to days since the 1st of
  startmonth, starting with days = 1 for the 1st day of startmonth.}

\usage{date2days(date, dateformat, startmonth)
}

\arguments{
  \item{date}{date character string with format as specified in dateformat}
  \item{dateformat}{string specifying format of date, e.g. "mm/dd/yyyy",
    "yyyy-dd-mm", "dd--mm" etc.}
  \item{startmonth}{integer between 1 and 12 indicating month from which
    to start counting.}}

\value{Returns an integer = number of days since 1st of startmonth,
  starting with 1 = 1st of startmonth. }

\author{Bo T. Bonnevie}

\examples{
date2days("01/01/2010", "dd/mm/yyyy", 1)
date2days("01-01-2010", "dd-mm-yyyy", 2)
date2days("2008/06/01", "yyyy/mm/dd", 8)    # year has no effect
}
