day2date <- function (dayno, startmonth)
{
     dno <- 0
     mm <- 1
     dayno <- round(dayno)
     ndays <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
     nmon <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
         "Aug", "Sep", "Oct", "Nov", "Dec")
     if (startmonth != 1) {
         for (i in 1:(startmonth - 1)) {
             dno <- dno + ndays[i]
         }
	}
     dno <- dno + dayno
	if (dno > 365) {
	    dno <- dno - 365
     }
     while (dno > (ndays[mm])) {
         dno <- dno - ndays[mm]
         mm <- mm + 1
     }
     print(paste(attributes(dno), "is", dno, nmon[mm]), quote = FALSE)
}
