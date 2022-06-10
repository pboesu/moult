ms2pfmg <- function (ms, fm, split = "") {
  pfmg <- numeric()
  
  mscore <- c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5")
  ps <- c(0, 0.0625, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 0.9375, 1)
  
  totalmass <- sum(fm)
  j <- length(fm)
  
  for (i in 1:length(ms)) {
    if(is.numeric(ms[i]))
      scores <- format(ms[i], scientific = FALSE, trim = TRUE)
    else
      scores <- ms[i]
    
    scores <- unlist(strsplit(scores, split = split))
    if (any(scores > 5)) warning("Some moult scores > 5, will result in NA \n")
    ind <- match(scores, mscore)     # index
    ind.fm <- ps[ind]   # individual feather mass grown
    pm <- sum(ind.fm * fm) / totalmass
    pfmg <- c(pfmg, pm)
  }
  return(pfmg)
}
