## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(moult)

data(sanderlings)
m2 <- moult(MIndex ~ Day, data = sanderlings) 
summary(m2)

## ---- fig.height = 8, fig.width = 8, warning = FALSE--------------------------
set.seed(123)
boot.ci <- confint(m2, B = 1000, pch = 19, cex = 0.3, las = 1)
boot.ci$bootstrap.percentile.ci
boot.ci$bootstrap.vcov
boot.ci$bootstrap.SE

