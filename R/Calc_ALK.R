#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param OTH.TOT XXX
#' @param kfile XXX
#' @param aldur XXX
#' @param lengd XXX

Calc.ALK <-
  function (OTH.TOT, kfile, aldur, lengd)
  {
    ALK.TOT <- OTH.TOT/apply(OTH.TOT, 1, sum)
    ALK.TOT[is.na(ALK.TOT)] <- 0
    i <- apply(ALK.TOT, 1, sum) == 0
    if (any(i)) {
      i1 <- 1:length(i)
      i1 <- i1[i]
      x <- smooth.spline(kfile$lengd, kfile$aldur, df = 2)
      meana <- round(predict.smooth.spline(x, lengd[i])$y)
      i <- meana > max(aldur)
      if (any(i))
        meana[i] <- max(aldur)
      i <- meana < min(aldur)
      if (any(i))
        meana[i] <- min(aldur)
      for (j in 1:length(i1)) ALK.TOT[i1[j], as.character(meana[j])] <- 1
    }
    return(ALK.TOT)
  }
