#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param OTH.KYNTH XXX
#' @param OTH.TOT XXX
#' @param ALK.TOT XXX
#' @param kfile XXX
#' @param aldur XXX
#' @param lengd XXX

Calc.ALK.KYNTH <-
  function (OTH.KYNTH, OTH.TOT, ALK.TOT, kfile, aldur, lengd)
  {
    ALK.KYNTH <- OTH.KYNTH/apply(OTH.TOT, 1, sum)
    ALK.KYNTH[is.na(ALK.TOT)] <- 0
    i <- apply(OTH.TOT, 1, sum) == 0
    if (any(i)) {
      tmp <- kfile[, c("Kynth", "lengd")]
      pred.data <- data.frame(Kynth = rep(0, length(lengd)),
                              lengd = lengd)
      x <- glm(Kynth ~ lengd, data = tmp, family = binomial)
      kynth.lengd <- predict(x, pred.data, type = "response")
      i1 <- 1:length(i)
      i1 <- i1[i]
      for (j in i1) ALK.KYNTH[j, ] <- ALK.TOT[j, ] * kynth.lengd[j]
    }
    return(ALK.KYNTH)
  }
