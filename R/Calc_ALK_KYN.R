#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param OTH.KYN XXX
#' @param OTH.TOT XXX
#' @param ALK.TOT XXX
#' @param kfile XXX
#' @param aldur XXX
#' @param lengd XXX
#'
Calc.ALK.KYN <-
  function (OTH.KYN, OTH.TOT, ALK.TOT, kfile, aldur, lengd)
  {
    ALK.KYN <- OTH.KYN/apply(OTH.TOT, 1, sum)
    ALK.KYN[is.na(ALK.TOT)] <- 0
    i <- apply(OTH.TOT, 1, sum) == 0
    if (any(i)) {
      tmp <- kfile[, c("Kyn", "lengd")]
      pred.data <- data.frame(Kyn = rep(0, length(lengd)),
                              lengd = lengd)
      x <- glm(Kyn ~ lengd, data = tmp, family = binomial)
      kyn.lengd <- predict(x, pred.data, type = "response")
      i1 <- 1:length(i)
      i1 <- i1[i]
      for (j in i1) ALK.KYN[j, ] <- ALK.TOT[j, ] * kyn.lengd[j]
    }
    print(length(ALK.KYN))
    return(ALK.KYN)
  }
