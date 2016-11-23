#' Title
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @param ALK XXX
#' @param LDIST XXX
#'
#' @return XXX
#' @export
#'
Calc.fj.per.station <-
  function (ALK, LDIST)
  {
    FjPerAldur <- LDIST$LDIST.ALLS %*% ALK$ALK.TOT
    LengdSinnumFjPerAldur <- LDIST$LENGD.SINNUM.FJOLDI %*% ALK$ALK.TOT
    Lengd2SinnumFjPerAldur <- LDIST$LENGD2.SINNUM.FJOLDI %*%
      ALK$ALK.TOT
    if (!is.na(match("ALK.KYNTH", names(ALK))))
      KynthFjPerAldur <- LDIST$LDIST.ALLS %*% ALK$ALK.KYNTH
    tmpwt <- matrix(LDIST$MEANWT.ALLS, nrow(LDIST$LDIST.ALLS),
                    ncol(LDIST$LDIST.ALLS), byrow = T)
    tmpwt[is.na(tmpwt)] <- 0
    WtPerAldur <- ((LDIST$LDIST.ALLS * tmpwt) %*% ALK$ALK.TOT)/FjPerAldur
    utkoma <- list(FjPerAldur = FjPerAldur, WtPerAldur = WtPerAldur)
    utkoma$MeanlePerAldur <- LengdSinnumFjPerAldur/FjPerAldur
    utkoma$LengdSinnumFjPerAldur <- LengdSinnumFjPerAldur
    utkoma$Lengd2SinnumFjPerAldur <- Lengd2SinnumFjPerAldur
    if (!is.na(match("ALK.KYNTH", names(ALK)))) {
      utkoma$KynthFjPerAldur <- LDIST$LDIST.ALLS %*% ALK$ALK.KYNTH
      utkoma$KynthWtPerAldur <- (LDIST$LDIST.ALLS * tmpwt) %*%
        ALK$ALK.KYNTH/utkoma$KynthFjPerAldur
    }
    attr(utkoma, "Stodvar") <- attr(LDIST, "Stodvar")
    attr(utkoma, "aldur") <- attr(ALK, "aldur")
    class(utkoma) <- c("aldurdre.station", "list")
    return(utkoma)
  }
