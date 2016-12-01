#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param ALK XXX
#' @param LDIST XXX

Calc.fj <-
  function (ALK, LDIST)
  {
    utkoma <- list()
    utkoma$FjPerAldur <- LDIST$LDIST.ALLS %*% ALK$ALK.TOT
    utkoma$ALD <- ALK$ALK.TOT * LDIST$LDIST.ALLS
    i <- is.na(LDIST$MEANWT.ALLS)
    if (any(i))
      LDIST$MEANWT.ALLS[i] <- 0
    utkoma$BiomassPerAldur <- (LDIST$LDIST.ALLS * LDIST$MEANWT.ALLS) %*%
      ALK$ALK.TOT
    utkoma$WtPerAldur <- utkoma$BiomassPerAldur/utkoma$FjPerAldur
    tmp <- (LDIST$LDIST.ALLS * LDIST$MEANLE)
    tmp[is.na(tmp)] <- 0
    utkoma$MeanlePerAldur <- tmp %*% ALK$ALK.TOT/utkoma$FjPerAldur
    if (!is.na(match("KYNTH", names(attributes(ALK))))) {
      utkoma$FjKynthPerAldur <- LDIST$LDIST.ALLS %*% ALK$ALK.KYNTH
      utkoma$KynthHlutfall <- utkoma$FjKynthPerAldur/utkoma$FjPerAldur *
        100
      utkoma$ALD.KYNTH <- ALK$ALK.KYNTH * LDIST$LDIST.ALLS
      utkoma$BiomassKynthPerAldur <- (LDIST$LDIST.ALLS * LDIST$MEANWT.ALLS) %*%
        ALK$ALK.KYNTH
      tmp <- (LDIST$LDIST.ALLS * LDIST$MEANLE)
      tmp[is.na(tmp)] <- 0
      utkoma$MeanleKynthPerAldur <- tmp %*% ALK$ALK.KYNTH/utkoma$FjKynthPerAldur
      utkoma$WtKynthPerAldur <- utkoma$BiomassKynthPerAldur/utkoma$FjKynthPerAldur
    }
    if (!is.na(match("KYN", names(attributes(ALK))))) {
      utkoma$FjKynPerAldur <- LDIST$LDIST.ALLS %*% ALK$ALK.KYN
      utkoma$KynHlutfall <- utkoma$FjKynPerAldur/utkoma$FjPerAldur *
        100
      utkoma$ALD.KYN <- ALK$ALK.KYN * LDIST$LDIST.ALLS
      utkoma$BiomassKynPerAldur <- (LDIST$LDIST.ALLS * LDIST$MEANWT.ALLS) %*%
        ALK$ALK.KYN
      # utkoma$WtKynPerAldur <- BiomassKynPerAldur/FjKynPerAldur # Einar: should be:
      utkoma$WtKynPerAldur <- utkoma$BiomassKynPerAldur/utkoma$FjKynPerAldur
    }
    return(utkoma)
  }
