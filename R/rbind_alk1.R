#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param fj1 XXX
#' @param fj2 XXX
rbind.alk1 <-
  function (fj1, fj2)
  {
    if (length(fj1) == 0)
      return(fj2)
    else if (length(fj2) == 0)
      return(fj1)
    else {
      fj3 <- list()
      fj3$FjPerAldur <- fj1$FjPerAldur + fj2$FjPerAldur
      fj3$BiomassPerAldur <- fj1$BiomassPerAldur + fj2$BiomassPerAldur
      fj3$ALD <- fj1$ALD + fj2$ALD
      fj3$WtPerAldur <- fj3$BiomassPerAldur/fj3$FjPerAldur
      if (!is.na(match("MeanlePerAldur", names(fj1)))) {
        tmp <- (fj1$MeanlePerAldur * fj1$FjPerAldur)
        tmp[is.na(tmp)] <- 0
        tmp1 <- (fj2$MeanlePerAldur * fj2$FjPerAldur)
        tmp1[is.na(tmp1)] <- 0
        fj3$MeanlePerAldur <- (tmp + tmp1)/fj3$FjPerAldur
      }
      if (!is.na(match("FjKynthPerAldur", names(fj1))))
        fj3$FjKynthPerAldur <- fj1$FjKynthPerAldur + fj2$FjKynthPerAldur
      if (!is.na(match("KynthHlutfall", names(fj1))))
        fj3$KynthHlutfall <- (fj1$KynthHlutfall * fj1$FjPerAldur +
                                fj2$KynthHlutfall * fj2$FjPerAldur)/fj3$FjPerAldur
      if (!is.na(match("ALD.KYNTH", names(fj1))))
        fj3$ALD.KYNTH <- fj1$ALD.KYNTH + fj2$ALD.KYNTH
      if (!is.na(match("MeanleKynthPerAldur", names(fj1)))) {
        tmp <- (fj1$MeanleKynthPerAldur * fj1$FjKynthPerAldur)
        tmp[is.na(tmp)] <- 0
        tmp1 <- (fj2$MeanleKynthPerAldur * fj2$FjKynthPerAldur)
        tmp1[is.na(tmp1)] <- 0
        fj3$MeanlePerAldur <- (tmp + tmp1)/fj3$FjKynthPerAldur
      }
      if (!is.na(match("BiomassKynthPerAldur", names(fj2)))) {
        fj3$BiomassKynthPerAldur <- fj1$BiomassKynthPerAldur +
          fj2$BiomassKynthPerAldur
        fj3$WtKynthPerAldur <- fj3$BiomassKynthPerAldur/fj3$FjKynthPerAldur
      }
      return(fj3)
    }
  }
