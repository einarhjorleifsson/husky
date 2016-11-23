#' Title
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @param lengdir XXX
#' @param numer XXX
#' @param tegund XXX
#' @param kyn XXX
#' @param kynth XXX
#' @param okynth.gildi XXX
#' @param kynth.gildi XXX
#' @param lengd XXX
#' @param Stodvar XXX
#' @param talid XXX
#' @param std.toglengd XXX
#' @param lengd.thyngd.data XXX
#' @param stodvar.col XXX
#'
#' @return XXX
#' @export
#'
MakeLdistbyStation <-
  function (lengdir, numer, tegund, kyn = F, kynth = F, okynth.gildi = 1,
            kynth.gildi = 2, lengd = seq(4.5, 130.5, by = 5), Stodvar,
            talid, std.toglengd = NULL, lengd.thyngd.data, stodvar.col = c("synis.id",
                                                                           "lat", "lon", "strata", "newstrata", "area", "toglengd"))
  {
    medalle <- (lengd[-length(lengd)] + lengd[-1])/2
    if (!missing(lengd.thyngd.data))
      WT <- T
    else WT <- F
    col.names <- c("lengd", "synis.id", "fjoldi")
    if (missing(lengdir))
      lengdir <- lesa.lengdir(Stodvar$synis.id, tegund, col.names)
    else lengdir <- lengdir[!is.na(match(lengdir$synis.id, Stodvar$synis.id)),
                            ]
    col.names <- c(col.names, "fj.alls")
    lengdir <- lengdir[!is.na(cut(lengdir$lengd, lengd)), ]
    if (talid) {
      if (missing(numer))
        numer <- lesa.numer(Stodvar$synis.id, tegund)
      else numer <- numer[!is.na(match(numer$synis.id, Stodvar$synis.id)),
                          ]
      lengdir <- Skala.med.toldum(lengdir, numer)
    }
    else lengdir$fj.alls <- lengdir$fjoldi
    if (!is.null(std.toglengd))
      lengdir <- Skala.med.toglengd(lengdir, Stodvar)
    n <- length(medalle)
    m <- nrow(Stodvar)
    tmp.data <- data.frame(synis.id = Stodvar$synis.id, lengd = rep(medalle[1],
                                                                    m), fjoldi = rep(0, m), fj.alls = rep(0, m), kynth = rep(0,
                                                                                                                             m), kyn = rep(1, m), wt = rep(0, m))
    tmp.data <- rbind(tmp.data, data.frame(synis.id = rep(lengdir$synis.id[1],
                                                          n), lengd = medalle, fjoldi = rep(0, n), fj.alls = rep(0,
                                                                                                                 n), kynth = rep(0, n), kyn = rep(1, n), wt = rep(0, n)))
    for (i in 1:ncol(tmp.data)) tmp.data[, i] <- as.numeric(tmp.data[,
                                                                     i])
    tmp.data$kyn[1] <- 2
    tmp.data$kynth[1] <- 1
    tmp.data <- tmp.data[, col.names]
    tmp.data <- rbind(tmp.data, lengdir[, names(tmp.data)])
    tmp.data$lenfl <- cut(tmp.data$lengd, lengd)
    LDIST.ALLS <- tapply(tmp.data$fj.alls, list(tmp.data$synis.id,
                                                tmp.data$lenfl), sum)
    LDIST.ALLS[is.na(LDIST.ALLS)] <- 0
    utkoma <- list()
    class(utkoma) <- c("ldiststation", "list")
    utkoma$LDIST.ALLS <- LDIST.ALLS
    LENGD.SINNUM.FJOLDI <- tapply(tmp.data$fj.alls * tmp.data$lengd,
                                  list(tmp.data$synis.id, tmp.data$lenfl), sum)
    LENGD.SINNUM.FJOLDI[is.na(LENGD.SINNUM.FJOLDI)] <- 0
    utkoma$LENGD.SINNUM.FJOLDI <- LENGD.SINNUM.FJOLDI
    LENGD2.SINNUM.FJOLDI <- tapply(tmp.data$fj.alls * tmp.data$lengd^2,
                                   list(tmp.data$synis.id, tmp.data$lenfl), sum)
    LENGD2.SINNUM.FJOLDI[is.na(LENGD2.SINNUM.FJOLDI)] <- 0
    utkoma$LENGD2.SINNUM.FJOLDI <- LENGD2.SINNUM.FJOLDI
    attr(utkoma, "lengd") <- lengd
    attr(utkoma, "medalle") <- medalle
    s.id <- as.numeric(dimnames(LDIST.ALLS)[[1]])
    i <- match(s.id, Stodvar$synis.id)
    attr(utkoma, "Stodvar") <- Stodvar[i, stodvar.col]
    if (WT) {
      tmp.data <- Reikna.thyngd(tmp.data, lengd.thyngd.data)
      MEANWT.ALLS <- tapply(tmp.data$fj.alls * tmp.data$wt,
                            tmp.data$lenfl, sum)/tapply(tmp.data$fj.alls, tmp.data$lenfl,
                                                        sum)
      utkoma$MEANWT.ALLS <- MEANWT.ALLS
      class(utkoma) <- c(class(utkoma), "weight")
    }
    return(utkoma)
  }
