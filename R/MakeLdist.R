#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param tegund XXX
#' @param kyn XXX
#' @param kynth XXX
#' @param okynth.gildi XXX
#' @param kynth.gildi XXX
#' @param lengd XXX
#' @param Stodvar XXX
#' @param lengd.thyngd.data XXX
#' @param lengdir XXX
#' @param numer XXX
#' @param talid XXX
#' @param afli XXX
#' @param std.toglengd XXX
#' @param trollbreidd XXX
#' @param area XXX

MakeLdist <-
  function (tegund, kyn = F, kynth = F, okynth.gildi = 1, kynth.gildi = 2,
            lengd, Stodvar, lengd.thyngd.data, lengdir, numer, talid = T,
            afli = NULL, std.toglengd = NULL,
            trollbreidd, area)  # EH added
  {
    medalle <- (lengd[-length(lengd)] + lengd[-1])/2
    col.names <- c("lengd", "synis.id", "fjoldi")
    if (kynth)
      col.names <- c(col.names, "kynthroski")
    if (kyn)
      col.names <- c(col.names, "kyn")
    if (!missing(lengd.thyngd.data))
      WT <- T
    else WT <- F
    if (missing(lengdir))
      lengdir <- fjolst::lesa.lengdir(Stodvar$synis.id, tegund, col.names)
    else if (!missing(Stodvar))
      lengdir <- lengdir[!is.na(match(lengdir$synis.id, Stodvar$synis.id)),
                         ]
    lengdir <- lengdir[!is.na(cut(lengdir$lengd, lengd)), ]
    if (talid) {
      if (missing(numer))
        numer <- fjolst::lesa.numer(Stodvar$synis.id, tegund)
      lengdir <- fjolst::Skala.med.toldum(lengdir, numer)
    }
    else lengdir$fj.alls <- lengdir$fjoldi
    if (is.na(match("kyn", names(lengdir))))
      lengdir$kyn <- rep(1, nrow(lengdir))
    if (is.na(match("kynth", names(lengdir))))
      lengdir$kynth <- rep(1, nrow(lengdir))
    lengdir$wt <- rep(1, nrow(lengdir))
    if (!is.null(afli)) {
      lengdir <- Reikna.thyngd(lengdir, lengd.thyngd.data)
      lengdir <- Skala.med.afla(lengdir, afli)
    }
    if (!is.null(std.toglengd))
      lengdir <- Skala.med.toglengd(lengdir, trollbreidd, area)
    col.names <- c(col.names, "fj.alls")
    tmp.data <- data.frame(synis.id = rep(lengdir$synis.id[1],
                                          length(medalle)), lengd = medalle, fjoldi = rep(0, length(medalle)),
                           fj.alls = rep(0, length(medalle)), kynth = rep(0, length(medalle)),
                           kyn = rep(1, length(medalle)), wt = rep(0, length(medalle)))
    tmp.data$kyn[1] <- 2
    tmp.data$kynth[1] <- 1
    tmp.data <- tmp.data[, col.names]
    tmp.data <- rbind(tmp.data, lengdir[, names(tmp.data)])
    tmp.data$lenfl <- cut(tmp.data$lengd, lengd)
    tmp.data <- tmp.data[!is.na(tmp.data$lenfl), ]
    utkoma <- list()
    if (WT) {
      tmp.data <- Reikna.thyngd(tmp.data, lengd.thyngd.data)
      attr(utkoma, "WT") <- T
    }
    else attr(utkoma, "WT") <- F
    Calc.tot.ldist <- function(tmp.data) {
      utkoma <- list()
      utkoma$LDIST.MAELT <- as.vector(tapply(tmp.data$fjoldi,
                                             tmp.data$lenfl, sum))
      utkoma$LDIST.ALLS <- as.vector(tapply(tmp.data$fj.alls,
                                            tmp.data$lenfl, sum))
      utkoma$MEANWT.ALLS <- as.vector(tapply(tmp.data$fj.alls *
                                               tmp.data$wt, tmp.data$lenfl, sum)/tapply(tmp.data$fj.alls,
                                                                                        tmp.data$lenfl, sum))
      utkoma$MEANLE <- as.vector(tapply(tmp.data$fj.alls *
                                          tmp.data$lengd, tmp.data$lenfl, sum)/tapply(tmp.data$fj.alls,
                                                                                      tmp.data$lenfl, sum))
      return(utkoma)
    }
    utkoma <- as.vector(Calc.tot.ldist(tmp.data))
    if (kyn) {
      i <- tmp.data$kyn == 1
      utkoma$LDIST.MALE <- as.vector(Calc.tot.ldist(tmp.data[i,
                                                             ]))
      i <- tmp.data$kyn == 2
      utkoma$LDIST.FEMALE <- as.vector(Calc.tot.ldist(tmp.data[i,
                                                               ]))
    }
    if (kyn && kynth)
      class(utkoma) <- c("ldist", "list")
    attr(utkoma, "medalle") <- medalle
    attr(utkoma, "get.alk.call") <- match.call()
    return(utkoma)
  }

