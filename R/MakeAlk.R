#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param kvarnir XXX
#' @param tegund XXX
#' @param kyn XXX
#' @param kynth  XXX
#' @param okynth.gildi XXX
#' @param kynth.gildi XXX
#' @param lengd XXX
#' @param aldur XXX
#' @param Stodvar XXX
#' @param FilterAldurLengd XXX
#' @param NAOkynthAge XXX

MakeAlk <- function (kvarnir, tegund, kyn = F, kynth = F, okynth.gildi = 1,
                     kynth.gildi = c(2:4), lengd, aldur, Stodvar, FilterAldurLengd = T,
                     NAOkynthAge = 3) {
  medalle <- (lengd[-length(lengd)] + lengd[-1])/2
  if (missing(kvarnir))
    kvarnir <- fjolst::lesa.kvarnir(Stodvar$synis.id, tegund, c("kyn",
                                                                "kynthroski"))
  else if (is.na(match("kynthroski", names(kvarnir))))
    kvarnir$kynthroski <- rep(1, nrow(kvarnir))
  if (is.na(match("kyn", names(kvarnir))))
    kvarnir$kyn <- rep(1, nrow(kvarnir))
  kvarnir <- kvarnir[!is.na(kvarnir$aldur), ]
  kvarnir <- kvarnir[!is.na(match(kvarnir$aldur, aldur)), ]
  if (FilterAldurLengd) {
    kvarnir <- kv.filter(kvarnir, tegund)
  }
  if (kyn)
    kvarnir <- kvarnir[!is.na(match(kvarnir$kyn, c(1, 2))),
                       ]
  if (kynth) {
    i <- kvarnir$aldur <= NAOkynthAge & is.na(kvarnir$kynthroski)
    kvarnir$kynthroski[i] <- okynth.gildi
    kvarnir <- kvarnir[!is.na(kvarnir$kynthroski) & !is.na(match(kvarnir$kynthroski,
                                                                 c(kynth.gildi, okynth.gildi))), ]
    i <- !is.na(match(kvarnir$kynthroski, okynth.gildi))
    kvarnir$Kynth[i] <- 0
    i <- !is.na(match(kvarnir$kynthroski, kynth.gildi))
    kvarnir$Kynth[i] <- 1
  }
  else kvarnir$Kynth <- rep(0, nrow(kvarnir))
  if (is.na(match("fjoldi", names(kvarnir))))
    kvarnir$fjoldi <- rep(1, nrow(kvarnir))
  kvarnir$lenfl <- cut(kvarnir$lengd, breaks = lengd)
  kvarnir <- kvarnir[!is.na(match(kvarnir$aldur, aldur)) &
                       !is.na(kvarnir$lenfl), ]
  n <- length(medalle)
  n1 <- length(aldur)
  tmp.data <- data.frame(aldur = rep(aldur[1], n), lengd = medalle,
                         Kynth = rep(0, n), kyn = rep(1, n), fjoldi = rep(0, n))
  tmp.data <- rbind(tmp.data, data.frame(aldur = aldur, lengd = rep(medalle[1],
                                                                    n1), Kynth = rep(1, n1), kyn = rep(2, n1), fjoldi = rep(0,
                                                                                                                            n1)))
  tmp.data <- rbind(tmp.data, kvarnir[, names(tmp.data)])
  tmp.data$lenfl <- cut(tmp.data$lengd, breaks = lengd)
  OTH.TOT <- tapply(tmp.data$fjoldi, list(tmp.data$lenfl, tmp.data$aldur),
                    sum)
  OTH.TOT[is.na(OTH.TOT)] <- 0
  ALK.TOT <- Calc.ALK(OTH.TOT, kvarnir, aldur, lengd)
  if (kynth) {
    i <- tmp.data$Kynth == 1
    OTH.KYNTH <- tapply(tmp.data$fjoldi[i], list(tmp.data$lenfl[i],
                                                 tmp.data$aldur[i]), sum)
    OTH.KYNTH[is.na(OTH.KYNTH)] <- 0
    ALK.KYNTH <- Calc.ALK.KYNTH(OTH.KYNTH, OTH.TOT, ALK.TOT,
                                kvarnir, aldur, lengd)
  }
  if (kyn) {
    i <- tmp.data$kyn == 2
    OTH.FEMALE <- tapply(tmp.data$fjoldi[i], list(tmp.data$lenfl[i],
                                                  tmp.data$aldur[i]), sum)
    OTH.FEMALE[is.na(OTH.FEMALE)] <- 0
    ALK.FEMALE <- Calc.ALK.KYN(OTH.FEMALE, OTH.TOT, ALK.TOT,
                               kvarnir, aldur, lengd)
  }
  utkoma <- list()
  utkoma$ALK.TOT <- ALK.TOT
  utkoma$OTH.TOT <- OTH.TOT
  class(utkoma) <- c("ALK1", "list")
  if (kynth) {
    utkoma$ALK.KYNTH <- ALK.KYNTH
    utkoma$OTH.KYNTH <- OTH.KYNTH
    class(utkoma) <- c(class(utkoma), "ALK.KYNTH")
    attr(utkoma, "KYNTH") <- T
  }
  if (kyn) {
    utkoma$ALK.FEMALE <- ALK.FEMALE
    utkoma$OTH.FEMALE <- OTH.FEMALE
    utkoma$ALK.MALE <- utkoma$ALK.TOT - utkoma$ALK.FEMALE
    utkoma$OTH.MALE <- utkoma$OTH.TOT - utkoma$OTH.FEMALE
    attr(utkoma, "KYN") <- T
    class(utkoma) <- c(class(utkoma), "ALK.KYN")
  }
  attr(utkoma, "lengd") <- lengd
  attr(utkoma, "aldur") <- aldur
  attr(utkoma, "medalle") <- medalle
  attr(utkoma, "get.alk.call") <- match.call()
  return(utkoma)
}
