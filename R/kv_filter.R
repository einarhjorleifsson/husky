#' @title XXX
#'
#' @description XXX
#'
#' location: file:/net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param kv XXX
#' @param tegund XXX
#' @param maxage XXX
#' @param vikmork XXX
kv.filter <-
  function (kv, tegund, maxage, vikmork=pax::vikmork)
  {
    all.teg <- c(1, 2, 3)
    if (is.na(match(tegund, all.teg))) {
      print("G\xf6gn ekki til fyrir tegund")
      return(kv)
    }
    minl <- vikmork[[as.character(tegund)]]$minl
    maxl <- vikmork[[as.character(tegund)]]$maxl
    kv <- kv[!is.na(kv$aldur) & !is.na(kv$lengd), ]
    nullgr <- kv[kv$aldur == 0 & kv$lengd < maxl[1], ]
    kv1 <- kv[1:2, ]
    for (i in 1:length(maxl)) {
      tmp <- kv[kv$aldur == i, ]
      if (nrow(tmp) > 0)
        tmp <- tmp[tmp$lengd > minl[i] & tmp$lengd < maxl[i],
                   ]
      if (nrow(tmp) > 0)
        kv1 <- rbind(kv1, tmp)
    }
    if (nrow(nullgr) > 0)
      kv1 <- rbind(kv1, nullgr)
    kv1 <- kv1[-c(1:2), ]
    return(kv1)
  }
