#' @title XXX
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @export
#'
#' @param lengdir XXX
#' @param l.th.gogn XXX
Reikna.thyngd <-
  function (lengdir, l.th.gogn)
  {
    if (is.data.frame(l.th.gogn)) {
      lengdir <- plyr::join(lengdir, l.th.gogn[, c("lengd", "oslaegt")],
                            "lengd")
      n <- ncol(lengdir)
      names(lengdir)[n] <- "wt"
      lengdir <- lengdir[!is.na(lengdir$wt), ]
    }
    else lengdir$wt <- l.th.gogn[1] * lengdir$lengd^l.th.gogn[2]
    return(lengdir)
  }
