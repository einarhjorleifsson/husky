#' @title XXX
#'
#' @description XXX
#'
#' @export
#'
#' @param lengdir XXX
#' @param Stodvar XXX
#' @param std.toglengd XXX
Skala.med.toglengd <-
  function (lengdir, Stodvar, std.toglengd)
  {
    if (is.na(match("fj.alls", names(lengdir)))) {
      cat("D\xe1lkur fj.alls ver\xf0ur a\xf0 vera til \xed lengdir")
      return(invisible())
    }
    i <- is.na(Stodvar$toglengd) | Stodvar$toglengd > 0
    if (any(i))
      Stodvar$toglengd[i] <- std.toglengd
    tmp <- plyr::join(lengdir, Stodvar[, c("synis.id", "toglengd")],
                      "synis.id")
    lengdir$fj.alls <- lengdir$fj.alls * std.toglengd/lengdir$toglengd
    return(lengdir)
  }
