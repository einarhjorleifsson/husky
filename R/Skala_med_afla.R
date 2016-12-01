#' @title XXX
#'
#' @description XXX
#'
#' location "/net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData"
#'
#' @export
#'
#' @param lengdir XXX
#' @param afli XXX
Skala.med.afla <-
  function (lengdir, afli)
  {
    rat <- afli/sum(lengdir$wt * lengdir$fjoldi/1000)
    lengdir$fj.alls <- lengdir$fjoldi * rat
    return(lengdir)
  }
