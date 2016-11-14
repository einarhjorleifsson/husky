#' @title XXX
#'
#' @description XXX
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
