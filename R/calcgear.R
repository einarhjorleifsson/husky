#' @title XXX
#'
#' @description XXX
#'
#' @export
#'
#' @param sfile XXX
#' @param gearlist XXX
calcgear <-
  function (sfile,gearlist=pax::gearlist)
  {
    j <- match(sfile$veidarfaeri, gearlist$veidarfaeri)
    sfile$vf <- gearlist$geartext[j]
    return(sfile)
  }
