#' @title XXX
#'
#' @description XXX
#'
#' @export
#'
#' @param object XXX
#' @param x XXX
#' @param deriv XXX
predict.smooth.spline <-
  function (object, x, deriv = 0)
  {
    if (missing(x)) {
      if (deriv == 0)
        return(object[c("x", "y")])
      else x <- object$x
    }
    fit <- object$fit
    if (is.null(fit))
      stop("not a valid smooth.spline object")
    else predict(fit, x, deriv)
  }
