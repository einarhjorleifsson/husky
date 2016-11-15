apply.shrink <-
  function (X, INDICES, FUN = NULL, names, ...)
  {
    if (missing(FUN))
      stop("No function to apply to data given (missing argument FUN)")
    if (!is.list(INDICES))
      INDICES <- list(INDICES)
    len.data <- length(X)
    all.indices <- rep(0, len.data)
    for (i in rev(INDICES)) {
      if (length(i) != len.data)
        stop("Data and all indices must have same length")
      i <- as.factor(i)
      all.indices <- all.indices * length(levels(i)) + (as.vector(unclass(i)) -
                                                          1)
    }
    all.indices <- all.indices + 1
    INDICES <- as.data.frame(INDICES)
    INDICES <- INDICES[match(sort(unique(all.indices)), all.indices,
                             nomatch = 0), ]
    if (is.character(FUN))
      FUN <- getFunction(FUN)
    else if (mode(FUN) != "function") {
      farg <- substitute(FUN)
      if (mode(farg) == "name")
        FUN <- getFunction(farg)
      else stop(paste("\"", farg, "\" is not a function", sep = ""))
    }
    X <- split(X, all.indices)
    X.apply <- lapply(X, FUN, ...)
    numb.FUN.value <- length(X.apply[[1]])
    if (numb.FUN.value == 1)
      X.apply <- data.frame(X = unlist(X.apply))
    else X.apply <- data.frame(matrix(unlist(X.apply), ncol = numb.FUN.value,
                                      byrow = T, dimnames = list(NULL, names(X.apply[[1]]))))
    X.apply <- cbind(INDICES, X.apply)
    if (!missing(names))
      names(X.apply) <- names
    return(X.apply)
  }
