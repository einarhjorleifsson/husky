
join <-
  function (x, y, column, name.x, name.y, outer.join = T, set = NA)
  {
    if (!missing(name.x)) {
      if (!missing(column))
        name.x <- unique(c(name.x, column))
      x <- x[, name.x]
    }
    if (!missing(name.y)) {
      if (!missing(column))
        name.y <- unique(c(name.y, column))
      y <- y[, name.y]
    }
    if (missing(column)) {
      d1 <- row.names(x)
      d2 <- row.names(y)
    }
    else {
      d1 <- x[, column]
      ind <- match(names(y), column)
      ind <- ind[!is.na(ind)]
      if (length(ind) == 0)
        d2 <- row.names(y)
      else d2 <- y[, column]
    }
    if (outer.join) {
      y0 <- matrix(set, nrow(x), ncol(y))
      y1 <- as.data.frame(y0)
      for (i in 1:ncol(y0)) y1[, i] <- I(y0[, i])
      names(y1) <- names(y)
      ind <- match(d1, d2)
      index <- c(1:length(ind))
      index <- index[!is.na(ind)]
      ind <- ind[index]
      y1[index, ] <- y[ind, ]
      outcome <- cbind(x, y1)
    }
    else {
      ind <- match(d1, d2)
      ind1 <- c(1:length(ind))
      ind1 <- ind1[!is.na(ind)]
      ind <- ind[!is.na(ind)]
      x <- (x[ind1, ])
      y <- (y[ind, ])
      outcome <- cbind(x, y)
    }
    if (!missing(column)) {
      i <- match(names(outcome), paste(column, "1", sep = ""))
      i1 <- c(1:length(i))
      i1 <- i1[!is.na(i)]
      if (length(i1) > 0) {
        outcome <- outcome[, -i1]
      }
      i <- match(names(outcome), column)
      i1 <- c(1:length(i))
      i1 <- i1[!is.na(i)]
      if (length(i1) > 1) {
        i1 <- i1[-1]
        outcome <- outcome[, -i1]
      }
    }
    return(outcome)
  }
