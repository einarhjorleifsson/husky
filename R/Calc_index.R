#' Calc.index
#'
#' @description XXX
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/GETALK/.RData
#'
#' @param sfile XXX
#' @param colname XXX
#' @param strata.list XXX
#' @param std.toglengd XXX
#' @param trollbreidd XXX
#' @param combine.output XXX
#' @param use.rallarea XXX
#' @param leidretta.fyrir.toglengd XXX
#' @param z XXX
#' @param std.cv XXX
#' @param mintoglengd XXX
#' @param maxtoglengd XXX
#'
#' @export
#'
Calc.index <-
  function (sfile, colname, strata.list, std.toglengd = 4, trollbreidd = 17,
            combine.output = Std.aggregation, use.rallarea = T, leidretta.fyrir.toglengd = T,
            z, std.cv = 1, mintoglengd, maxtoglengd)
  {
    sfile$strata <- sfile$newstrata
    if (!missing(z)) {
      sfile$outcome <- z
      colname <- "outcome"
    }
    if (missing(mintoglengd))
      mintoglengd <- std.toglengd/2
    if (missing(maxtoglengd))
      maxtoglengd <- std.toglengd * 2
    std.area <- std.toglengd * trollbreidd/1852
    if (!is.na(match("toglengd", names(sfile))) && leidretta.fyrir.toglengd) {
      i <- is.na(sfile$toglengd)
      if (any(i))
        sfile$toglengd[i] <- std.toglengd
      i <- sfile$toglengd < mintoglengd
      if (any(i))
        sfile$toglengd[i] <- mintoglengd
      i <- sfile$toglengd > maxtoglengd
      if (any(i))
        sfile$toglengd[i] <- maxtoglengd
      sfile[, colname] <- sfile[, colname] * std.toglengd/sfile$toglengd
    }
    if (use.rallarea)
      areas <- attributes(STRATAS)$rall.area
    else areas <- attributes(STRATAS)$area
    Names <- attributes(STRATAS)$name
    if (!missing(strata.list)) {
      for (i in 1:length(strata.list)) {
        areas[strata.list[[i[1]]]] <- sum(areas[strata.list[[i]]])
        j <- !is.na(match(sfile$strata, strata.list[[i]]))
        if (any(j))
          sfile$strata[j] <- strata.list[[i]][1]
      }
    }
    tmp6 <- apply.shrink(sfile[, colname], sfile$strata, mean)
    tmp7 <- apply.shrink(sfile[, colname], sfile$strata, sdev)
    tmp8 <- apply.shrink(rep(1, nrow(sfile)), sfile$strata, sum)
    result <- data.frame(strata = tmp8[, 1], mean = tmp6[, 2],
                         sdev = tmp7[, 2], count = tmp8[, 2])
    names(result) <- c("strata", "mean", "sdev", "count")
    i <- result$count == 1
    if (any(i))
      result$sdev[i] <- result$mean[i] * std.cv
    result$area <- areas[result$strata]/1.854^2
    result$se <- result$sdev/sqrt(result$count)
    result$cv <- result$se/result$mean
    result$total <- result$mean * result$area/std.area
    Res.names <- Names[result$strata]
    aggr.output <- data.frame(matrix(0, length(combine.output),
                                     6))
    names(aggr.output) <- c("mean", "se", "cv", "count", "area",
                            "total")
    row.names(aggr.output) <- names(combine.output)
    for (i in 1:length(combine.output)) {
      j <- !is.na(match(result$strata, combine.output[[i]]))
      j1 <- c(1:length(j))
      j1 <- j1[j]
      if (length(j1) > 0)
        aggr.output[i, ] <- Combine.strata(result, j1)
    }
    aggr.output$area <- round(aggr.output$area)
    aggr.output$mean <- round(aggr.output$mean, 3)
    aggr.output$se <- round(aggr.output$se, 3)
    aggr.output$cv <- round(aggr.output$cv, 3)
    aggr.output$se <- round(aggr.output$se, 4)
    aggr.output$total <- round(aggr.output$total, 1)
    result$area <- round(result$area)
    result$mean <- round(result$mean, 3)
    result$sdev <- round(result$sdev, 3)
    result$cv <- round(result$cv, 3)
    result$se <- round(result$se, 4)
    result$total <- round(result$total, 1)
    return(list(result = result, Res.names = Res.names, aggr.output = aggr.output))
  }
