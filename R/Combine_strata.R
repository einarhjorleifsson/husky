#' Sums mean, standard error, cv by aggregated area
#'
#' Function used both in \code{\link{Calc.index}}.
#'
#' location: /net/hafkaldi/export/u2/reikn/Splus5/SMB/GEOMETRY.NEW/.RData
#'
#' @aliases index.aggregate.combine
#' @param result Result
#' @param index Index
#'
Combine.strata <- function (result, index) {

  if (!missing(index)) result <- result[index, ]
  Mean <- sum(result$mean * result$area)/sum(result$area)
  Sum <- sum(result$mean * result$area)
  total <- sum(result$total)
  totalarea <- sum(result$area)
  i <- !is.na(result$sdev)
  tmpsum <- sum(result$mean[i] * result$area[i])
  Calc.sdev <- sqrt(sum(result$sdev[i]^2 * result$area[i]^2/result$count[i])/sum(result$area[i])^2)
  Sdev <- Calc.sdev * Sum/tmpsum
  return(data.frame(mean = Mean, se = Sdev, cv = Sdev/Mean,
                    count = sum(result$count), area = totalarea, total = total))
}
