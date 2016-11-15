#' @title Calculate length based indices
#'
#' @description Calculates abundance and biomass survey indices based on length
#' classes for a particular species in a given year.
#'
#' Not really a husky function - but a script turned into a function
#'
#' @param lengths \emph{data.frame} containing length measurement for a given
#' species in a given survey year.
#' @param stations \emph{data.frame} containg station survey station informations.
#' @param sex \emph{boolean} A flag indicating if indices should also be
#' compiled by sex
#' @param lenClass \emph{numerical vector} containing length classes to
#' compile the indices for
#' @param yr \emph{numerical} value for year
#' @param lwcoeff \emph{numerical vector} of length 2, containing parameter
#' a and b of the length weight relationship.
#' @param ... Other parameters to pass to Calc.index

bioIndex <- function(lengths,stations,sex,lenClass,yr,lwcoeff, ...) {
  for( i in 1:length(lenClass)) {
    L <- lengths

    # A.1 total number of fish of a length class by station
    le <- L[L$lengd==lenClass[i],]
    if(nrow(le) > 0) {
      x <- apply.shrink(le$fj.alls,le$synis.id,sum)
      names(x) <- c("synis.id","fj")
      st <- join(stations[,c("newstrata","toglengd","synis.id")],x,"synis.id",set=0)
    } else {
      st <- stations[,c("newstrata","toglengd","synis.id")]
      st$fj <- 0
    }
    # A.2 total number of fish of a length class by station and sex
    if(sex) {
      le <- L[L$lengd==lenClass[i] & L$kyn==1 & !is.na(L$kyn),]
      if(nrow(le) > 0) {
        x <- apply.shrink(le$fj.alls,le$synis.id,sum)
        names(x) <- c("synis.id","fjhaenga")
        st <- join(st,x,"synis.id",set=0)
      } else {
        st$fjhaenga <- 0
      }
      le <- L[L$lengd==lenClass[i] & L$kyn==2 & !is.na(L$kyn),]
      if(nrow(le) > 0) {
        x <- apply.shrink(le$fj.alls,le$synis.id,sum)
        names(x) <- c("synis.id","fjhrygna")
        st <- join(st,x,"synis.id",set=0)
      } else {
        st$fjhrygna <- 0
      }
    }

    # B. Biomass GREATER than a certain length class
    le <- L[L$lengd >= lenClass[i],]
    if(nrow(le) > 0) {
      x <- apply.shrink(le$bio,le$synis.id,sum)
      names(x) <- c("synis.id","bioge")
      st <- join(st,x,"synis.id",set=0)
    } else {
      st$bioge <- 0
    }

    # C. Abundance LESS than a certain length class
    le <- L[L$lengd <= lenClass[i],]
    if(nrow(le) > 0) {
      x <- apply.shrink(le$fj.alls,le$synis.id,sum)
      names(x) <- c("synis.id","fjle")
      st <- join(st,x,"synis.id",set=0)
    } else {
      st$fjle <- 0
    }

    # D. Calculate the indices
    visit <- Calc.index(st,"fj")
    B0 <- Calc.index(st,"bioge")
    N0 <- Calc.index(st,"fjle")
    if(sex) {
      N0sex1 <- Calc.index(st,"fjhaenga")
      N0sex2 <- Calc.index(st,"fjhrygna")
    }

    base0 <- visit$result[,c("strata","total","cv")]
    names(base0)[2:3] <- c("fj","cv.fj")
    base0$bio.staerri <- B0$result[,"total"]
    base0$cv.bio.staerri <- B0$result[,"cv"]
    base0$fj.minni <- N0$result[,"total"]
    base0$cv.fj.minni <- N0$result[,"cv"]

    if(sex) {
      base0$fj.haenga <- N0sex1$result[,"total"]
      base0$cv.fj.haenga <- N0sex1$result[,"cv"]
      base0$fj.hrygna <- N0sex2$result[,"total"]
      base0$cv.fj.hrygna <- N0sex2$result[,"cv"]
    }

    aggr0 <- visit$aggr.output[,c("total","cv")]
    names(aggr0) <- c("fj","cv.fj")
    aggr0$bio.staerri <- B0$aggr.output[,"total"]
    aggr0$cv.bio.staerri <- B0$aggr.output[,"cv"]
    aggr0$fj.minni <- N0$aggr.output[,"total"]
    aggr0$cv.fj.minni <- N0$aggr.output[,"cv"]

    if(sex) {
      aggr0$fj.haenga <- N0sex1$aggr.output[,"total"]
      aggr0$cv.fj.haenga <- N0sex1$aggr.output[,"cv"]
      aggr0$fj.hrygna <- N0sex2$aggr.output[,"total"]
      aggr0$cv.fj.hrygna <- N0sex2$aggr.output[,"cv"]
    }

    aggr0$svaedi <- dimnames(aggr0)[[1]]
    aggr0$svaedisnr <- 1:nrow(aggr0)
    dimnames(aggr0)[[1]] <- 1:nrow(aggr0)
    base0$lengd <- lenClass[i]
    aggr0$lengd <- lenClass[i]
    base0$ar <- yr
    aggr0$ar <- yr

    if(i == 1 ) {
      base <- base0
      aggr <- aggr0
    } else {
      base <- rbind(base,base0)
      aggr <- rbind(aggr,aggr0)
    }
  } # end of lenClass loop

  aggr$bio <- aggr$fj*lwcoeff[1]*aggr$lengd^lwcoeff[2]/1e3
  base$bio <- base$fj*lwcoeff[1]*base$lengd^lwcoeff[2]/1e3

  outdata <- list(base=base,aggr=aggr)
  return(outdata)
}
