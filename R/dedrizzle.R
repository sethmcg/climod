##' Removes excess drizzle from model data.
##' 
##' Climate models typically generate too much drizzle.  Given modeled
##' and observational precipitation, \code{dedrizzle} finds a
##' threshold value that will cause the wet/dry ratio in the model
##' data to match the wet/dry ratio in the observations, then sets all
##' model values below that threshold to zero.
##'
##' \code{dedrizzle} takes two vectors of model data, one for the
##' current period and one for the future period.  It finds the
##' threshold based on the current data, and floors data below the
##' threshold for both current and future data.  If the observational
##' data has more wet (non-zero) values than the model data, the
##' threshold will be zero.
##' 
##' Dedrizzle floors all its inputs at zero before calculation.  It
##' ignores NA values when calculating the wet/dry equalization
##' threshold.
##'
##' @param x A list of vectors named "obs", "cur", and "fut".  The
##' vectors do not need to be the same length.
##'
##' @param att Logical: whether to set attributes "pwet" and
##' "threshold" on the returned list to record the wet/dry ratio and
##' threshold value, respectively.  (Defaults to FALSE.)
##'
##' @return A list containing obs, cur, and fut with drizzle removed.
##'
##' @examples
##' obs <- c(rep(0,10), rep(1,10), rep(NA,3))
##' cur <- c(seq(10)/10, NA)
##' fut <- seq(15)/10
##' namelist(obs,cur,fut)
##' dedrizzle(namelist(obs, cur, fut))
##' 
##' @export



dedrizzle <- function(x, att=FALSE){
    
    ## floor x at zero
    x <- lapply(x, pmax, 0)    

    ## discard NA values
    y <- lapply(x, function(z){z[!is.na(z)]})
    
    ## calculate observed wet/dry fraction
    pwet <- sum(y$obs > 0) / length(y$obs)
    
    ## find threshold equalizing model wet/dry with obs
    threshold <- sort(y$cur)[(1-pwet)*length(y$cur)]
    
    ## set values at or below threshold to zero
    x$cur[x$cur <= threshold] <- 0
    x$fut[x$fut <= threshold] <- 0

    if(att){
        x@pwet <- pwet        
        x@threshold <- threshold
    }
    
    return(x)
}

## Unit test:
## x <- c(rep(0,10), seq(10))  ## length 20, 50% dry
## y <- c(rep(0,10), seq(20))  ## length 30, 33% dry
## z <- c(rep(0,15), seq(45)   ## length 60, 25% dry

## obs/cur/fut  theta   future
## x / y / z => 5       1/3 dry
## x / z / y => 15      5/6 dry
## y / z / x => 5       3/4 dry
## y / x / z => 0       1/4 dry
## z / x / y => 0       1/3 dry
## z / y / x => 0       1/2 dry

## use att=TRUE in the unit tests.

## Unit test: x <- c(x,rep(5,NA))

## unit test: y[1:5] <- -2

## Unit test: shuffle all the vectors first

## Unit test: what happens when fut is NULL?
