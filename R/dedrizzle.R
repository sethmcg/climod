##' Removes excess drizzle from model data.
##' 
##' Given three vectors of precipitation (current and future model
##' data, plus observations), dedrizzle finds a threshold value that
##' will equalize the wet/dry ratio in the model data to match the
##' wet/dry ratio in the observations, and sets all model values below
##' that threshold to zero.
##'
##' Dedrizzle floors all its inputs at zero before calculation.  It
##' ignores NA values when calculating the wet/dry equalization
##' threshold.
##'
##' @param x A list of vectors named "obs", "cur", and "fut".
##'
##' @return A list with obs, cur, and fut with drizzle removed.
##'
##' @examples
##' obs <- c(rep(0,10), rep(1,10), rep(NA,3))
##' cur <- c(seq(10)/10, NA)
##' fut <- seq(15)/10
##' namelist(obs,cur,fut)
##' dedrizzle(namelist(obs, cur, fut))
##' 
##' @export



dedrizzle <- function(x){
    
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

    return(x)
}

