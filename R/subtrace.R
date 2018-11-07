##' Substitute precip values below the trace threshold.
##' 
##' The \code{subtrace} function ensures that the lower tail of the
##' precipitation distribution is smooth (unclustered) and has a
##' reasonable range by replacing values below the trace threshold
##' with values that decay smoothly in density in log space.
##' 
##' Historically, observed precipitation in amounts below 1/100 of an
##' inch (0.254 mm) has generally been unmeasurable, and such amounts
##' are recorded as "trace".  Modeled precipitation can have values
##' recorded down to machine precision.  Bias correction transfer
##' functions can be distorted by clustering of modeled values that
##' are far below the trace threshold, but these values have no
##' importance.  This problem is exacerbated when the lower tail of
##' the distribtion covers a much wider range than the bulk and the
##' upper tail.
##'
##' In univariate bias correction, this problem can be addressed by
##' setting values below the trace threshold to zero (done implicitly
##' by use of \code{dedrizzle} to eliminate excess drizzle), but for
##' multivariate bias correction, it is desirable for all
##' precipitation values to be defined and non-degenerate.  This
##' function imposes those conditions on its input data.
##'
##' The trace threshold defaults to the AMS definition of 1/100 inch
##' (0.254 mm).  This value can be scaled by the grid spacing for
##' area-averaged values.  The default gridscale is 4 km, which is
##' used by a number of gridded observational products.  If a
##' different gridscale is provided, the trace threshold is divided by
##' gridscale^2 / 4^2, which results in a new threshold represening
##' the smallest possible area-average precip value (trace in one
##' location, zero everywhere else) that could result if the gridded
##' obs dataset were aggregated up to coarser resolution.
##' 
##' \code{subtrace} replaces values below the trace threshold with
##' values that are evenly-spaced from the threshold down to the lower
##' limit in log space.  The lower limit is chosen so that it has the
##' same distance from the trace threshold as the upper limit in log
##' space.
##' 
##' @param x A vector of values to substitute subtrace values on.
##' Units presumed to be mm by default values of other arguments.
##'
##' @param traceval The threshold of observable values.  Defaults to
##' 0.254 mm.
##'
##' @param gridscale The grid scale for area average values.  Traceval
##' is scaled by the ratio of gridcell areas.  Defaults to 4 km.
##'
##' @param upper The range from traceval to upper defines the range
##' that substituted values will be distributed over.  Defaults to
##' 1000 mm.
##' 
##' @param att Logical: whether to set attribute "traceval" on the
##' returned vector to record the (scaled) trace threshold.  (Defaults
##' to TRUE.)
##'
##' @return A vector with values below the trace threshold replaced.
##'
##' @examples
##' set.seed(222)
##' x <- rgamma(200, shape=0.2, scale=10)
##' theta <- 0.254
##' print(sum(x < theta))
##' y <- subtrace(x, upper=max(x))
##' 
##' plot(density(log10(x)))
##' lines(density(log10(y)), col="blue")
##' 
##' abline(col="gray", v=log10(c(1000, theta, min(y))))
##' 
##' rug(log10(x))
##' rug(log10(y[y < theta]), col="blue")
##' 
##' @export


subtrace <- function(x, traceval=0.254, gridscale=4, upper=1000, att=TRUE){

    ## scale trace threshold by grid scale
    traceval <- traceval * 16/gridscale^2

    ## find sub-trace values to replace
    subind <- x < traceval
    N <- sum(subind)
    
    ## find lower bound for replacements
    delta <- log(upper) - log(traceval)
    lower <- exp(log(traceval) - delta)

    ## replacement values evenly distributed in log space
    newvals <- exp(seq(from=log(lower), to=log(traceval), length.out=N))    

    y <- x
    y[subind] <- sample(newvals, N)

    if(att){
        y@traceval <- traceval
    }
    
    return(y)
}
