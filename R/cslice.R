##' Climatological slicing of cyclic data
##'
##' Creates an object that indexes climatological moving 
##' windows for cyclic data.
##'
##' DETAILS GO HERE
##'
##' Allows for irregularly spaced data and non-standard calendars
##' 
##' The size and number of the inner and outer windows can be
##' specified using any two of the four parameters \code{num},
##' \code{ratio}, \code{inner}, and \code{outer}.  By default, the
##' \code{inner} and \code{outer} are assumed to have units of days,
##' but other units should work as long as they are consistent with
##' the value of \code{year}.
##'
##' @param time A vector of times, encoded as time elapsed since some
##' start date.
##'
##' @param num The number of windows.
##'
##' @param ratio The length of the outer window as a multiple of the
##' length of the inner window.
##'
##' @param inner The length of the inner window.
##'
##' @param outer The length of the outer window.
##'
##' @param year The length of the year.
##'
##' @param names A vector of names for the windows.
##'
##' @return An object of class 'cslice' containing lists of indices
##' associated with moving windows across multiple years and a list of
##' the parameters specifying the object.
##'
##' @seealso \code{\link{slice}}
##' 
##' @export

# add plot, print, str? methods


cslice <- function(time,
                   num = ceiling(year / inner),
                   ratio = outer / inner,
                   inner = year / num,
                   outer = inner * ratio,
                   year  = yearlength(time),
                   names = NULL
                   ){

#    if(outer < inner){ stop("outer must be >= inner")}
    ## check for:
    ## inner > outer
    ## num > year? well, no, could have subdaily.  Or diff units.
    ## all finite, positive
    
    N = length(time)
    
    span <- diff(range(time))
    cycles <- span / year
  
    cs <- list(
        time   = time,
        inner  = vector(mode="list", length=num),
        outer  = vector(mode="list", length=num),
        params = namelist(span, cycles, year, ratio, inner, outer)
        )

    class(cs) <- "cslice"
    
    names(cs$inner) <- names
    names(cs$outer) <- names

    time <- (time - time[1])

    for(i in 1:num){
        cs$inner[[i]] <- which( (time - inner*(i-1)) %% year < inner)
        cs$outer[[i]] <- which( (time - inner*(i-1/2) + outer/2 ) %% year < outer)
    }    
    return(cs)
}

## examples for cslice probably just use plot/print methods
 
# 
# ## examples
# library(ncdf4)
#  
# nc <- nc_open("test/tmax.gcm.cur.nc")
# time <- ncvar_get(nc, "time")
# time@calendar <- ncatt_get(nc, "time", "calendar")$value
# tmax <- ncvar_get(nc,"tmax")[1,1,]
# 
# 
# ## inner slice example: separate data into 12 "monthly" windows
# 
# months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# monthly <- cslice(time, num=12, ratio=1, names=months)
# 
# mondata <- slice(tmax, monthly)
# montime <- slice(time, monthly)
#  
# xr <- c(0,365*5)
# yr <- range(unlist(tmax))
# col <- topo.colors(12); names(col) <- months
# plot(NA, xlim=xr, ylim=yr, xlab="time", ylab="tmax")
# for(i in months){
#     points(montime[[i]],mondata[[i]],col=col[i])
# }
# abline(h=sapply(mondata,mean),col=col)
# 
# 
# ## outer slice example: normalizing daily data based on 30-day moving
# ## window vs normalizing monthly
# 
# mnorm <- lapply(mondata, normalize)
# monanom <- unslice(mnorm, monthly)
# 
# 
# daily <- cslice(time, inner=1, outer=30)
# ddata <- slice(tmax, daily, outer=TRUE)
# 
# ndata <- lapply(ddata, normalize)
# dayanom <- unslice(ndata, daily, outer=TRUE)
# 
# doy <- rep(1:365,round(monthly$param$cycles))[1:length(tmax)]
# plot(doy, monanom, xlim=c(1,30))
# points(doy, dayanom, pch='.', col="red", cex=3)
# 
# ## Enh. Not the greatest example.  May do
# 
