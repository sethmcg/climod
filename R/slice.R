##' Climatological slicing of timeseries
##'
##' \code{slice}, \code{subslice}, and \code{unslice} are used to
##' separate a timeseries into climatological windows, then
##' reconstruct the timeseries from the windowed data.
##'
##' \code{slice} uses a \code{\link{cslice}} object to separate a
##' timeseries into climatological windows.  \code{unslice}
##' reconstructs the timeseries from the (possibly modified) windowed
##' data.  \code{subslice} 
##'
##' @param x A vector of values to be sliced into windows.
##'
##' @param s Sliced data -- a list of vectors, each corresponding to
##' a different climatological window.
##'
##' @param how A \code{cslice} object defining the windows to use for
##' slicing.
##'
##' @param outer Logical: if TRUE, uses the outer windows from
##' \code{how}, which can overlap.  If FALSE (the default), uses inner
##' windows, which cover the data completely without overlapping.
##' 
##' @name slice

## Dummy object to make the documentation work

NULL


##' @rdname slice
##' @usage slice(x, how, outer=FALSE)
##' @export

slice <- function(x, how, outer=FALSE){
    if(outer){
        index <- how$outer
    } else {
        index <- how$inner
    }
    lapply(index, function(i,y){return(y[i])}, x)
}


##' @rdname slice
##' @usage subslice(s, how)
##' @export

subslice <- function(s, how){
    inout <- mapply(match, how$inner, how$outer, SIMPLIFY=FALSE)
    mapply('[', s, inout, SIMPLIFY=FALSE)
}


##' @rdname slice
##' @usage unslice(s, how)
##' @export

unslice <- function(s, how){
    if(all(sapply(s, length) == sapply(how$outer, length))){
        s <- subslice(s, how)
    }
    if(!all(sapply(s, length) == sapply(how$inner, length))){
        stop("Sliced data lengths don't match cslice object")
    }
    result <- c()
    for(i in 1:length(s)){
        result[how$inner[[i]]] <- s[[i]]
    }
    return(result)
}


plot.cslice <- function(cs, inner.args=NULL, ...){
#    itime = slice(cs$time, cs, outer=FALSE)
#    otime = slice(cs$time, cs, outer=TRUE)
#    mplot(otime, cs$outer, pch=pcho, ...)
#    mapply(points, itime, cs$inner, inner.args)
}


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
