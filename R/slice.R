##' Climatological slicing of timeseries
##'
##' \code{slice}, \code{subslice}, and \code{unslice} are used to
##' separate a timeseries into climatological windows, then
##' reconstruct the timeseries from the windowed data.
##'
##' \code{slice} uses a \code{\link{cslice}} object to separate a
##' timeseries into climatological windows.  \code{unslice}
##' reconstructs the timeseries from the (possibly modified) sliced
##' data.  \code{subslice} takes data that has been sliced into outer
##' windows and extracts the inner windows from each slice.
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
    if(how$param$split){
      index <- unlist(index)
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

