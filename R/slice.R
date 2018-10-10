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
##' Whether the sliced data is further segmented by year depends on
##' whether the \code{cslice} object is segmented (i.e., called with
##' \code{split=TRUE}).  The \code{subslice} function also accepts a
##' \code{split} argument to override the default.
##'
##' If provided with a matrix or array of values instead of a vector,
##' these functions operate on the first dimension.
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
##' @param split Whether to segment subsliced data by year.  Defaults
##' to the \code{split} parameter of \code{how}.
##' 
##' @name slice

## Dummy object to make the documentation work

NULL


########################
## helper functions for dealing with array data:

ndim <- function(x){ length(dim(x)) }

## Metafunction that returns a function to index first dimension of x
indexfn <- function(x){
    if(ndim(x) < 2) {
        ## special case for 1-dim x since dim() of a vector is NULL
        function(i, ...){x[i]}
    } else if (ndim(x) == 2) {
        ## special case for 2 dimensions since `[` is faster than do.call()
        function(i, ...){x[i,]}
    } else {
        ## index first dimension of x, regardless of how many there are
        function(i, ...){do.call("[",list(x,i,lapply(dim(x)[-1],seq)))}
    }    
}

ixfn <- function(i, x){ indexfn(x)(i) }



########################
##' @rdname slice
##' @usage slice(x, how, outer=FALSE)
##'
##' @return a list of vectors, one for each climatological window.
##' 
##' @export

slice <- function(x, how, outer=FALSE){
    if(outer){
        index <- how$outer
    } else {
        index <- how$inner
    }
    rapply(index, indexfn(x), how="replace")
}



########################
## More helper functions

## Depth of a nested list
depth <- function(obj, tally=0){
    if(!is.list(obj)){
        return(tally)
    } else {
        return(max(unlist(lapply(obj, depth, tally=tally+1))))    
    }
}

## Converts a list of arrays into a single array
desegmentize <- function(x){
    if(ndim(x[[1]]) < 2) {
        do.call(c, x)
    } else {
        do.call(rbind, x)
    }
}




## Splits an array into a list of arrays at bdys
segmentize <- function(x, bdy){
    indf <- indexfn(x)
    lower <- c(0, bdy)+1
    upper <- c(bdy, ifelse(ndim(x) > 1, dim(x)[1], length(x)))
    mapply(function(a,b){indf(a:b)}, lower, upper, SIMPLIFY=FALSE)
}


########################
##' @rdname slice
##' @usage subslice(s, how, split=how$params$split)
##'
##' @return a list of vectors, one for each climatological window.
##' 
##' @export

## Subslice gets a "split" argument to override the default because
## its most common use case is being called as part of unslice.  Since
## unslice needs the inner slice to be unsegmented, it would be
## pointless to desegementize the data in order to subslice, resegment
## it to return, and then immediately desegmentize it again in the
## calling function.

subslice <- function(s, how, split=how$params$split){
    iind <- lapply(how$inner, unlist)
    oind <- lapply(how$outer, unlist)
    inout <- mapply(match, iind, oind, SIMPLIFY=FALSE)

    ## desegmentize if needed
    if(depth(s) > 1) { s <- lapply(s, desegmentize) }

    ## subslicing
    s <- mapply(ixfn, inout, s, SIMPLIFY=FALSE)

    ## resegmentize if needed
    if(split){
        ibdy <- lapply(inout, function(x){which(diff(x) > 1)})
        s <- mapply(segmentize, s, ibdy)
    }    
    return(s)
}

########################
## Yet more helper functions

len <- function(x){
    if(ndim(x) < 2){
        length(x)
    } else {
        dim(x)[1]
    }
}


########################
##' @rdname slice
##' @usage unslice(s, how)
##' 
##' @return a single timeseries vector.
##' 
##' @export

unslice <- function(s, how, template=NULL, odim=dim(template), onames=dimnames(template)){
    
    ## if s in an outer slice, convert it to an inner slice
    if(identical(rapply(s, len, how="replace"),
                 rapply(how$outer, len, how="replace"))){
        s <- subslice(s, how, split=FALSE)
    }

    ## desegmentize if needed
    if (depth(how$inner) > 1){
        iind <- lapply(how$inner, desegmentize)
    } else {
        iind <- how$inner
    }
    
    if(!identical(sapply(s, len), sapply(iind, len))){
        stop("Sliced data lengths don't match cslice object")
    }
    
    ## easy case - vectors
    if(length(odim) < 2){
        result <- c()
        for(i in 1:length(s)){
            result[iind[[i]]] <- s[[i]]
        }
    } else {
        ## matrices
        result <- array(dim=odim, data=NA)
        if (length(odim) == 2){
            for(i in 1:length(s)){
                result[iind[[i]],] <- s[[i]]
            }
            if(!is.null(onames)){
                dimnames(result) <- onames
            }
        } else {
            ## general arrays are hard and I'm tired
            stop("unslicing arrays with dimension > 2 not supported yet.")
        }
    }
    return(result)
}

