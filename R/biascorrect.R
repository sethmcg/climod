##' Bias correction of climate model output using KDDM.
##'
##' This is the description of the function.  It's the first full
##' paragraph, and should briefly describe what the function does.
##'
##' Third and subsequent paragraphs are details: a long section shown
##' after the argument that goes into detail about how the function
##' works.
##'
##' @param bcdata The data to be bias-corrected.  A list with elements
##' named "obs", "cur", and "fut".
##'
##' @param norm The type of normalization to use.  See
##' \code{\link{normalize}} for more details.
##'
##' @param dmap Logical; if TRUE, returns the \code{\link{distmap}}
##' object generated during bias correction as an element named "dmap"
##' in the returned list.  N.B.: This option defaults to FALSE because
##' the distmap object is significantly larger than the inputs and
##' outputs of this function combined.
##' 
##' @param ... Arguments to pass to \code{\link{distmap}}.
##' 
##' @return The return value of the function.
##'
##' @examples
##' x <- "example code"
##' y <- "more example code"
##' paste(x, y)
##' 
##' @export

biascorrect <- function(bcdata, norm="zscore", dmap=FALSE, ...){

    ## normalize the three data components
    nbcd <- lapply(bcdata, normalize, norm)

    ## construct distribution mapping
    mapping <- distmap(nbcd$cur, nbcd$obs, ...)

    ### don't subslice yet...

    ## drop obs data
    fixme <- nbcd[-(names(nbcd)=="obs")]

    ## apply KDDM transfer function
    fixed <- lapply(fixme, function(x){predict(mapping, x)})


    ### hey, normalization atts carry over now!
    
    
    ## denormalize bias-corrected data
    result <- list()
    result$obs <- bcdata$obs
    result$cur <- denormalize(fixed$cur, match=nbcd$obs)
    result$fut <- denormalize(fixed$fut, match=nbcd$obs, adjust=nbcd$cur)

    if(dmap){
        result$distmap <- mapping
    }
    
    copyatts(bcdata, result)
    return(result)
}

