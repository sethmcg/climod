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
##' @importFrom stats predict
##' 
##' @export

biascorrect <- function(bcdata, norm="zscore", dmap=FALSE, ...){

    ## normalize the three data components

    if(norm=="boxcox"){
        gamma <- lapply(lapply(bcdata, unlist), maxent)
        nbcd <- list()
        nbcd$obs <- rapply(bcdata$obs, normalize, how="replace", norm=norm, gamma=gamma$obs)    
        nbcd$cur <- rapply(bcdata$cur, normalize, how="replace", norm=norm, gamma=gamma$cur)
        nbcd$fut <- rapply(bcdata$fut, normalize, how="replace", norm=norm, gamma=gamma$cur)
        # yes, using gamma$cur for the future case is intentional
    } else {    
        nbcd <- rapply(bcdata, normalize, how="replace", norm=norm)
    }
    
    ## construct distribution mapping
    mapping <- distmap(unlist(nbcd$cur), unlist(nbcd$obs), ...)

    ## apply KDDM transfer function
    ## note: could skip predict(obs), but need its norm atts for next step
    fixed <- rapply(nbcd, function(x){predict(mapping, x)}, how="replace")

    ## extract & average normalization attributes
    adj <- rapply(fixed, attributes, how="replace")
    adj <- lapply(adj, renest)
    adj <- lapply(adj, function(x){lapply(x, function(y){unlist(y)})})
    adj <- lapply(adj, function(x){x$norm<-NULL; x})
    adj <- rapply(adj, mean, how="replace")
    adj <- renest(adj)

    ## calculate adjustments for denormalization
    shift <- NA
    scale <- NA
    pscale <- NA

    if(norm == "zscore"){
        shift  <- adj$mu$obs - adj$mu$cur
        scale  <- adj$sigma$obs   / adj$sigma$cur
    }

    if(norm == "boxcox"){
        shift <- adj$gamma$obs - adj$gamma$cur
        pscale <- 1
    }

    ## denormalize bias-corrected data
    result <- list()
    result$obs <- bcdata$obs

    result$cur <- rapply(fixed$cur, denormalize, how="replace",
                         shift=shift, scale=scale, pscale=pscale)
    result$fut <- rapply(fixed$fut, denormalize, how="replace",
                         shift=shift, scale=scale, pscale=pscale)
    
    if(dmap){
        result$distmap <- mapping
    }
    
#    copyatts(bcdata, result)
    return(result)
}

