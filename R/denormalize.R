utils::globalVariables(c("norm","from","to","scale","power","mean","sd"))
## The overloaded use of @ involves non-standard evaluation that makes
## devtools::check() (or rather, codetools::checkUsagePackage()) think
## that the attribute names are undefined global variable.  Declaring
## them to be globalVariables suppresses the NOTE that would otherwise
## result.  (It actually only complains about "from","to", and
## "power", because the others all correpond to functions that come
## from base:: or have been imported elsewhere, but I list them all
## here to make it clearer what's going on.)

##' Undo the normalization of a vector of values
##'
##' Reverses the transformation applied to a vector of values by the
##' function \code{normalize}, using the type and parameters of the
##' transformation as recorded in the vector's attributes.  See
##' \code{\link{normalize}} for details on the transformations.
##'
##' The ability to undo a normalization is useful in bias correction
##' of climate model output.  A bias correction will adjust the values
##' of the climate model output to match the statistical distribution
##' of observational data for a historical period, and then apply the
##' corresponding adjustment to model output for a future period.
##' Often, the performance of the bias correction will be improved if
##' the various inputs are first normalized in some way.  This
##' requires the normalization to be reversed after the bias
##' correction has been applied.
##'
##' After bias-correcting current period model output, the
##' normalization should be reversed in a way that matches the
##' obervations, not the input model data.  The \code{match} argument
##' is used for this purpose.  Similarly, after bias-correcting future
##' period model output, the normalization needs to be reversed so
##' that the result matches the observations, adjusted by the
##' difference between the original current and future model outputs.
##' The \code{adjust} argument is used for this purpose.
##'
##' The "power" transformation raises the data to an arbitrary power.
##' When undoing this transformation, the data is first floored at
##' zero to avoid problems with negative inputs.
##' 
##' @param x A normalized vector
##' 
##' @param match A normalized vector whose parameters should be used to
##' reverse the normalization.
##'
##' @param adjust A normalized vector; the parameters of \code{match}
##' are adjusted by the difference or ratio (as appropriate) of the
##' paramters of \code{adjust} and \code{x} when reversing the
##' normalization.
##'
##' @examples
##'
##' obs <- rgamma(10000, shape=5, scale=3)
##' cur <- rgamma(10000, shape=6, scale=2)
##' fut <- rgamma(10000, shape=3, scale=4)
##'
##' data <- namelist(obs, cur, fut)
##' ndata <- lapply(data, normalize, norm="power")
##'
##' dmap <- distmap(ndata$cur, ndata$obs)
##' ndata$bcc <- predict(dmap)
##' ndata$bcf <- predict(dmap,ndata$fut)
##'
##' N <- length(ndata)
##' mplot(lapply(ndata,density), type="l")
##' legend("topleft",names(ndata),lwd=1,lty=seq(N),col=seq(N))
##'
##' denorm <- lapply(ndata[1:3], denormalize)
##' denorm$bcc <- denormalize(ndata$bcc, match=ndata$obs)
##' denorm$bcf <- denormalize(ndata$bcf, match=ndata$obs, adjust=ndata$cur)
##'
##' N <- length(denorm)
##' mplot(lapply(denorm,density), type="l")
##' legend("topright",names(denorm),lwd=1,lty=seq(N),col=seq(N)) 
##'
##' @seealso \code{\link{normalize}}
##'
##' @importFrom scales rescale
##' 
##' @export


denormalize <- function(x, match=NULL, adjust=NULL){

    if(is.null(match) & ! is.null(adjust)){
        stop("argument 'match' required if 'adjust' is provided")
    }

    if(!is.null(match)){
        norm <- match@norm
    } else {
        norm <- x@norm
    }

    if(!norm %in% c("range", "zscore", "power")){
        stop(paste("unknown normalization",norm))
    }
     
    if(norm == "range"){
        if(!is.null(match)){
            out.range <- match@from
            if(!is.null(adjust)){
                out.range <- out.range + x@from - adjust@from
            }
        } else {
            out.range <- x@from
        }
        
        result <- rescale(x, from=x@to, to=out.range)
        result@range <- NULL
    }

    if(norm == "zscore"){
        if(!is.null(match)){
            out.mean <- match@mean
            out.sd   <- match@sd
            if(!is.null(adjust)){
                out.mean <- out.mean + x@mean - adjust@mean
                out.sd   <- out.sd   * x@sd   / adjust@sd
            }
        } else {
            out.mean <- x@mean
            out.sd   <- x@sd
        }
              
        result <- x * out.sd + out.mean
        result@mean <- NULL
        result@sd   <- NULL
    }

    if(norm=="power"){
        if(!is.null(match)){
            out.scale <- match@scale
            out.power <- match@power
            if(!is.null(adjust)){
                out.scale <- out.scale * x@scale / adjust@scale
                out.power <- out.power * x@power / adjust@power
                if(x@power != adjust@power){
                    warning("Arguments x and adjust were normalized using norm=\"power\" with different exponents.")
                }
            }
        } else {
            out.scale <- x@scale
            out.power <- x@power
        }

        ## floor before exponentiation to avoid NaN (or worse, if power is 1/even...)
        x <- pmax(x,0)
        
        if (out.power == 0){
          result <- exp(x)
        } else {
          result <- x^(1/out.power)
        }
        result <- result * out.scale
	result@power <- NULL
	result@scale <- NULL
    }
    
    result@norm <- NULL
    return(result)  
}


### Copyright 2015 Univ. Corp for Atmos. Research
### Author: Seth McGinnis, mcginnis@ucar.edu
