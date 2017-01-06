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
##' When reversing the normalization, correction factors can be
##' supplied to adjust the normalization parameters to compensate for
##' bias.  The \code{shift} and \code{scale} arguments apply additive
##' and multiplicative adjustments, respectively, to the data.  The
##' \code{pscale} argument applies a multiplicative adjustment to the
##' exponent used in power normalization.
##'
##' The "power" transformation raises the data to an arbitrary power.
##' When undoing this transformation, the data is first floored at
##' zero to avoid problems with negative inputs.
##' 
##' @param x A normalized vector
##' 
##' @param shift An adjustment factor added to the denormalized data.
##' Only used for "zscore" and "range" normalizations.  Adjusts the
##' mean of zscore-normalized data.  For range-normalized data, a
##' 2-element vector can be supplied to adjust both ends of the range.
##'
##' @param scale A multiplicative adjustment factor applied to the
##' denormalized data.  Only used for the "zscore" normalization.
##' Adjusts the variance of zscore-normalized data.
##'
##' @param pscale A multiplicative adjustment factor applied to the
##' exponent when denormalizing power-transformed data.
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
##' ndata$bcc <- predict(dmap, ndata$cur)
##' ndata$bcf <- predict(dmap, ndata$fut)
##' 
##' par(mfrow=c(2,1))
##' 
##' N <- length(ndata)
##' mplot(lapply(ndata,density), type="l")
##' legend("topleft",names(ndata),lwd=1,lty=seq(N),col=seq(N))
##' 
##' denorm <- lapply(ndata[1:3], denormalize)
##' adjust <- ndata$obs@scale / ndata$cur@scale
##' denorm$bcc <- denormalize(ndata$bcc, SCALE=adjust)
##' denorm$bcf <- denormalize(ndata$bcf, SCALE=adjust)
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


denormalize <- function(x, shift=0, scale=1, pscale=1){

    norm <- x@norm
    
    if(!norm %in% c("range", "zscore", "power")){
        stop(paste("unknown normalization",norm))
    }
     
    if(norm == "range"){

        out.range <- x@from + shift
        
        result <- rescale(x, from=x@to, to=out.range)
        result@range <- NULL
    }

    if(norm == "zscore"){

        out.mean <- x@mean + shift
        out.sd   <- x@sd   * scale
           
        result <- x * out.sd + out.mean
        result@mean <- NULL
        result@sd   <- NULL
    }

    if(norm=="power"){

        out.shift <- x@shift + shift
        out.power <- x@power * pscale

        ## floor before exponentiation to avoid NaN (or worse, if 1/power is even...)
        x <- pmax(x,0)
        
        if (out.power == 0){
          result <- exp(x)
        } else {
          result <- (out.power * x + 1)^(1/out.power)
        }
        result <- result - out.shift
	result@power <- NULL
	result@shift <- NULL
    }
        
    result@norm <- NULL
    return(result)  
}


### Copyright 2015 Univ. Corp for Atmos. Research
### Author: Seth McGinnis, mcginnis@ucar.edu
