##' Construct a mapping between the distributions of two sets of
##' values.
##'
##' Distribution mapping adjusts the individual values of a dataset, x,
##' such that its statistical distribution matches that of a second
##' dataset,y. This is accomplished by converting data values to
##' probabilities using the CDF of the first distribution, and then
##' from probabilities back to data values using the CDF of the second
##' distribution.
##'
##' The \code{distmap} function constructs this mapping and a transfer
##' function that implements it non-parametrically, by integrating
##' kernel density estimation of the PDFs using the trapezoid rule.  It
##' uses monotone Hermite splines for the construction to guarantee
##' monotonicity of the transfer function.
##' 
##' @param x A vector of values whose distribution is to be mapped.
##'
##' @param y A vector of values whose distribution is the target of
##' the mapping.
##'
##' @param densfun A kernel density estimation function, such as
##' \code{stats::density} or \code{KernSmooth::bkde}.  This function
##' must return an object that can be passed as input to
##' \code{\link{pdf2cdf}}.
##' 
##' @param pgrid A vector of probability values used for the mapping
##' between CDFs.  It need not be regular but must be strictly
##' increasing.
##'
##' @param na.rm Logical; if TRUE (default), remove NA values before
##' constructing distribution mapping.
##' 
##' @param ... Additional arguments to densfun.
##'
##' @return A list of class \code{distmap}, with the following elements:
##' 
##' x,y: The input x and y data.
##'
##' pgrid: The vector of probabilities used to construct the mapping.
##'
##' xpdf,ypdf: The estimated PDFs of x and y.
##'
##' xq,yq: The quantiles of x and y corresponding to the probabilities
##' in pgrid.
##'
##' transfer: A function that will transform x to have the same
##' distribution as y.
##' 
##'
##' @examples
##' library(nor1mix)
##' set.seed(222)
##' 
##' x <- rnorMix(1e6, norMix(mu=c(-3,2),sigma=c(1,2),w=c(2,1)))
##' y <- rnorMix(1e6, norMix(mu=c(-2,2)))
##' dmap <- distmap(x,y)
##' z <- predict(dmap)
##' 
##' title = "PDFs; z = x mapped to match y"
##' plot(NA, type="n", xlim=c(-10,10), ylim=c(0,0.33), main=title)
##' lines(density(x), col="blue",   lwd=2)
##' lines(density(y), col="black", lwd=3)
##' lines(density(z), col="red",  lwd=2, lty=2)
##' legend("topright", c("x","y","z"),  col=c("blue","black","red"), lty=c(1,1,2))
##'
##' dev.new()
##' 
##' pp <- pnorm(seq(-2,2))
##' yr <- c(-5,5)
##' 
##' par(mfrow=c(3,1))
##' plot(y[1:100],type="b",ylim=yr,xlab="")
##' abline(h=quantile(y,pp), lty=c(4,3,2,3,4))
##' 
##' plot(x[1:100],type="b",col="red", ylim=yr, xlab="")
##' abline(h=quantile(x,pp), col="red", lty=c(4,3,2,3,4))
##' 
##' plot(z[1:100],type="b",col="blue", ylim=yr, xlab="")
##' abline(h=quantile(z,pp), col="blue", lty=c(4,3,2,3,4))
##'
##' 
##' @importFrom KernSmooth bkde
##' @importFrom stats splinefun ppoints
##'
##' @export



distmap <- function(x, y, densfun=bkde, pgrid=ppoints(1000), na.rm=TRUE, ...){

    if(na.rm){
        x <- x[!(is.na(x))]
        y <- y[!(is.na(y))]
    }
    
    stopifnot(is.atomic(x), all(is.finite(x)),
              is.atomic(y), all(is.finite(y)),            
              is.atomic(pgrid), all(is.finite(pgrid)),
              all(pgrid >= 0), all(pgrid <= 1),
              all(diff(pgrid) > 0)
              )
   
    xpdf <- densfun(x, ...)
    ypdf <- densfun(y, ...)
    
    xcdf <- pdf2cdf(xpdf)
    ycdf <- pdf2cdf(ypdf)
    
    xq <- splinefun(xcdf$y, xcdf$x, method="monoH.FC")(pgrid)
    yq <- splinefun(ycdf$y, ycdf$x, method="monoH.FC")(pgrid)

    transfer <- splinefun(xq, yq, method="monoH.FC")

    result <- list(x=x, y=y, pgrid=pgrid, xpdf=xpdf, ypdf=ypdf,
                   xq=xq, yq=yq, transfer=transfer)
                   
    class(result) <- "distmap"
    return(result)
}

### Copyright 2015 Univ. Corp for Atmos. Research
### Author: Seth McGinnis, mcginnis@ucar.edu
