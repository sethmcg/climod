##' Reversibly normalize a vector of values
##'
##' Normalizes or standardizes a vector of values, recording the
##' details of the transformation in attributes of the result so that
##' the transformation can be reversed later on.
##'
##' Six types of normalization are supported, chosen by the value of
##' the parameter \code{norm}:
##' 
##' \code{"zscore"}: subtract the mean (mu) and divide by the standard
##' deviation (sigma).  The result will have mean zero and unit
##' standard deviation.  Also known as normalizing residuals,
##' calculating the z-score, centering and scaling, or normalizing the
##' 2nd central moment.  If the sample mean and sample variance are
##' used (as in the default arguments), this is technically
##' 'studentizing' the data.
##'
##' \code{"boxcox"}: apply the Box-Cox power transformation, which
##' raises the data to a power and scales it appropriately such that
##' it's continuous down to a power of 0, which becomes the log.  This
##' transformation will stabilize the variance for highly-skewed data.
##'
##' \code{"log"}: take the (natural) log of the data.  Equivalent to
##' boxcox with lambda=0 and gamma=0.
##'
##' \code{"scale"}: divide by the mean of the absolute value of the
##' data (mu2).  This normalization is appropriate for data where
##' zero is a natural limit to the range of values.
##'
##' \code{"range"}: scale the data from the provided range
##' [lower,upper] to the range [0,1].
##' 
##' \code{"identity"}: passes the data through unchanged.  This option
##' is sometimes useful for testing and development.
##' 
##' Each normalization by default uses parameters based on the sample
##' statistics of the input, but these parameters can be overridden.
##'
##' Appropriate normalizations for various climate variables used in
##' impacts analyis are as follows:
##'
##' \itemize{
##' \item{Temperatures:}{zscore.}
##'
##' \item{Precipitation:}{log or boxcox.  Although the Box-Cox
##' transformation works well for normalizing observed precipitation,
##' which has a cutoff at some trace threshold, there are problems
##' using it to normalize model data for distribution mapping.}
##'
##' \item{Winds:}{zscore for directional wind, but scale for wind
##' speed.}
##'
##' \item{Humidity:}{scale for specific humidity, range for relative
##' humidity.}
##'
##' \item{solar radiation:}{scale.}
##' }
##' 
##' @param x A vector of values to be transformed.
##'
##' @param norm The type of normalization to apply.  See Details for
##' more information.  Type names can be abbreviated.
##'
##' @param mu For zscore normalization, the mean.  Defaults to the
##' sample mean.
##'
##' @param sigma For zscore normalization, the standard deviation.
##' Defaults to the sample standard deviation.
##'
##' @param gamma An offset to the data for the Box-Cox transform
##' normalization.  Must be greater than minus the minimum value,
##' otherwise NA values will result.  Defaults to zero.
##'
##' @param lambda For the Box-Cox transform normalization, the
##' exponent of the power transform. Defaults to zero.
##'
##' @param mu2 For scale normalization, the mean.  Defaults to the
##' sample mean of the absolute value of x.
##'
##' @param lower For range normalization, the lower bound of the input
##' range.  Defaults to 0.
##'
##' @param upper For range normalization, the upper bound of the input
##' range.  Defaults to 100.
##'
##' @param clip If TRUE (default), negative values will be clipped to
##' zero for scale normalization, and values outside the input range
##' will be clipped to the input range for range normalization.
##' 
##' @return Both \code{normalize} and \code{denormalize} return a
##' vector of values.  \code{normalize} adds attributes to the vector
##' that record the type of normalization and the parameters used;
##' \code{denormalize} removes these attributes from its output.
##'
##' @examples
##'
##' # z-score normalization
##' x <- rnorm(10000, mean=3, sd=2)
##' y <- normalize(x)
##' mplot(lapply(list(x,y),density), type="l")
##' str(y)
##' 
##' # Box-Cox normalization
##' x <- rgamma(10000, shape=3, rate=4)
##' y <- normalize(x, "boxcox")
##' mplot(lapply(list(x,y),density), type="l")
##' str(y)
##' 
##'
##' @seealso \code{\link{denormalize}}
##'
##' @importFrom stats sd
##' 
##' @export


normalize <- function(x,
                      norm="zscore",
                      mu=base::mean(x, na.rm=TRUE),
                      sigma=stats::sd(x, na.rm=TRUE),
                      lambda=0,
                      gamma=0,
                      mu2=base::mean(base::abs(x), na.rm=TRUE),
                      lower=0,
                      upper=100,
                      clip=TRUE){

    norm.names <- c("zscore", "boxcox", "log", "scale", "range", "identity")
    n <- norm.names[pmatch(norm, norm.names)]
    if(is.na(n)){
      stop(paste("unknown normalization",norm))
    } else {
      norm <- n
    }

    result <- NA
  
    if(norm == "zscore"){
        result       <- (x - mu)/sigma
        result@mu    <- mu
        result@sigma <- sigma
    }

    if(norm=="boxcox"){
        x <- x + gamma     
        if (lambda == 0){
            result <- log(x)
        } else {
            result <- (x^lambda - 1) / lambda
        }
        result@lambda <- lambda
        result@gamma  <- gamma
    }

    if(norm=="log"){
        result <- log(x)
    }

    if(norm=="scale"){
        if(clip){ x[x < 0] <- 0 }
        result <- x / mu2
        result@mu2  <- mu2
        result@clip <- clip
    }

    if(norm=="range"){
        if(clip){
            x[x < lower] <- lower
            x[x > upper] <- upper
        }
        result <- (x - lower) / (upper - lower)
        result@lower <- lower
        result@upper <- upper
        result@clip  <- clip
    }
    
    if(norm=="identity"){
      result <- x
    }
    
    result@norm <- norm
    return(result)
}



### Copyright 2015 Univ. Corp for Atmos. Research
### Author: Seth McGinnis, mcginnis@ucar.edu
