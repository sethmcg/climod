##' Reversibly normalize a vector of values
##'
##' Normalizes or standardizes a vector of values, recording the
##' details of the transformation in attributes of the result so that
##' the transformation can be reversed later on.
##'
##' Three types of normalization are supported, chosen by the value of
##' the parameter \code{norm}:
##' 
##' \code{"zscore"}: subtract the mean and divide by the standard
##' deviation.  The result will have mean zero and unit standard
##' deviation.  Also known as normalizing residuals, calculating the
##' z-score, centering and scaling, or normalizing the 2nd central
##' moment.  If the sample mean and sample variance are used (as in
##' the default arguments), this is technically 'studentizing' the
##' data.
##'
##' \code{"range"}: subtract the minimum and divide by the range
##' (max-min).  The default is to transform the data to have min 0 and
##' max 1.
##'
##' \code{"power"}: divide by the scale (variance/mean) and raise to a
##' power.  This transformation will stabilize the variance for
##' highly-skewed data.
##'
##' Each normalization by default uses parameters based on the sample
##' statistics of the input, but these parameters can be overridden.
##' 
##' 
##' @param x A vector of values to be transformed.
##'
##' @param norm The type of normalization to apply.  See Details for
##' more information.  Type names can be abbreviated.
##'
##' @param mean For zscore normalization, the mean.  Defaults to the
##' sample mean.
##'
##' @param sd For zscore normalization, the standard deviation.
##' Defaults to the sample standard deviation.
##'
##' @param from For range normalization, the input range.  Defaults to
##' range(x).
##'
##' @param to For range normalization, the output range.  Defaults to
##' c(0,1).
##'
##' @param scale For power normalization, the scaling factor to divide
##' by before the power transform.  Defaults to the sample variance
##' divided by the sample mean, which is the scale parameter of the
##' gamma and related distributions.
##'
##' @param power For power normalization, the exponent of the power
##' transform. Defaults to 0.25, which is empirically known to be
##' well-suited for working with precipitation.
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
##' # range normalization
##' x <- sample(2:9, 5, replace=TRUE)
##' y <- normalize(x, norm="range", from=c(1,10))
##' range(x)
##' range(y)
##' str(y)
##' 
##' # power normalization
##' x <- rgamma(10000, shape=3, rate=4)
##' y <- normalize(x, "power")
##' mplot(lapply(list(x,y),density), type="l")
##' str(y)
##' 
##'
##' @seealso \code{\link{denormalize}}
##'
##' @importFrom scales rescale
##' @importFrom stats var sd
##' 
##' @export


## TODO: consider adding car::bcPower and car::yjPower (box-cox and
## yeo-johnson) options for power normalization


normalize <- function(x,
                      norm="zscore",
                      mean=base::mean(x, na.rm=TRUE),
                      sd=stats::sd(x, na.rm=TRUE),
                      from=range(x, na.rm=TRUE),
                      to=c(0,1),
                      scale=stats::var(x, na.rm=TRUE)/base::mean(x, na.rm=TRUE),
                      L2=1,
                      power=0.25){

    norm.names <- c("range", "zscore", "power", "boxcox", "log")
    n <- norm.names[pmatch(norm, norm.names)]
    if(is.na(n)){
      stop(paste("unknown normalization",norm))
    } else {
      norm <- n
    }

    result <- NA
  
    if(norm == "range"){
        result <- rescale(x, to=to, from=from)
        result@from <- from
        result@to   <- to
       }

    if(norm == "zscore"){
        result <- (x - mean)/sd
        result@mean <- mean
        result@sd   <- sd
    }

    if(norm == "log"){
        result <- log(x)
    }
    
    if(norm=="power"){
        result <- x / scale
        if (power == 0){
          result <- log(result)
        } else {
          result <- result^power
        }
        result@scale <- scale
        result@power <- power
    }

    if(norm=="boxcox"){
        result <- x + L2
        if (power == 0){
          result <- log(result)
        } else {
          result <- result^power
        }
        result@L2    <- L2
        result@power <- power
    }

    result@norm <- norm
    return(result)
}



### Copyright 2015 Univ. Corp for Atmos. Research
### Author: Seth McGinnis, mcginnis@ucar.edu
