##' Trims extreme values from a vector that exceed Chauvenet's criterion.
##'
##' To prevent overfitting of the distributional tails when
##' constructing a \code{\link{distmap}} object, it is useful to trim
##' very extreme values from the input datasets.  This function trims
##' values whose probability is less than 1/N (greater than 1-1/N),
##' where N is the number of data points + 1/2, under the assumption that
##' the data is gaussian.
##'
##' The probability bounds are calculated using \code{qnorm} with the
##' sample mean and standard deviation.
##'
##' @param x A vector of values to trim.
##'
##' @param upper Logical; trim values in the upper tail.  Defaults to TRUE.
##'
##' @param lower Logical; trim values in the lower tail.  Defaults to TRUE.
##'
##' @return A trimmed vector of values.
##'
##' @examples
##' x <- c(9, 10, 10, 10, 11, 50)
##' chauvenet(x)
##' chauvenet(c(x, -20), upper=FALSE)
##' @export

chauvenet <- function(x, upper=TRUE, lower=TRUE){

  n <- length(x) + 1/2
  mu  <- mean(x)
  sigma <- sd(x)

  if(lower){ x <- x[!(x < qnorm(1/n, mu, sigma))] }
  if(upper){ x <- x[!(x > qnorm(1-1/n, mu, sigma))] }

  return(x)
}
