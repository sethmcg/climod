##' Plot multiple data series
##'
##' Given a list of objects, extracts the x and y elements from each
##' object, joins them together in a matrix, and passes them to
##' matplot.
##'
##' This is a convenience function used primarily for plotting
##' multiple PDFs (from, e.g., \code{stats::density} or
##' \code{KernSmooth::bkde}) on the same graph.
##' 
##' Each object must be a list with elements named "x" and "y", and
##' all x and y elements must be the same length.  (The elements are
##' merged into matrices via calls to sapply.)  These restrictions
##' make the function not terribly useful in the general case, but it
##' works fine on PDFs.
##'
##' @param list A list of objects, each with elements named 'x' and 'y'.
##'
##' @param ... Arguments to \code{matplot}
##'
##' @seealso matplot
##'
##' @examples
##'
##' mu <- rnorm(5, sd=5)
##' sigma <- rgamma(5, shape=10, scale=0.25)
##' x <- mapply(function(m,s){rnorm(10000,m,s)}, mu, sigma, SIMPLIFY=FALSE)
##' d <- lapply(x, density)
##' mplot(d, type="l")
##'
##' @export


mplot <- function(list, ...){
    x <- sapply(list, '[[', "x")
    y <- sapply(list, '[[', "y")
    matplot(x, y, ...)
}
