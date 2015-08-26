##' Copy all attributes from one variable to another
##'
##' This function sets the attributes of the argument \code{to} to be
##' the same as the attributes of \code{from}.  Any existing
##' attributes are replaced.
##'
##' @param from The variable to copy attributes from.
##'
##' @param to The variable whose attribute are to be set.
##'
##' @return A copy of \code{to} with the attributes of \code{from}.
##'
##' @examples
##' x <- 3
##' x@@foo <- "things"
##' y <- 5
##' y@@bar <- "more things"
##' y@@baz <- "other things"
##' z <- copyatts(y,x)
##' z
##' 
##' @export

copyatts <- function(from, to){
    attributes(from) -> attributes(to)
    return(to)
}

