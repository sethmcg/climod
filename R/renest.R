##' Invert the nesting of a list-of-lists
##'
##' Given a nested list, invert the nesting so that the inner list
##' becomes the outer list and vice-versa.  The structures of the
##' elements of the inner lists are preserved.  This function is more
##' easily demonstrated than explained; see examples.
##' 
##' All inner lists must have the same length.  Inner lists that are
##' longer than the first inner list will be truncated; inner lists
##' that are shorter than the first inner list will cause an error.
##' Names of inner lists after the first are ignored.
##' 
##' @param lol a nested list.
##'
##' @return A nested list
##'
##' @examples
##'
##' 
##' @export


renest <- function(lol) {
    apply(do.call(rbind, lol), 2, as.list)
}


f <- list(a=list(x=1, y=2, z=3), b=list(x="one", y="two", z="three"))
g <- list(x=list(a=1, b="one"), y=list(a=2, y="two"), z=list(a=3,b="three"))


print(str(f))
print(str(g))          
print(identical(g, renest(f)))


f1 <- list(a=list(x=1, y=2, z=3),
          b=list(x=identity, y=seq(2), z=lapply(f,length)))

          
g1 <- list(x=list(a=1, b=identity),
          y=list(a=2, b=seq(2)),
          z=list(a=3, b=lapply(f,length)))
        
h1 <- renest(f1)

print(str(f))
print(str(g))          
print(identical(g, h))



fextra <- list(a=list(x=1, y=2, z=3),
          b=list(x="one", y="two", z="three", d="dummy"))

str(fextra)
str(renest(fextra))

fmissing <- list(a=list(x=1, y=2, z=3, e="extra"),
          b=list(x="one", y="two", z="three"))

str(fmissing)
str(renest(missing))

forder <- list(a=list(x=1, y=2, z=3),
          b=list(x="one", z="three", y="two"))

str(forder)
str(renest(forder))



fouter <- list(list(x=1, y=2, z=3),
          list(x="one", z="three", y="two"))

str(fouter)
str(renest(fouter))



finner <- list(a=list(1, 2, 3),
          b=list("one", "two", "three"))

str(finner)
str(renest(finner))



fboth <- list(list(1, 2, 3),
          list("one", "two", "three"))

str(fboth)
str(renest(fboth))



### Copyright 2015 Univ. Corp for Atmos. Research
### Author: Seth McGinnis, mcginnis@ucar.edu
