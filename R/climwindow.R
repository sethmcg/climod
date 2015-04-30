##' generate climatology moving window
##'
##' creates moving window indexes for climatology iwth inner and outer
##'
##' 


climwindow <- function(time,
                       num = ceiling(year / inner),
                       ratio = outer / inner,
                       inner = year / num,
                       outer = inner * ratio,
                       year  = yearlength(time)
                       ){

    ## check for:
    ## inner > outer
    ## num > year
    ## all finite, positive
    
    N = length(time)

    span <- diff(range(time))
    cycles <- span / year
  
    window <- list(    
        inner  = vector(mode="list", length=num),
        outer  = vector(mode="list", length=num),
        io     = vector(mode="list", length=num),
        params = namelist(span, cycles, year, ratio, inner, outer)
        )


    time <- (time - time[1])

    for(i in 1:num){
        window$inner[[i]] <- which( (time - inner*(i-1)) %% year < inner)
        window$outer[[i]] <- which( (time - inner*(i-1/2) + outer/2 ) %% year < outer)
        window$io[[i]]    <- match(window$inner[[i]], window$outer[[i]])
    }
    
    return(window)
}



### needed by apply.window.mask, below

indexof <- function(index,x){
    return(x[index])
}

