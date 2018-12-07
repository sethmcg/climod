library(devtools)
load_all("~/climod")
suppressMessages(library(ncdf4))
suppressMessages(library(KernSmooth))

## Call as: Rscript bc.gyro.R var obs cur fut cout fout


vars <- c("huss", "prec", "rsds", "tmax", "tmin", "uas", "vas")

## Note: consider transforming tmax/tmin to tmed/dtr to guarantee tmax >= tmin

# For testing
args <- c("huss",
          "big7-obs/huss.METDATA.22i.armsite.nc",
          "big7-raw/huss.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.armsite.nc",
          "big7-raw/huss.rcp85.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.armsite.nc",
          "big7-gyro/huss.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.gyro.armsite.nc",
          "big7-gyro/huss.rcp85.GFDL-ESM2M.RegCM4.day.NAM-22i.gyro.armsite.nc"
          )

## comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)


####################
### Argument parsing

varstring <- args[1]

infiles   <- list()
outfiles  <- list()

vsub <- function(pat, reps, x, ...){
        result <- list()
    for(r in reps){
        result <- c(result, gsub(pat, r, x, ...))
    }
    names(result) <- reps
    result
}


infiles[["obs"]] <- vsub(varstring, vars, args[2])
infiles[["cur"]] <- vsub(varstring, vars, args[3])
infiles[["fut"]] <- vsub(varstring, vars, args[4])

outfiles[["cur"]] <- vsub(varstring, vars, args[5])
outfiles[["fut"]] <- vsub(varstring, vars, args[6])

if(length(args) > 6){
    warning("Ignoring superfluous arguments: ", paste(args[-(1:6)], collapse=" "))
}


print(basename(outfiles$fut[[1]]))


#####################
### Ingest input data

innc <- rapply(infiles, nc_ingest, how="replace")

rawdata <- renest(innc)
for(v in vars){
    rawdata[[v]] <- lapply(rawdata[[v]], `[[`, v)
}



## bailout: if any input is all bad, all outputs are NA

allbad <- rapply(rawdata, function(x){all(!is.finite(x))})

## Copy infile to outfile, set all values of vname to NA
nacopy <- function(infile, outfile, vname){
    file.copy(infile, outfile, overwrite=TRUE, copy.mode=FALSE)
    fout <- nc_open(outfile, write=TRUE)
    ncvar_put(fout, vname, ncvar_get(fout, vname)+NA)
    nc_close(fout)
}

bailout <- function(msg, infiles=infiles, outfiles=outfiles, vars=vars){
    warning("Bailing out: ", msg)
    for(i in c("cur","fut")){
        for(v in vars){
            nacopy(infiles[[i]][[v]], outfiles[[i]][[v]], v)
        }
    }
    q(save="no")    
}

if(any(allbad)){ bailout("all values bad in ", paste(names(allbad)[allbad], collapse=", ")) }

########################
### Get time coordinates
### Todo: change this to autochecking for time in external sharded file, as in bc.kddm.R 

time <- lapply(innc, function(x){lapply(x, `[[`, "time")})

## check that they all have the same temporal coverage

time <- rapply(time, alignepochs, epoch="days since 1950-01-01", how="replace")

samevec <- function(x){
    diff(range(x)) == 0
}

checktimeok <- function(tlist){
    tmat <- matrix(unlist(tlist), ncol=length(tlist), byrow=FALSE)
    all(apply(tmat, 1, samevec))
}

oktime <- sapply(time, checktimeok)

if(!all(oktime)){
    bailout(paste("Times not identical in", paste(names(oktime)[!oktime], collapse = ", ")))
}

## Only need one copy once we've verified they're all the same

time <- lapply(time, `[[`, 1)


#### norm & detrending code moved from here


###################################
### Window data for bias correction

#indata <- dtdata


### Set up moving windows

### generate climatology moving window index arrays
### 36 ~10-day moving windows w/ 30-day outer pool

### (360 is important because we want to be consistent from model to
### model and HadGEM uses a 360-day calendar, which will give you
### empty days if you use 365 windows.)

### Normaly I'd do 360 1-day windows, but with the extra iterations
### it's slow enough that a factor of 10 speedup is valuable.


## width of moving window
mwinwidth = 30
Nwin <- 36

cwin <- lapply(time, cslice, outer=mwinwidth, num=Nwin, split=FALSE)

wintime <- renest(mapply(slice, time, cwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE))

## window data using outer window
windata <- list()
for(v in vars){
#    windata[[v]] <- renest(mapply(slice, indata[[v]], cwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE))
    windata[[v]] <- renest(mapply(slice, rawdata[[v]], cwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE))
}



########################
### Bias correction code

### Random rotation matrix of rank k
### Adapted from Alex Cannon's MBC library under GPL
randrot <- function(k) {
  qrand <- qr(matrix(rnorm(k * k), ncol=k))
  dR <- diag(qr.R(qrand))
  rot <- qr.Q(qrand) %*% diag(dR/abs(dR))
  return(rot)
}

## list to matrix
L2M <- function(L){
    do.call(cbind,L)
}

## matrix to list
M2L <- function(M){
    lapply(seq(ncol(M)), function(i){M[,i]})
}

## apply rotation matrix
rotate <- function(dlist, R){
    M2L(L2M(dlist) %*% R)
}


## rank tranformation for marginal-copula decomposition
rankxform <- function(x){rank(x)/(length(x)+1)}


## add tails to "pin" the transfer function outside the [0,1] copula
## range.  This is bit of a hack; the really right thing to do is to
## use beta function kernels for the PDF estimation, but I haven't
## figured out how to implement that yet.  This works well enough for
## now.  Adding tails that are the same size as x ensures that the
## pinning starts at the same percentiles when we build the
## distribution mapping.

addtails <- function(x, lower=c(-1,0), upper=c(1,2)){
    lotail <- seq(min(lower), max(lower), along=x)
    hitail <- seq(min(upper), max(upper), along=x)
    c(lotail, x, hitail)
}


## copula bias correction for innermost loop
cbc <- function(ocf, ...){
    mapping <- distmap(addtails(ocf$cur), addtails(ocf$obs), ...)
    result <- lapply(ocf, function(x){predict(mapping, x)})
    result$obs <- ocf$obs
    return(result)
}


## Pitie's n-pdft algorithm
npdft <- function(copula, niter=600, verbose=TRUE, progress=10, ...){
#npdft <- function(copula, niter=10, verbose=TRUE, progress=1, vnames=vars, ...){

    k <- length(vnames)

    x <- copula
    
    for(i in 1:niter){
        ## random rotation matrix
        R <- randrot(k)
        
        ## rotate copula
        y <- lapply(x, rotate, R)
        
        ## bias-correct rotated copula
        z <- renest(lapply(renest(y), cbc, ...))
        
        ## undo rotation
        x <- lapply(z, rotate, t(R))
        
        if(verbose && i %% progress == 0 ){cat(".")}
    }

    lapply(x, function(xx){names(xx) <- vnames; xx})
}




## Calculate trends for obs, cur, & fut.  Need to do a broken-stick
## (continuous piecewise linear) regression for model data or
## discontinuity between trends eats the climate change signal.

ocftrend <- function(ocf, time){
    result <- list()
    result$obs <- predict(lm(ocf$obs ~ time$obs))

    ## fragile; should check that cur & fut are properly aligned
    knot <- mean(max(time$cur), min(time$fut))
    t1 <- c(time$cur, time$fut)
    t2 <- pmax(t1 - knot, 0)
    x <- c(ocf$cur, ocf$fut)
    model <- lm(x ~ t1 + t2)

    result$cur <- predict(model)[t1 <= knot]
    result$fut <- predict(model)[t1 > knot]

    return(result)
}




## Copula gyration bias correction: Do a marginal-copula decomposition
## using a rank transformation, bias-correct the copula alone using
## Pitie's nPDF-t algorithm, then convert the corrected copula values
## back through the observed marginals.  Conceptually, this basically
## injects the entire npdft rotation sequence into the middle of the
## KDDM mapping.

## Input data.ocf.mv if a nested list with ocf outside, variables inside

cogyro <- function(data.ocf.mv, ...){

    k <- length(vars)

    copula <- rapply(data.ocf.mv, rankxform, how="replace")

    copula <- npdft(copula, ...)
    
    ## Some contraction happens in the n-pdft process; fix it by
    ## reapplying the rank transformation.

    ## Probably this is happening because I'm using a gaussian kernel with pinned
    ## tails for the distribution mapping, and the correct way to go about
    ## it is to use a beta kernel limited to [0,1], but I don't have time
    ## to figure out how to do that right now.

    copula <- rapply(copula, rankxform, how="replace")

    ## Now convert rank data back to values via the obs mapping

    result <- list()
    for(v in vars){
        bctrim <- ifelse(v == "prec", ptrim, identity)
        mapping <- distmap(copula$obs[[v]], data.ocf.mv$obs[[v]], trim=bctrim)
        result[[v]] <- lapply(renest(copula)[[v]], function(x){predict(mapping, x)})
        result[[v]]$obs <- data.ocf.mv$obs[[v]]
    }

    result <- renest(result)
}





## multivariate bias correction

## data.mv.ocf is a nested list with variables outer, ocf inner. Time is ocf.

mvbc <- function(data.mv.ocf, time, gridscale=25, ...){

    vars <- names(data.mv.ocf)

    ## save the original for later
    praw <- data.mv.ocf$prec
    
    ## Apply subtrace replacement & log norm to precip
    ## Subtrace uses runif(); set RNG seed first
    data.mv.ocf$prec <- lapply(data.mv.ocf$prec, subtrace, gridscale=gridscale)
    data.mv.ocf$prec <- lapply(data.mv.ocf$prec, log)

    ## Calculate trend vectors
    trends <- lapply(data.mv.ocf, ocftrend, time)    

    ## Detrend data.  (Note: calculating the trend doesn't get folded
    ## into removing it, because we need to put it back at the end of
    ## the bias correction.

    dtdata <- list()
    for(v in vars){
        dtdata[[v]] <- mapply(`-`, data.mv.ocf[[v]], trends[[v]], SIMPLIFY=FALSE)
    }

    ## Do multivariate copula gyration bias correction
    bcdata <- renest(cogyro(renest(dtdata), ...))


    ## calculate bias deltas beween obs and cur trends.  

    ## Figure out the overlap period:
    tstart  <- max(min(time$obs), min(time$cur))
    tfinish <- min(max(time$obs), max(time$cur))

    ## This is not exact because of calendar issues, but good enough for now
    obsmask <- tstart <= time$obs & time$obs <= tfinish
    curmask <- tstart <= time$cur & time$cur <= tfinish

    ## calculate bias delta
    delta <- c()
    for(v in vars){
        delta[v] <- mean(trends[[v]]$cur[curmask]) - mean(trends[[v]]$obs[obsmask])
    }

    ## remember, *subtract* delta, don't add it    

    ## retrend data & subtract trend bias

    redata <- list()
    for(v in vars){
        redata[[v]] <- mapply(`+`, bcdata[[v]], trends[[v]], SIMPLIFY=FALSE)
        for(p in c("cur","fut")){
            redata[[v]][[p]] <- redata[[v]][[p]] - delta[v]
        }
        redata[[v]]$obs <- data.mv.ocf[[v]]$obs
    }

    ## final fixup steps    
    fixdata <- redata

    ## ensure that rsds & huss are non-negative
    fixdata$huss <- lapply(fixdata$huss, pmax, 0)
    fixdata$rsds <- lapply(fixdata$rsds, pmax, 0)

    
    ## Fix up precip, which takes extra work    
    p <- fixdata$prec
    
    ## apply trace threshold to get wet/dry proportions right
    traceval <- p$obs@traceval

    ## percentage of non-trace precip in obs
    pwet <- sum(praw$obs >= traceval) / length(praw$obs)
    
    
    ## do a final bias-correction pass applied only to non-trace values
    
    ## bc mapping for non-trace-precip
    ntp <- lapply(p, function(x){x[rank(x) > pwet*length(x)]})
    pmapping <- distmap(ntp$cur, ntp$obs, trim=ptrim)
    
    ## we don't care if it changes sub-trace values, because we're
    ## about to discard them.  prec$obs also gets mangled, but we
    ## overwrite it down below.
    fixdata$prec <- lapply(fixdata$prec, function(x){predict(pmapping, x)})
    
    ## undo log
    fixdata$prec <- lapply(fixdata$prec, exp)
    
    ## retrace
    fixdata$prec <- lapply(fixdata$prec, function(x){x[x < traceval] <- 0; x})

    ## put obs back to original
    fixdata$prec$obs <- praw$obs
    
    ## [if this doesn't work, use ptrace, subtract zero, add traceval]
    
    ## scale mean precip
    pmean <- lapply(fixdata$prec, mean)
    pscale <- pmean$obs / pmean$cur
    fixdata$prec$cur <- fixdata$prec$cur * pscale
    fixdata$prec$fut <- fixdata$prec$fut * pscale

    ## done!
    return(fixdata)
}




#########################
### Apply bias-correction

## reorder to [Nwin, vars, ocf]
datatobc <- renest(windata)


## set RNG seed; used by subtrace() and cogyro() within mvbc()
set.seed(222)


## apply bias-corrections by window
bcdata <- list()
for(i in 1:Nwin){
    cat(i)
    bcdata[[i]] <- mvbc(datatobc[[i]], wintime[[i]], densfun=bkde)
    cat("\n")
}


## rearrange and collate inner windows back into timeseries
unwindata <- renest(bcdata)
outdata <- list()
for(v in vars){
    outdata[[v]] <- mapply(unslice, renest(unwindata[[v]]), cwin, SIMPLIFY=FALSE)
}


save.image("testing.Rdata")

#stop()

#########################
### Write results to file

## write results out to netcdf file

outgest <- function(ifiles, ofiles, datalist){
    for(v in vars){
        file.copy(ifiles[[v]], ofiles[[v]], overwrite=TRUE, copy.mode=FALSE)
        fout <- nc_open(ofiles[[v]], write=TRUE)
        ncvar_put(fout, v, datalist[[v]])
        
        ## bias correction metadata gets added to the final result, not
        ## intermediate output
        nc_close(fout)
    }
}

for(i in c("cur","fut")){
    outgest(infiles[[i]], outfiles[[i]], renest(outdata)[[i]])
}

