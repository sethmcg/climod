## Multivariate bias correction script using Cannon's MBC package

suppressMessages(library(devtools))
suppressMessages(load_all("~/climod"))
suppressMessages(library(ncdf4))
suppressMessages(library(zoo))
suppressMessages(library(MBC))

## "ratio sequence": variables that are adjusted with multiplicative
## scaling instead of additive shifting.

## Note: arguably, we should transform tasmax/tasmin to tasmed/dtr,
## with dtr=TRUE

rvar <- c(huss=TRUE, prec=TRUE, rsds=FALSE, tmax=FALSE, tmin=FALSE,
          uas=FALSE, vas=FALSE)

vars <- names(rvar)

## scale variables to range of 10s; needed for precip trace=0.1 not to flatten huss to 0.
scaling <- c(huss=1000, prec=1, rsds=1/10, tmax=1, tmin=1, uas=10, vas=10)


#######

# For testing
args <- c("huss",
          "mbcn/raw/huss.METDATA.NAM-44i/slice/x150/data.x150.y040.nc",
          "mbcn/raw/huss.hist.CanESM2.RCA4.day.NAM-44i.raw/slice/x150/data.x150.y040.nc",
          "mbcn/raw/huss.rcp85.CanESM2.RCA4.day.NAM-44i.raw/slice/x150/data.x150.y040.nc",
          "mbcn/mbcn-METDATA/huss.hist.CanESM2.RCA4.day.NAM-44i.mbcn-METDATA/slice/x150/data.x150.y040.nc",
          "mbcn/mbcn-METDATA/huss.rcp85.CanESM2.RCA4.day.NAM-44i.mbcn-METDATA/slice/x150/data.x150.y040.nc"
          )

## comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)


varstring <- args[1]

infiles   <- list()
outfiles  <- list()

vsub <- function(pat, reps, x, ...){
    result <- c()
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


print(basename(outfiles[["fut"]][1]))

## early bailout: if any input is all bad, all outputs are NA
bail <- function(iflist = infiles, oflist = outfiles){
  for(i in c("cur","fut")){
    ifiles <- iflist[[i]]
    ofiles <- oflist[[i]]
    for(v in names(ifiles)){
      file.copy(ifiles[v], ofiles[v], overwrite=TRUE, copy.mode=FALSE)
      fout <- nc_open(ofiles[v], write=TRUE)
      ncvar_put(fout, v, ncvar_get(fout, v)+NA)        
      nc_close(fout)
    }
  }
  q(save="no")
}


NAfill <- function(x){
  x <- zoo::na.approx(x, na.rm=FALSE)
  x <- zoo::na.locf(x, na.rm=FALSE)
  x <- zoo::na.locf(x, na.rm=FALSE, fromLast=TRUE)
  if(any(is.na(x))){stop("NA values remain after filling")}
  x
}

#  N <- length(x)
#  if(is.na(x[1])) x[1] <- x[2]
#  if(is.na(x[N])) x[N] <- x[N-1]
#  bad <- is.na(x)
#  if(any(bad)){
#    x <- zoo::na.approx(x)
#    ibad <- which(bad)
#    x[ibad] <- (x[ibad-1] + x[ibad+1])/2
#  }
#  return(x)



ingest <- function(files){
    result <- c()
    for(v in names(files)){
        fin <- files[v]
        nc <- nc_open(fin)
        x <- ncvar_get(nc, v)
        if(all(!is.finite(x)) || all(x == 0)){bail()}
        x <- x * scaling[v]
        x <- NAfill(x)
        result <- cbind(result, x)
    }
    colnames(result) <- vars
    return(result)
}

indata <- lapply(infiles, ingest)


#######

gettime <- function(ncfile){
    nc <- nc_open(ncfile)
    time <- ncvar_get(nc, "time")
    time@calendar <- ncatt_get(nc, "time", "calendar")$value
    time@units <- ncatt_get(nc, "time", "units")$value
    return(time)
}


## for disassembled data
timefiles <- sapply(infiles, `[`, 1)
names(timefiles) <- names(infiles)
timefiles <- gsub("slice/.*nc", "time.nc", timefiles)

times <- lapply(timefiles, gettime)



N <- 12

monwin <- lapply(times, cslice, num=N, ratio=2, split=FALSE)

bcdata <- renest(mapply(slice, indata, monwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE))




fits <- list()

cat(date(),"\n")
for(i in 1:N){
    cat(i, "of", N, ": ")
    bcme <- bcdata[[i]]
    fits[[i]] <- MBCn(bcme$obs, bcme$cur, bcme$fut, ratio.seq=rvar, trace=0.1) 
}
cat(date(),"\n")

rfit <- renest(fits)

outdata <- list()
outdata$cur <- unslice(rfit$mhat.c, monwin$cur, template=indata$cur)
outdata$fut <- unslice(rfit$mhat.p, monwin$fut, template=indata$fut)


## write results out to netcdf file

outgest <- function(ifiles, ofiles, datarray){
    for(v in names(ifiles)){
        file.copy(ifiles[v], ofiles[v], overwrite=TRUE, copy.mode=FALSE)
        fout <- nc_open(ofiles[v], write=TRUE)
        ncvar_put(fout, v, datarray[,v]/scaling[v])
        
        ## bias correction metadata gets added to the final result, not
        ## intermediate output
        nc_close(fout)
    }
}

for(i in c("cur","fut")){
    outgest(infiles[[i]], outfiles[[i]], outdata[[i]])
}

