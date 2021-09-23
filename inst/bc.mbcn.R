## Multivariate bias correction script using Cannon's MBC package

suppressMessages(library(devtools))
suppressMessages(load_all("~/climod"))
suppressMessages(library(ncdf4))
suppressMessages(library(zoo))
suppressMessages(library(MBC))

## TODO: add option to transform tmax/tmin to tmed = mean(tmax,tmin)
## and dtr = diff(tmax,tmin) (and undo it on readout)

## TODO: add option to declare scaling from command-line

## "ratio sequence": variables that are adjusted with multiplicative
## scaling instead of additive shifting.

rvar <- c(huss=TRUE, prec=TRUE, rsds=FALSE, tmax=FALSE, tmin=FALSE,
          uas=FALSE, vas=FALSE)

## scale variables to range of 10s; needed for precip trace=0.1 not to flatten huss to 0.
scaling <- c(huss=1000, prec=1, rsds=1/10, tmax=1, tmin=1, uas=10, vas=10)


#######

## Call as: Rscript bc.mbcn.R varlist obs cur fut cout fout [seed]

# For testing
args <-c("huss,prec,rsds,tmax,tmin,uas,vas",
         "data/obs/huss.obs.gridMET.day.NAM-44i.nc",
         "data/raw/huss.hist.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.raw.nc",
         "data/raw/huss.rcp85.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.raw.nc",
         "data/fix/huss.hist.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.mbcn-gridMET.nc",
         "data/fix/huss.rcp85.MPI-ESM-MR.CRCM5-UQAM.day.NAM-44i.mbcn-gridMET.nc"
         )
verbose <- TRUE


## comment out these lines for testing
args <- commandArgs(trailingOnly=TRUE)
verbose <- FALSE


vars <- unlist(strsplit(args[1], ","))
## first variable should be the one used in the filename args
v1 <- vars[1]

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


infiles[["obs"]] <- vsub(v1, vars, args[2])
infiles[["cur"]] <- vsub(v1, vars, args[3])
infiles[["fut"]] <- vsub(v1, vars, args[4])

outfiles[["cur"]] <- vsub(v1, vars, args[5])
outfiles[["fut"]] <- vsub(v1, vars, args[6])

## set RNG seed used by mbcn functions
if(length(args) == 7){
  set.seed(as.numeric(args[7]) %% 2^31)
} else {
  set.seed(222)
}

if(length(args) > 7){
    warning("Ignoring superfluous arguments: ", paste(args[-(1:7)], collapse=" "))
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



ingest <- function(files){
    result <- c()
    for(v in names(files)){
        fin <- files[v]
        nc <- nc_open(fin)
        x <- ncvar_get(nc, v)
        if(all(!is.finite(x)) || all(x == 0, na.rm=TRUE)){bail()}
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


## get times and convert to fractional year

times <- lapply(timefiles, gettime)

atime <- lapply(times, alignepochs, "days since 1900-01-01")
years <- lapply(atime, function(x){x/yearlength(x)+1900})
yfrac <- lapply(years, '%%', 1)



## construct indexing for ~monthly windowing, plus an outer window of
## 2 months on either side (3 months total) for better statistics.

N <- 12
mbounds <- seq(0, 1, length=N+1)

mbins <- lapply(yfrac, .bincode, breaks=mbounds, include.lowest=TRUE) 

innerwindow <- lapply(mbins, function(x){outer(x, seq(N), '==')})

before <- c(N, seq(N-1))
after  <- c(seq(N-1)+1, 1)

outerwindow <- lapply(innerwindow, function(x){x[,before] | x | x[,after]})


## TODO: automate and parameterize splitting model dataq into chunks
## that match the obs data in size (and start/end, for 'current'.)

## indexing to split data for BC into 4 x 38-year chunks, plus a
## 38-year chunk ("current") matching the obs period

period <- c("late20c", "turncty", "mid21c", "end21c", "current")
ystart <- c(1949, 1987, 2025, 2063, 1979)
yend   <- c(1986, 2024, 2062, 2101, 2017)
names(ystart) <- names(yend) <- period


## concatenate model data, then re-split

y <- c(years$cur, years$fut)

iper <- mapply(function(a, b){a<y & y<=b}, ystart, yend)
adata <- rbind(indata$cur, indata$fut)

## lapply is very slow; use loop instead
perdata <- list()
for(p in period){
    perdata[[p]] <- adata[iper[,p],]
}
perdata$obs <- indata$obs
                           

### quick check:
# par(mfrow=c(6,1))
# for(p in period){
#     plot(y[iper[,p]], perdata[[p]][,"tmax"], xlim=c(1948, 2102), pch='.', main=p)
# }
# plot(years$obs, perdata$obs[,"tmax"], xlim=c(1948, 2102), pch='.', main="obs")


## and do the same for the windowing indexes

ainnerwindow <- rbind(innerwindow$cur, innerwindow$fut)
aouterwindow <- rbind(outerwindow$cur, outerwindow$fut)

pinnerwin <- list()
pouterwin <- list()
for(p in period){
    pinnerwin[[p]] <- ainnerwindow[iper[,p],]
    pouterwin[[p]] <- aouterwindow[iper[,p],]
}

pinnerwin$obs <- innerwindow$obs
pouterwin$obs <- outerwindow$obs



## empty arrays to hold  results
bcdata <- lapply(perdata, function(x){x+NA})
bcdata$obs <- bcdata$current <- NULL


### then loop over both period and window

if(verbose) cat(date(),"\n")
for(p in head(period, -1)){
    ocf <- c("obs","current",p)
    if(verbose) cat(p, ":\n")
    for(i in 1:N){
        if(verbose) cat(i, "of", N, ": ")
        bcme <- mapply(function(D, mask){D[mask[,i],]},
                       perdata[ocf], pouterwin[ocf], SIMPLIFY=FALSE)
        result <- suppressWarnings(
            MBCn(bcme$obs, bcme$current, bcme[[p]],
                 ratio.seq=rvar, trace=0.1, silent=!verbose)
            )

        # indexes of inner window for month i, period p
        idex <- which(pinnerwin[[p]][,i])

        # indexes of inner window within outer window for month i, period p
        odex <- which(pinnerwin[[p]][pouterwin[[p]][,i],i])

        bcdata[[p]][idex,] <- result$mhat.p[odex,]
    }
}
if(verbose) cat(date(),"\n")

# run time for bc loop: 105 seconds on my mac 


## combine 4x38-year chunks back into cur & fut
odata <- do.call(rbind, bcdata)

outdata <- list()
outdata$cur <- head(odata, nrow(indata$cur))
outdata$fut <- tail(odata, nrow(indata$fut))



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
