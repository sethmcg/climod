#library(climod)
library(devtools)
load_all("~/climod")
suppressMessages(library(KernSmooth))
suppressMessages(library(ncdf4))
suppressMessages(library(quantreg))

## Call as: Rscript kddm.R var obs cur fut cout fout
## (used to call with --vanilla but on cheyenne it can't find pkgs...)


# For testing
args <- c("prec",
          "obs/prec.obs.livneh.ftlogan.nc",
          "raw/prec.hist.HadGEM2-ES.WRF.ftlogan.nc",
          "raw/prec.rcp85.HadGEM2-ES.WRF.ftlogan.nc",
          "prec.test.cur.nc",
          "prec.test.fut.nc",
          "prec.test.Rdata"
          )

## comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)


varname <- args[1]

infiles   <- c()
outfiles  <- c()

infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfiles["cur"] <- args[5]
outfiles["fut"] <- args[6]


## If the 7th argument is defined, it's the name of a save file for
## data slices and distmap objects from the observed average min and
## max in the seasonal cycle, the slice containing the peak (single
## greatest value) for cur and fut, and the year end (index 1).
saveslice <- !is.na(args[7])
if(saveslice){
  savefile <- args[7]
}

print(basename(outfiles["fut"]))
        
nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
indata <- lapply(nc,"[[",varname)


## which norm to use depends on the variable
## use standard_name, which is less idiosyncratic than varname

stdname <- (indata$obs)@standard_name
isprec <- grepl("precip", stdname)
istemp <- grepl("temperature", stdname)

standard_norms <- c(
    "precip"                    = "log",
    "relative_humidity"         = "range",
    "welling.*wave_flux_in_air" = "scale",
    "wind_speed"                = "scale",
    "air_pressure"              = "zscore",
    "ward_wind"                 = "zscore",
    "air_temperature"           = "zscore"
    )
norm <- standard_norms[sapply(names(standard_norms), grepl, stdname)]


## if inputs are bad (all NA or obs always 0), output is all NA
if(any(sapply(indata, function(x){all(!is.finite(x))})) ||
   all(indata$obs == 0)){

  for(i in c("cur","fut")){
    file.copy(infiles[i], outfiles[i], overwrite=TRUE, copy.mode=FALSE)
    fout <- nc_open(outfiles[i], write=TRUE)
    ncvar_put(fout, varname, ncvar_get(fout,varname)*NA)
    nc_close(fout)
  }  
  
  quit(save="no", status=0, runLast=FALSE)
}


## Normally there's a time coordinate variable in each file.  But if
## the data has been sharded for parallelism, it gets pulled out into
## its own file that lives alongside all the shards.  Here we attempt
## to determine which is the case and automagically do the right
## thing.

time <- lapply(nc, "[[", "time")

## When sharded, the timefile is named "time.nc"
timefiles <- paste(sapply(infiles, dirname), "time.nc", sep="/")
names(timefiles) <- names(infiles)

for(t in names(time)){
  if(is.null(time[[t]])){    
    ## Note: we do this stuff by hand instead of using nc_ingest
    ## because the time files have no non-dimension vars in them.
    timenc <- nc_open(timefiles[t])
    time[[t]] <- ncvar_get(timenc, "time")
    time[[t]]@units <- ncatt_get(timenc, "time", "units")$value
    time[[t]]@calendar <- ncatt_get(timenc, "time", "calendar")$value
  }   
}


### Need to match up units and epochs
time <- lapply(time, alignepochs, "days since 1950-01-01")


### data structures for storage of results
outdata <- lapply(indata,"+",NA)

## width of moving window
mwinwidth = 30


### generate climatology moving window index arrays
### 360 ~1-day moving windows w/ 30-day outer pool

### (360 is important because we want to be consistent from model to
### model and HadGEM uses a 360-day calendar, which will give you
### empty days if you use 365 windows.)

### Temperature data needs to be segemented into years and each year
### de/normalized separately to handle the trend.  Precip data needs
### to be pooled across years for dedrizzling and to ensure there are
### enough non-zero values to be able to fit the boxcox normalization.
### Other variables don't generally show dramatic trends the way
### temperature does, so we leave them pooled across years.

cwin <- lapply(time, cslice, outer=mwinwidth, num=360, split=istemp)

    
##############################
## Do all the bias corrections

## window data using outer window
wind <- mapply(slice, indata, cwin,
               MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)

if(isprec){

  ## dedrizzle by slice
  ddz <- renest(lapply(renest(wind), dedrizzle))
  ## unzero
  ddu <- rapply(ddz, unzero, how="replace")
  ## dummy list on innermost nesting to match segmented slice structure for temp
  datatobc <- rapply(renest(ddu), list, how="replace")

  bctrunc <- TRUE
  bctrim  <- function(x){ptrim(x, p=0.05)}
  densfun <- akde

} else {

    if(istemp){
        ## just invert list nesting
        datatobc <- renest(wind)
    } else {
        ## dummy inner list to match temp structure as with prec
        datatobc <- rapply(renest(wind), list, how="replace")
    }

    bctrunc <- FALSE
    bctrim  <- NULL
    if(stdname == "wind_speed"){ bctrim <- ptrim }
    densfun <- bkde
}



## bias-correct each window

fixdata <- lapply(datatobc, biascorrect, norm, dmap=saveslice,
                  truncate=bctrunc, trim=bctrim, densfun=densfun,
                  minobs=300)


if(saveslice){
  ## Have to separate out the distmaps, or unslicing breaks
  fixdata <- renest(fixdata)
  dmaps <- fixdata$distmap
  fixdata$distmap <- NULL
  fixdata <- renest(fixdata)
}

## rearrange back to sliced conformation
if(isprec){
    ## rezero data
    fixdata <- rapply(fixdata, rezero, how="replace")
}
if(!istemp){
    ## get rid of the dummy inner list
    fixdata <- lapply(fixdata, function(x){lapply(x, unlist)})
}
## re-invert
bc <- renest(fixdata)

## collate inner windows back into timeseries
outdata <- mapply(unslice, bc, cwin, SIMPLIFY=FALSE)


## Save distmaps and data at peak slices

if(saveslice){

  oraw <- lapply(wind, function(x){lapply(x, unlist)})
  iraw <- mapply(subslice, wind, cwin, split=FALSE, SIMPLIFY=FALSE)
  ifix <- mapply(subslice, bc,   cwin, split=FALSE, SIMPLIFY=FALSE)
  
  peaks <- c()
  peaks["obsmin"]   <- which.min(sapply(oraw$obs, mean, na.rm=TRUE))
  peaks["obsmax"]   <- which.max(sapply(oraw$obs, mean, na.rm=TRUE))
  peaks["rawcurpk"] <- which.max(sapply(iraw$cur, max,  na.rm=TRUE))
  peaks["rawfutpk"] <- which.max(sapply(iraw$fut, max,  na.rm=TRUE))
  peaks["fixcurpk"] <- which.max(sapply(ifix$cur, max,  na.rm=TRUE))
  peaks["fixfutpk"] <- which.max(sapply(ifix$fut, max,  na.rm=TRUE))

  jdays <- round(peaks/length(fixdata)*365.25)
  dates <- format(as.Date(jdays, origin="1950-01-1"), format="%b %d")

  dmaps <- dmaps[peaks]
  names(dmaps) <- names(peaks)

  ## save times (as fractional year) for outer segmented data
  
  toutseg <- mapply(slice, time, cwin,
                    MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)
  
  ## convert time to year (
  for(i in names(toutseg)){
    toutseg[[i]] <- rapply(toutseg[[i]], how="replace",
                           function(x){x/yearlength(time[[i]]@calendar)+1950})
  }
                         
  tslice <- renest(toutseg)[peaks]; names(tslice) <- names(peaks)
  
  
  ## not needed:
  ##  raw outer continuous <- renest(oraw)[peaks]
  ##  fix outer continuous <- renest(ofix)[peaks]

  ## (inner segemented doesn't really make sense)
  
  icraw <- renest(iraw)[peaks]; names(icraw) <- names(peaks)
  icfix <- renest(ifix)[peaks]; names(icfix) <- names(peaks)
  
  osraw <- renest(wind)[peaks]; names(osraw) <- names(peaks)
  osfix <- renest(bc)[peaks]  ; names(osfix) <- names(peaks)

  innerdata <- list(raw=icraw, fix=icfix)
  outerdata <- list(raw=osraw, fix=osfix)

  units <- indata$obs@units

  ## annual maxima
  iodata <- c(indata, outdata)
  iodata[4] <- NULL
  names(iodata) <- c("obs", "rawcur", "rawfut", "fixcur", "fixfut")
  iotime <- c(time,time)
  iotime[4] <- NULL
  names(iotime) <- c("obs", "rawcur", "rawfut", "fixcur", "fixfut")
  year <- lapply(iotime, function(x){floor(x/yearlength(x))})
  ydata <- mapply(split, iodata, year)
  bmax <- lapply(ydata, function(x){sapply(x, max, na.rm=TRUE, USE.NAMES=FALSE)})
  yyear <- lapply(year, function(x){unique(x)+1950})
  annmax <- lapply(mapply(cbind, yyear, bmax), as.data.frame)
  annmax <- lapply(annmax, `colnames<-`, c("year",varname))

  
  save(peaks, jdays, dates, dmaps, innerdata, outerdata, tslice, units,
       annmax, file=savefile)
}




## write results out to netcdf file

for(i in c("cur","fut")){

  file.copy(infiles[i], outfiles[i], overwrite=TRUE, copy.mode=FALSE)
        
  fout <- nc_open(outfiles[i], write=TRUE)
        
  ncvar_put(fout, varname, outdata[[i]])

  ## bias correction metadata gets added to the final result, not
  ## intermediate output
    
  nc_close(fout)
}

