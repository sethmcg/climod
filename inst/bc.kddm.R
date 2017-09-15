#library(climod)
library(devtools)
load_all("~/climod")
suppressMessages(library(KernSmooth))
suppressMessages(library(ncdf4))
suppressMessages(library(quantreg))

## Call as: Rscript --vanilla kddm.mpmd.R var obs cur fut cout fout

args <- commandArgs(trailingOnly=TRUE)

varname <- args[1]

infiles   <- c()
outfiles  <- c()

infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfiles["cur"] <- args[5]
outfiles["fut"] <- args[6]


norms <- c(prec="boxcox", tmax="zscore", tmin="zscore")

norm <- norms[varname]

print(basename(outfiles["fut"]))
        
nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
indata <- lapply(nc,"[[",varname)


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


## A separate time.nc file lives alongside each input file.

## Note: we use sub instead of dirname+paste to preserve the names on
## the array, and ncvar_get, etc. instead of nc_ingest because the
## time files have no non-dimension vars in them.

timefiles <- sub("data\\.x.*\\.y.*\\.nc", "time.nc", infiles)
time <- list()
for(t in names(timefiles)){
  timenc <- nc_open(timefiles[t])
  time[[t]] <- ncvar_get(timenc, "time")
  time[[t]]@units <- ncatt_get(timenc, "time", "units")$value
  time[[t]]@calendar <- ncatt_get(timenc, "time", "calendar")$value
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

cwin <- lapply(time, cslice, outer=mwinwidth, num=360, split=(varname!="prec"))

    
##############################
## Do all the bias corrections

## window data using outer window
wind <- mapply(slice, indata, cwin,
               MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)

if(varname == "prec"){
  ## dedrizzle by slice
  ddz <- renest(lapply(renest(wind), dedrizzle))
  ## unzero
  ddu <- rapply(ddz, unzero, how="replace")
  ## dummy list on innermost nesting to match segmented slice structure for temp
  datatobc <- rapply(renest(ddu), list, how="replace")
} else {
  ## just invert list nesting
  datatobc <- renest(wind)
}


## bias-correct each window
fixdata <- lapply(datatobc, biascorrect, norm, truncate=(varname=="prec"))


## rearrange back to sliced conformation
if(varname == "prec"){
  ## rezero data, get rid of the dummy inner list, and re-invert
  rez <- rapply(fixdata, rezero, how="replace")  
  bc <- renest(lapply(rez, function(x){lapply(x, unlist)}))
} else {
  ## just re-invert
  bc <- renest(fixdata)
}


## collate inner windows back into timeseries
outdata <- mapply(unslice, bc, cwin, SIMPLIFY=FALSE)


## write results out to netcdf file

for(i in c("cur","fut")){

  file.copy(infiles[i], outfiles[i], overwrite=TRUE, copy.mode=FALSE)
        
  fout <- nc_open(outfiles[i], write=TRUE)
        
  ncvar_put(fout, varname, outdata[[i]])

  ## bias correction metadata gets added to the final result, not
  ## intermediate output
    
  nc_close(fout)
}
