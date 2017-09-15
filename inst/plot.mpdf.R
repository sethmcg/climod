#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript --vanilla mpdf.R label obs cur fut out varname

args <- commandArgs(trailingOnly=TRUE)

label <- args[1]
        
infiles <- c()
infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfile <- args[5]
v <- args[6]

# label <- "tmax rcp85 GFDL-ESM2M RegCM4 NAM-22i kddm-livneh boulder 2036-2065"
# 
# infiles <- c()
# infiles["obs"] <- "obs/tmax.livneh.boulder.nc"
# infiles["cur"] <- "ts/tmax.hist.HadGEM2-ES.WRF.day.NAM-22i.kddm-livneh.boulder.nc"
# infiles["fut"] <- "ts/fut/tmax.rcp85.HadGEM2-ES.WRF.day.NAM-22i.kddm-livneh.boulder.2036-2065.nc"
# 
# outfile <- "test.png"
# v <- "tmax"


nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
data <- lapply(nc,"[[", v)

units <- data$cur@units
    
time <- lapply(nc,"[[","time")
time <- lapply(time, alignepochs, "days since 1950-01-01")


## Log precip
if(v == "prec"){
  units <- paste0("log(",units,")")
  suppressWarnings(data <- lapply(data, log10))
}


## generate monthly climatology window index arrays
dwin <- lapply(time, cslice, ratio=1, num=12, split=FALSE)

## slice data
cdata <- mapply(slice, data, dwin, SIMPLIFY=FALSE)


## drop non-finite & sub-trace values
trace <- -7   # log scale
cdata <- rapply(cdata, function(x){x[is.finite(x) & x>trace]}, how="replace")

## rearrange
mdata <- renest(cdata)

## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")

## plotting limits
extend <- 3*stats::bw.nrd(unlist(mdata))
    xr <- extend*c(-1,1)+range(unlist(mdata))

## plotting

png(outfile, units="in", res=120, width=10, height=7)

par(mfrow=c(3,4), oma=c(0,0,3,0), mar=c(2.5,1.5,3,1.5)+0.1)

mdata <- renest(cdata)
names(mdata) <- month.abb

for(m in month.abb){
  mplot(lapply(mdata[[m]], akde, min.x=min(xr), max.x=max(xr)),
        yaxt="n", col=ocf, type='l', lty=1, main=m, xlab=units, ylab="")  
}

lloc <- ifelse(v=="prec", "topleft", "topright")
legend(lloc, names(ocf), col=ocf, lty=1, seg.len=1)

mtext(paste(label, "PDFs (",units,")"), line=1, outer=TRUE)
    
dev.off()
