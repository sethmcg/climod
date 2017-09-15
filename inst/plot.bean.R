#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)
library(beanplot)

## Call as: Rscript --vanilla bean.R label obs cur fut out varname

args <- commandArgs(trailingOnly=TRUE)

label <- args[1]
        
infiles <- c()
infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfile <- args[5]
v <- args[6]

# 
# label <- "prec rcp85 GFDL-ESM2M RegCM4 NAM-22i raw boulder 2036-2065"
# 
# infiles <- c()
# infiles["obs"] <- "obs/prec.livneh.boulder.nc"
# infiles["cur"] <- "ts/prec.hist.GFDL-ESM2M.RegCM4.day.NAM-22i.raw.boulder.nc"
# infiles["fut"] <- "ts/fut/prec.rcp85.HadGEM2-ES.WRF.day.NAM-22i.raw.boulder.2036-2065.nc"
# 
# outfile <- "raw.png"
# v <- "prec"
# 
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


## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")

## bean element colors
bcol <- lapply(ocf, function(x){c(adjustcolor(x,0.5),x,x,x)})

## bean elements to plot
what <- c(FALSE, TRUE, TRUE, FALSE)


yr <- range(unlist(cdata))

## plotting

png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(3,1), oma=c(0,0,3,0), mar=c(3,5,3,2)+0.1)

for(i in names(data)){

  beanplot(cdata[[i]], what=what, col=bcol[[i]], ylim=yr,
           names=month.abb, ylab=units, main=paste(i,v,"PDF"))
}

mtext(label, line=1, outer=TRUE)
    
dev.off()


