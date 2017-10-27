#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)

## Call as: Rscript plot.intsp.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("rcp85 HadGEM2-ES RegCM4 birmingham",
          "obs/prec.obs.livneh.birmingham.nc",
          "save-test/prec.hist.HadGEM2-ES.RegCM4.birmingham.nc",
          "save-test/prec.rcp85.HadGEM2-ES.RegCM4.birmingham.nc",
          "test.intsp.png",
          "prec",
          "test.intsp.txt"
          )

## comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)

label <- args[1]
        
infiles <- c()
infiles["obs"] <- args[2]
infiles["cur"] <- args[3]
infiles["fut"] <- args[4]

outfile <- args[5]

v <- args[6]

txtfile <- args[7]


## color palette
cmap <- c(obs="black", cur="blue", fut="red")

        
nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
data <- lapply(nc,"[[",v)
units <- data$cur@units

xlab <- paste0("intensity (", units, ")")

time <- lapply(nc,"[[","time")


## The first cslice window starts at the epoch, so setting it to 12-01
## instead of 01-01 gives DJF/MAM/JJA/SON seasons, rather than
## JFM/AMJ/JAS/OND.

time <- lapply(time, alignepochs, "days since 1949-12-01")


### generate climatology moving window index arrays
### 4 seasonal windows
swin <- lapply(time, cslice, ratio=1, num=4, split=FALSE)
sdata <- renest(mapply(slice, data, swin, SIMPLIFY=FALSE))
names(sdata) <- c("DJF","MAM","JJA","SON")


## Data from one gridcell is much noisier than data over a region;
## need to trim extreme values and use P-B (bigger) binning instead of
## FD (default) to keep the noise from distorting things.

## Drop zero values
unzdata <- rapply(sdata, how="replace", function(x){x[x>0]})

## Trim upper 99%ile values
trim <- function(x){x[x<quantile(x,0.99,na.rm=TRUE)]}
zdata <- rapply(unzdata, trim, how="replace")

## Panofsky-Brier recommendation for number of histogram bins
PB <- function(x){
  pretty(x, n=5*log10(length(x)))
}



## Empty metrics dataframe

metrics <- data.frame(infile = "dummy", period="seas", analysis="intsp",
                      pislope=0,
                      stringsAsFactors=FALSE)

## plotting

png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(2,2), oma=c(0,0,3,0), mgp=c(2,1,0), mar=c(4,3.5,2.5,1))

for(s in names(zdata)){

  p <- zdata[[s]]

  h <- lapply(p, hist, breaks=PB, plot=FALSE)  

  intsp <- lapply(h,
                  function(x){
                    freq <- log10(x$counts/sum(x$counts))
                    good <- is.finite(freq)
                    list(int=x$mids[good], freq=freq[good])
                  })
  
  mplot(intsp, x="int", y="freq", col=cmap, type="p", pch=1,
        ylab="log fraction", xlab=xlab, main=s)
  
  ## Fit a line
  fits <- lapply(intsp, function(x){lm(freq ~ int, x)})

  for(f in names(fits)){
    abline(fits[[f]], col=cmap[f])
    metrics <- rbind(metrics, list(infiles[f], s, "intsp", coef(fits[[f]])[2]))
  }
  
  abline(v=quantile(p$obs,0.5),col="gray")  
}

mtext(label, line=1, outer=TRUE)


## centered legend outside plots

par(fig=c(0,1,0,1), mar=c(0,0,0,0), oma=c(0,0,0,0), new=TRUE)
plot(c(0,1), c(0,1), type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)
legend(0.49, 0.51, names(cmap), col=cmap, pch=1)

dev.off()


## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
