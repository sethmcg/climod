#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)
library(beanplot)
library(moments)

## Call as: Rscript plot.bean.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("tmin rcp85 HadGEM2-ES WRF birmingham",
          "obs/tmin.obs.livneh.birmingham.nc",
          "save-test/tmin.hist.HadGEM2-ES.WRF.birmingham.nc",
          "save-test/tmin.rcp85.HadGEM2-ES.WRF.birmingham.nc",
          "test.bean.png",
          "tmin",
          "test.bean.txt"
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

nc <- lapply(infiles, nc_ingest)

## extract variables of interest from the netcdf objects
data <- lapply(nc,"[[", v)

units <- data$cur@units
    
time <- lapply(nc,"[[","time")
time <- lapply(time, alignepochs, "days since 1950-01-01")


## generate monthly climatology window index arrays
cwin <- lapply(time, cslice, ratio=1, num=12, split=FALSE)

## slice data
cdata <- mapply(slice, data, cwin, SIMPLIFY=FALSE)


## Log precip (for plotting only)
if(v == "prec"){
  units <- paste0("log(",units,")")
  suppressWarnings(ddata <- rapply(cdata, log10, how="replace"))
  ## trace <- -7   # log scale
  ## ddata <- rapply(ddata, function(x){x[x>trace]}, how="replace")
  ## # Not clipping trace because we want proper min value in metrics
} else {
  ddata <- cdata
}


ddata <- rapply(ddata, function(x){x[is.finite(x)]}, how="replace")


## color palette for bias-correction
ocf <- c(obs="black", cur="blue", fut="red")

## bean element colors
bcol <- lapply(ocf, function(x){c(adjustcolor(x,0.5),x,x,x)})

## bean elements to plot
what <- c(FALSE, TRUE, TRUE, FALSE)


yr <- range(unlist(ddata))

## plotting

png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(3,1), oma=c(0,0,3,0), mar=c(3,5,3,2)+0.1)

for(i in names(data)){

  beanplot(ddata[[i]], what=what, col=bcol[[i]], ylim=yr,
           names=month.abb, ylab=units, main=paste(i,v,"PDF"))
}

mtext(label, line=1, outer=TRUE)
    
dev.off()


## metrics

## Metrics come from raw (non-log precip) We want to operate on raw untrimmed data Don't use log data for 

# ## Re-trim in case NA values exist
# mdata <- rapply(cdata, function(x){x[is.finite(x)]}, how="replace")

summarize <- function(x){c(mean=mean(x, na.rm=TRUE),
                           stdev=sd(x, na.rm=TRUE),
                           skew=skewness(x, na.rm=TRUE),
                           min=min(x, na.rm=TRUE),
                           max=max(x, na.rm=TRUE))}


p <- c(10, 25, 50, 75, 90)/100
pname <- paste0("q", p*100)

m1 <- lapply(rapply(cdata, summarize, how="replace"), as.matrix)
m1 <- lapply(m1, t)
m2 <- lapply(rapply(ddata, quantile, probs=p, how="replace"), as.matrix)
m2 <- lapply(m2, function(x){rownames(x)<-pname; t(x)})


## Empty metrics dataframe

metrics <- data.frame(infile = "dummy", period="mon", analysis="bean",
                      mean=0, stdev=0, skew=0, min=0, max=0,
                      q10=0, q25=0, q50=0, q75=0, q90=0,
                      stringsAsFactors=FALSE)


for(p in names(data)){
  metrics <- rbind(metrics, data.frame(infile=infiles[p],
                                       period=month.abb,
                                       analysis="bean",
                                       m1[[p]],
                                       m2[[p]],
                                       row.names=NULL))
}



## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
