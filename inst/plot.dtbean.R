#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)
library(beanplot)
library(moments)
library(PCICt)
library(udunits2)

## Call as: Rscript plot.dtbean.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

# for testing
args <- c("tmax rcp85 HadGEM2-ES RegCM4 birmingham",
          "big7-obs/tmax.METDATA.22i.birmingham.nc",
          "big7-raw/tmax.hist.HadGEM2-ES.RegCM4.day.NAM-22i.raw.birmingham.nc",
          "big7-raw/tmax.rcp85.HadGEM2-ES.RegCM4.day.NAM-22i.raw.birmingham.nc",
          "test.dtbean.png",
          "tmax",
          "test.dtbean.txt"
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


## convert netcdf time to posix representation via PCICt
nctime2posix <- function(nctime){
    calendar <- nctime@calendar
    ustring  <- strsplit(nctime@units, " since ")[[1]]
    units    <- ustring[1]
    origin   <- ustring[2]
    
    ptime <- as.PCICt(ud.convert(nctime, units, "sec"), calendar, origin)
    as.POSIXlt(ptime)
}

ptime <- lapply(time, nctime2posix)


## Log precip (for plotting only)
if(v == "prec"){
  units <- paste0("log(",units,")")
  suppressWarnings(data <- lapply(data, log10, how="replace"))
  ## trace <- -7   # log scale
  ## data <- lapply(data, function(x){x[x>trace]}, how="replace")
  ## # Not clipping trace because we want proper min value in metrics
}
#data <- lapply(data, function(x){x[is.finite(x)]})


## beanplots are for a given month pooled by decade

month  <- lapply(ptime, function(x){month.abb[x$mon+1]})
decade <- lapply(ptime, function(x){ 10 * x$year %/% 10 + 1900 })
ocf    <- mapply(rep, names(time), lapply(time, length))

intermediate <- renest(namelist(ocf, decade, month, data))

df <- do.call(rbind, lapply(intermediate, data.frame, stringsAsFactors=FALSE))

df$ocf <- factor(df$ocf, levels=c("obs","cur","fut"))



## color palette for bias-correction
cmap <- c(obs="black", cur="blue", fut="red")

## bean element colors
bcol <- lapply(cmap, function(x){c(adjustcolor(x,0.5),x,x,x)})

## bean elements to plot
what <- c(FALSE, TRUE, TRUE, FALSE)

yr <- range(unlist(data))
xr <- range(unlist(decade))

## plotting

png(outfile, units="in", res=100, width=15, height=11)

par(mfrow=c(4,3), oma=c(0,0,3,0), mar=c(3,5,3,2)+0.1)

for(m in month.abb){
    plot(xr, yr, type="n", main=m, xlab="decade", ylab=data$obs@units)
    axis(side=4, labels=FALSE)
    for(i in rev(levels(df$ocf))){
        sdf <- df[df$month==m & df$ocf == i, ]
        dex <- sort(unique(sdf$decade))
        beanplot(data ~ decade, data=sdf, at=dex, what=what, col=bcol[[i]],
                 add=TRUE, axes=FALSE, maxwidth=8)
    }    
}

mtext(label, line=1, outer=TRUE)

dev.off



## metrics

## No metrics as yet

## Note: use original (non-log) for precip


# summarize <- function(x){c(mean=mean(x, na.rm=TRUE),
#                            stdev=sd(x, na.rm=TRUE),
#                            skew=skewness(x, na.rm=TRUE),
#                            min=min(x, na.rm=TRUE),
#                            max=max(x, na.rm=TRUE))}
# 
# 
# p <- c(10, 25, 50, 75, 90)/100
# pname <- paste0("q", p*100)
# 
# m1 <- lapply(rapply(cdata, summarize, how="replace"), as.matrix)
# m1 <- lapply(m1, t)
# m2 <- lapply(rapply(ddata, quantile, probs=p, how="replace"), as.matrix)
# m2 <- lapply(m2, function(x){rownames(x)<-pname; t(x)})
# 
# 
# ## Empty metrics dataframe
# 
# metrics <- data.frame(infile = "dummy", period="mon", analysis="bean",
#                       mean=0, stdev=0, skew=0, min=0, max=0,
#                       q10=0, q25=0, q50=0, q75=0, q90=0,
#                       stringsAsFactors=FALSE)
# 
# 
# for(p in names(data)){
#   metrics <- rbind(metrics, data.frame(infile=infiles[p],
#                                        period=month.abb,
#                                        analysis="bean",
#                                        m1[[p]],
#                                        m2[[p]],
#                                        row.names=NULL))
# }
# 
# 
# 
# ## write out metrics
# 
# metrics <- metrics[-1,]
# 
# write.table(format(metrics, trim=TRUE, digits=3),
#             file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
# 
