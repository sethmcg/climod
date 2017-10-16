#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
library(ncdf4)
suppressMessages(library(extRemes))


## Call as: Rscript plot.gev.R label obs cur fut png var txt

## png = name of output file for plotted figure
## txt = name of output file for plain-text table of metrics

args <- c("rcp85 HadGEM2-ES RegCM4 birmingham",
          "obs/prec.obs.livneh.birmingham.nc",
          "save-test/prec.hist.HadGEM2-ES.RegCM4.birmingham.nc",
          "save-test/prec.rcp85.HadGEM2-ES.RegCM4.birmingham.nc",
          "test.gev.png",
          "prec",
          "test.gev.txt"
          )

## Comment out this line for testing
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


time <- lapply(nc,"[[","time")

time <- lapply(time, alignepochs, "days since 1950-01-01")

### find block maxima
year <- lapply(time, function(x){floor(x/yearlength(x))})
ydata <- mapply(split, data, year)
bmax <- lapply(ydata, function(x){sapply(x, max, USE.NAMES=FALSE)})

yyear <- lapply(year, function(x){unique(x)+1950})
annmax <- lapply(mapply(cbind, yyear, bmax), as.data.frame)
annmax <- lapply(annmax, `colnames<-`, c("year","prec"))

fits <- lapply(bmax, fevd, units=units)


nsfits <- lapply(annmax, function(x){
  fevd(data=x, x=prec, annmax, units=units,
       location.fun=~year, scale.fun=~year)})


nsper <- c(50,100,200)
nserlev <- lapply(lapply(nsfits, erlevd, period=nsper), t)


png(outfile, units="in", res=120, width=7, height=7)

ylim <- c(0, max(sapply(fits, function(x){ci(x, return.period=500)[3]})))

par(mfrow=c(2,2), oma=c(0,0,3,0), mgp=c(2,1,0), mar=c(4,3.5,2.5,1))

for(f in names(fits)){
  plot(fits[[f]], type="rl", main=f, col=cmap[f], ylim=ylim)
  abline(h=pretty(ylim), v=c(2,5,10,20,50,100,200,500), col="gray", lwd=0.5)
}

mplot(annmax, x="year", y="prec", type="l", col=cmap, lty=1,
      xlab="year", ylab=units, main="Non-st'nary ERLs 50,100,200")


for(i in names(yyear)){
  matplot(yyear[[i]], nserlev[[i]], add=TRUE, type="l", lty=3, col=cmap[i])
}

mtext(label, line=1, outer=TRUE)
    
dev.off()



## Calculate

metrics <- data.frame(infile = "dummy", period="Xyr", analysis="gev",
                      rlevlo=0, rlevmid=0, rlevhi=0,
                      
                      stringsAsFactors=FALSE)

rps <- c(5, 10, 20, 50, 100)

for(p in names(fits)){

  m1 <- data.frame(infile=infiles[p], 
                   period=paste0(rps, "yr"),
                   analysis="gev",
                   row.names=NULL)
  
  m2 <- as.data.frame(unclass(ci(fits[[p]], return.period=rps)))
  rownames(m2) <- NULL
  colnames(m2) <- c("rlevlo", "rlevmid", "rlevhi")

  metrics <- rbind(metrics, cbind(m1, m2))
}



## write out metrics

metrics <- metrics[-1,]

write.table(format(metrics, trim=TRUE, digits=3),
            file=txtfile, quote=FALSE, sep="\t", row.names=FALSE)
