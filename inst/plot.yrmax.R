#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")
#library(ncdf4)

## Call as: Rscript --vanilla adj.R label obs cur fut out varname slice

args <- commandArgs(trailingOnly=TRUE)

#args <- c("rcp85 HadGEM2-ES WRF ftlogan",
#          "plot-test/dmaps/tmax.rcp85.HadGEM2-ES.WRF.ftlogan.Rdata",
#          "fig/adj.tmax.rcp85.HadGEM2-ES.WRF.ftlogan.png",
#          "tmax")

label <- args[1]

load(args[2])

outfile <- args[3]

v <- args[4]


png(outfile, units="in", res=120, width=7, height=7)

data <- lapply(annmax, `[[`, v)
rng <- range(unlist(data), na.rm=TRUE)

xlab <- paste("raw", v, paste0("(",units,")"))
ylab <- paste("BC",  v, paste0("(",units,")"))

plot(0,0, xlim=rng, ylim=rng, type="n", xlab=xlab, ylab=ylab, main=label)

abline(0,1)
points(data$rawfut, data$fixfut, col="red")
points(data$rawcur, data$fixcur, col="blue")
points(data$obs, data$obs)

legend("bottomright", c("obs","cur","fut"),
       col=c("black","blue","red"), pch=1)
       
dev.off()


