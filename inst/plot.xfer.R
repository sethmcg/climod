#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")

## Call as: Rscript plot.xfer.R label save png var

args <- c("rcp85 HadGEM2-ES WRF ftlogan",
          "save-test/dmaps/prec.rcp85.HadGEM2-ES.WRF.ftlogan.Rdata",
          "test.yrmax.png",
          "prec")

## Comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)



label <- args[1]

load(args[2])

outfile <- args[3]

v <- args[4]


png(outfile, units="in", res=120, width=7, height=11)

par(mfrow=c(3,2), oma=c(0,0,2,0), mar=c(3,3,2,1))

for (p in names(peaks)){
  plot(dmaps[[p]], xlab="", ylab="")
  title(paste0(p ,": ", dates[p]), line=1)
}

mtext(label, outer=TRUE)

dev.off()
