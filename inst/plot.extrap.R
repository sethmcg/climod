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


xlab <- paste("raw", v, paste0("(",units,")"))
ylab <- paste("BC",  v, paste0("(",units,")"))

cmap <- c(train="gray", cur="blue", fut="red")

png(outfile, units="in", res=120, width=7, height=7)

par(mfrow=c(2,2), oma=c(0,0,2,0), mar=c(4,4,3,2), mgp=c(2.5,1,0))

for (p in names(peaks)){

  od <- renest(outerdata)[[p]]
  id <- renest(innerdata)[[p]]

  pdata <- list(train = lapply(renest(od)$cur, unlist),
                cur = renest(id)$cur,
                fut = renest(id)$fut)

  mplot(pdata, x="raw", y="fix", xlab=xlab, ylab=ylab,
        pch=c(19,1,1), col=cmap, cex=c(1.5,1,1))

  abline(0,1)
  
  title(paste0(p ,": ", dates[p]), line=1)
}

mtext(label, outer=TRUE)


## centered legend outside plots

par(fig=c(0,1,0,1), new=TRUE)
plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)
legend("center", c("train","cur","fut"),
        pch=c(19,1,1), col=cmap, pt.cex=c(1.5,1,1))


dev.off()


