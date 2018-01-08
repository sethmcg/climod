#library(climod)
library(devtools)
#load_all("climod")
load_all("~/climod")

## Call as: Rscript plot.extrap.R label save png var

args <- c("prec rcp85 HadGEM2-ES WRF ftlogan",
          "v1.lof.bc/dmaps/prec.rcp85.HadGEM2-ES.WRF.ftlogan.Rdata",
          "test.extrap.png",
          "prec")

## Comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)



label <- args[1]

load(args[2])

outfile <- args[3]

v <- args[4]


xlab <- paste("raw", v, paste0("(",units,")"))
ylab <- paste("BC",  v, paste0("(",units,")"))

cmap <- c(train="gray", cur="blue", fut="red")

png(outfile, units="in", res=120, width=7, height=11)

par(mfrow=c(3,2), oma=c(1,0,2,0), mar=c(4,4,3,2), mgp=c(2.5,1,0))

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


## centered bottom legend outside plots

par(mfrow=c(1,1), fig=c(0,1,0,1), oma=c(0.5,0,0,0), mar=c(0,0,0,0), new=TRUE)
plot(0, 1, type="n", bty="n", xaxt="n", yaxt="n", ann=FALSE)
legend("bottom", c("train","cur","fut"), horiz=TRUE,
        pch=c(19,1,1), col=cmap, pt.cex=c(1.5,1,1))


dev.off()


