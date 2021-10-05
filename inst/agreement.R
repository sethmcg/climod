library(devtools)
load_all("~/climod")

## Evaluation of how much PDF agreement one should expect based on
## sample size

size <- round(10^seq(2, 5, length=16))

m <- length(size)
N <- 1000

set.seed(222)

p <- list()
t <- list()
for(s in 1:m){
    cat(size[s], "")
    x <- lapply(rep(size[s], N), rnorm)
    y <- lapply(rep(size[s], N), rnorm)
    p[[s]] <- mapply(pdfskill, x,y)
    t[[s]] <- mapply(tailskill, x,y)
}
cat("\n")


par(mfrow=c(2,1))

boxplot(p, xaxt="n", main=paste("PDF skill"),
        xlab="sample size", ylab="skill score" )
axis(side=1, at=seq(1,m,length=4), labels=10^(2:5))


boxplot(t, xaxt="n", main=paste("tail skill"), ylim=c(0,1),
        xlab="sample size", ylab="skill score" )
axis(side=1, at=seq(1,m,length=4), labels=10^(2:5))

