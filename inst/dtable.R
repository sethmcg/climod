library(DT)
library(RColorBrewer)
library(colorspace)
library(devtools)
load_all("~/Desktop/climod")

## Call as: Rscript dtable.R analysis indir outname

## analysis = type of analysis (pfit, gev, etc.)
## indir    = directory where metrics files (type.*.txt) are found
## outname  = name of output files (out.html & out_files/)

## Note: if you only want some of the files in the directory, just run
## it for everything and filter them out in the datatable.

## Note also: if outname is html/dirname, saving generates an error
## message and saves to html/html/dirname instead.

args <- c("pfit",
          "~/Desktop/fig/raw2/all/metrics",
          "test.pfit")

## Comment out this line for testing
args <- commandArgs(trailingOnly=TRUE)

analysis <- args[1]
indir    <- args[2]
outname  <- args[3]


outfile <- paste0(outname, ".html")
outdir  <- paste0(outname, "_files")
caption <- gsub(".", " ", basename(outname), fixed=TRUE)


## Read in all the analysis.*.txt files in indir and bind them
## together into a big dataframe

infiles <- list.files(indir, paste0(analysis,".*.txt"))

rawdata <- lapply(paste0(indir,"/",infiles), read.table,
                  header=TRUE, stringsAsFactors=FALSE)

rawdf <- do.call(rbind, rawdata)


## The first three columns are filename, period, and analysis.

## Remove duplicate rows based on those columns.  (Typically the obs
## row, which is duplicated across analyses of different models.)
## Have to base redundancy on only the first three (identifier)
## columns, because some analyses (e.g., GEV) don't produce
## numerically identical results from run to run.

ddf <- rawdf[!duplicated(rawdf[,1:3]),]

## Reform the dataframe: drop the analysis column and convert the
## filename column to categorical variables.

## Model filenames look like: dir/var.scen.GCM.RCM.loc.nc
## Obs data filenames look like: dir/var.obs.dataset.loc.nc

fields <- strsplit(basename(ddf$infile), ".", fixed=TRUE)

## Obs: drop last field (.nc)
obsfields <- lapply(fields, `[`, -5)

## Model: drop last field, combine RCM & GCM into dataset
modfields <- lapply(fields, function(x){c(x[1], x[2], paste0(x[4],'.',x[3]), x[5])})


## Use grepl instead of grep (returns logical vector instead of
## indices) her because there may be no obs at all, in which case grep
## would return 0, and array[-0] = array[0] = nothing.
obsmask <- grepl("obs", ddf$infile)
modmask <- !obsmask

catlist <- list()
catlist[obsmask] <- obsfields[obsmask]
catlist[modmask] <- modfields[modmask]


## Convert into a dataframe
categoricals <- as.data.frame(t(as.matrix(catlist)))
colnames(categoricals) <- c("variable", "scenario", "dataset", "location")


dframe <- cbind(categoricals, ddf[,c(-1,-3)])

## relevel scenario to put obs first in sorting
if("obs" %in% levels(dframe$scenario)){
    dframe$scenario <- relevel(dframe$scenario, ref="obs")
}

## drop nodata data points (used in test suite)
dframe <- dframe[!dframe$location == "nodata",]


## initial sort by variable, scenario, location, and dataset
sdf <- dframe[with(dframe, order(variable, scenario, location, dataset)), ]


## Drop obs for analyses where all metrics are relative to obs
if(analysis == "pfit" || analysis == "tmmd"){
    sdf <- sdf[!sdf$scenario == "obs",]
}



## Create the datatable
    
## Note: colors need to be either hex (which RColorBrewer and
## colorspace emit) or names for HTML; R color names don't work.


## Basic colorization for categoricals first

html <- datatable(sdf,
                  options = list(paging=FALSE),
                  rownames=FALSE,
                  caption=caption,
                  filter="top"
                  ) %>%
    formatStyle("scenario",
                backgroundColor=styleEqual(
                    c("obs","hist","rcp45","rcp85"),
                    c("#FFFFFF","#E4E4FF","#FFFFE4","#FFE4E4")
                    )) %>%
    formatStyle("dataset",
                backgroundColor=styleEqual(
                    levels(sdf$dataset),
                    brewer.pal(nlevels(sdf$dataset), "Pastel1")
                    )) %>%
    formatStyle("location",
                color=styleEqual(
                    levels(sdf$location),
                    rainbow_hcl(nlevels(sdf$location))
                    ))


## Red to blue bg color style for correlations, etc.
redblue <- styleInterval(seq(-0.95,0.95,0.1),
                         c(hsv(1,   (10:1)/10, 1),
                           "#FFFFFF",
                           hsv(2/3, (1:10)/10, 1)
                           ))

## Red to blue bg color style for percentages
pctredblue <- styleInterval(seq(2.5,97.5,5),
                         c(hsv(1,   (10:1)/10, 1),
                           "#FFFFFF",
                           hsv(2/3, (1:10)/10, 1)
                           ))

## Rainbow bg color style for days of the year
doyrainbow <- styleInterval(1:365, rainbow_hcl(366, c=50, l=100))



if(analysis == "gev"){

    ## background bars for the different return levels
    for(m in c("rlevlo", "rlevmid", "rlevhi")){
        html <- html %>% formatStyle(m, background=styleColorBar(sdf[[m]], "#CCCCCC"))
    }
}



if(analysis == "pfit"){

    ## range of MAD metrics for background bars
    vrange <- lapply(sdf[,c("freqmad","intmad","totmad")], range, na.rm=TRUE)
    vrange <- lapply(vrange, function(x){x[1]<-0;x})  ## set min for MAD to zero

    ## red-blue background for correlations
    for(m in c("freqcor", "intcor", "totcor")){
        html <- html %>% formatStyle(m, backgroundColor=redblue)
    }

    ## background bar from 0 to max for MADs
    for(m in c("freqmad","intmad","totmad")){
        html <- html %>% formatStyle(m, background=styleColorBar(vrange[[m]], "#CCCCCC"))
    }
}



if(analysis == "tmmd"){

    vrange <- lapply(sdf[,c("dmad", "dmaxval","dminval")], range)
    vrange$dmad[1] <- 0  ## set min for MAD to zero

    ## dmad, dminval, dmaxval: colorbars from min to max
    for(d in c("dmad", "dminval", "dmaxval")){
        html <- html %>%
            formatStyle(d, background=styleColorBar(vrange[[d]], "#CCCCCC"))
    }
    
    ## dpctpos: red-blue background 0:100
    html <- html %>% formatStyle("dpctpos", backgroundColor=pctredblue)

    ## minday & maxday: colorwheel
    for(d in c("dminday", "dmaxday")){
        html <- html %>% formatStyle(d, backgroundColor=doyrainbow)
    }
}
    


## Save to file

htmltools::save_html(html, file=outfile, lib=outdir, background="white")

