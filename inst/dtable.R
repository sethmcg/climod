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

#args <- c("pfit",
#          "~/Desktop/fig/raw2/all/metrics",
#          "test.pfit")
args <- c("mpdf",
          "lowbar/mpdf/metrics",
          "upper-atm.SGP.mpdf",
          "variable.level.scenario.GCM.RCM.drop.years.drop.lon.lat.drop",
          FALSE,
          TRUE)

## Comment out this line for testing
#args <- commandArgs(trailingOnly=TRUE)

analysis <- args[1]
indir    <- args[2]
outname  <- args[3]


## pattern for spitting filename to columns.  Components are used as
## column names, 'drop' columns are removed.

if(length(args)<4){
    defaultpattern <- TRUE
    colpattern <- "variable.scenario.GCM.RCM.location.drop"
} else {
    defaultpattern <- FALSE
    colpattern <- args[4]
}


## If true, then obs filenames have a dataset component instead of
## RCM/GCM analogs; RCM & GCM are merged to dataset to match.

if(length(args)<5){
    merge <- TRUE
} else {
    merge <- args[5]
}


## Read in all the analysis.*.txt files in indir and bind them
## together into a big dataframe

infiles <- list.files(indir, paste0(analysis,".*.txt"))

rawdata <- lapply(paste0(indir,"/",infiles), read.table,
                  header=TRUE, stringsAsFactors=FALSE)


## Would be nice to make these links instead of just filepaths, but I
## don't know how to get raw HTML to pass through all the DT
## processing.

## Add filenames of analysis pngs.  Need to do this there while
## length(links) matches length(rawdata)

pngs <- paste0(dirname(indir), "/", sub("txt","png",infiles))

xdata <- mapply(function(x,y){cbind(x,png=y)}, rawdata, pngs, SIMPLIFY=FALSE)
rawdf <- do.call(rbind, xdata)


## The first three columns are filename, period, and analysis.

## Remove duplicate rows based on those columns.  (Typically the obs
## row, which is duplicated across analyses of different models.)
## Have to base redundancy on only the first three (identifier)
## columns, because some analyses (e.g., GEV) don't produce
## numerically identical results from run to run.

ddf <- rawdf[!duplicated(rawdf[,1:3]),]

## Reform the dataframe: drop the analysis column and convert the
## filename column to categorical variables.


## extract filename components
fields <- strsplit(basename(ddf$infile), ".", fixed=TRUE)

cnames <- unlist(strsplit(basename(colpattern), ".", fixed=TRUE))

keep <- !grepl("drop", cnames)

categoricals <- as.data.frame(t(as.matrix(fields)),stringsAsFactors=TRUE)[,keep]
colnames(categoricals) <- cnames[keep]

## relevel scenario to put obs first in sorting
categoricals$scenario <- relevel(categoricals$scenario, ref="obs")


# ### OLD:
# 
# ## Model filenames look like: dir/var.scen.GCM.RCM.loc.nc
# ## Obs data filenames look like: dir/var.obs.dataset.loc.nc
# 
# fields <- strsplit(basename(ddf$infile), ".", fixed=TRUE)
# 
# ## Obs: drop last field (.nc)
# obsfields <- lapply(fields, `[`, -5)
# 
# ## Model: drop last field, combine RCM & GCM into dataset
# modfields <- lapply(fields, function(x){c(x[1], x[2], paste0(x[4],'.',x[3]), x[5])})
# 
# 
# ## Use grepl instead of grep (returns logical vector instead of
# ## indices) her because there may be no obs at all, in which case grep
# ## would return 0, and array[-0] = array[0] = nothing.
# obsmask <- grepl("obs", ddf$infile)
# modmask <- !obsmask
# 
# catlist <- list()
# catlist[obsmask] <- obsfields[obsmask]
# catlist[modmask] <- modfields[modmask]
# 
# 
# ## Convert into a dataframe
# categoricals <- as.data.frame(t(as.matrix(catlist)))
# colnames(categoricals) <- c("variable", "scenario", "dataset", "location")

## Add column with links to analyses
links <- ddf$infile


dframe <- cbind(categoricals, ddf[,c(-1,-3)])


## if period is monthly, change column name and sort correctly
if(all(dframe$period %in% month.abb)){
    colnames(dframe)[which(names(dframe) == "period")] <- "month"
    dframe$month <- factor(dframe$month, levels=month.abb)
}

## drop nodata data points (used in test suite)
if("location" %in% colnames(dframe)){
    dframe <- dframe[!dframe$location == "nodata",]
}

## Drop obs for analyses where all metrics are relative to obs
if(analysis %in% c("mpdf", "pfit", "tmmd")){
    dframe <- dframe[dframe$scenario != "obs",]
}

## drop any factor levels that are now unused
dframe <- droplevels(dframe)

# ## initial sort by categoricals
# #sdf <- dframe[with(dframe, order(variable, scenario, location, dataset)), ]
# sdf <- dframe[with(dframe, order(c(colnames(categoricals),"period"))),]


#for (v in c("variable", "scenario", "dataset", "location", "period")){
#    sdf[[v]] <- factor(sdf[[v]])
#}

##############################
## Datatable styling functions

## Note: colors need to be either hex (which RColorBrewer and
## colorspace emit) or names for HTML; R color names don't work.

##############################
## Give categorical columns background colors.

## May need to pull from different parts of multiple palettes in order
## to get enough distinct colors, plus we want pastel colors that are
## easy to read black text on top of, so use this function to get
## select chunks of colorbrewer palettes.  Also wraps the HCL palettes
## to have the same signature so that the catcolor function (which
## applies the styling) can pass arguments through.

## color palette functions for categorical data
pal <- function(N, pname, just=c("just","left","right","center"), reverse=FALSE){
    f <- ifelse(reverse, rev, I)
    if(pname %in% c("rainbow_hcl","heat_hcl","terrain_hcl","diverging_hcl")){
        f(do.call(pname, list(N)))
    } else {
        just <- match.arg(just)
        ## justified colormap
        if(just == "just"){
            return(f(brewer.pal(N, pname)))
        }
        np <- brewer.pal.info[pname,"maxcolors"]
        p <- brewer.pal(np, pname)
        ## left-justified colormap
        if(just == "left"){            
            return(f(p[1:N]))
        }
        ## right-justified colormap
        if(just == "right"){
            return(f(p[np-1:N]))
        }
        ## centered colormap
        ## (drop mid value if N odd & np even or vice-versa)
        if(just == "center"){
            x <- floor((np-N)/2)
            p <- p[(x+1):(np-x)]
            while(length(p) > N){
                p <- p[-(length(p)+1)/2]
            }
            return(f(p))
        }
    }
}


## Function to apply background color styling to categorical columns
catcolor <- function(cname, pname, ...){
    if(cname %in% colnames(sdf)){
        html <- html %>%
            formatStyle(cname,
                        backgroundColor=styleEqual(
                            levels(sdf[[cname]]),
                            pal(nlevels(sdf[[cname]]), pname, ...)
                            )
                        )        
    }
    return(html)
}

##############################
## Numerical columns get a few different styles, depending on
## characteristics of the metric.


##########
## Color palette functions for numerical data

## Background color: red-to-blue with white in the middle
rwbpal <- c(hsv(1, (10:1)/10, 1),
            "#FFFFFF",
            hsv(2/3, (1:10)/10, 1)
            )
## for -1:1 correlations
redbluecorr <- styleInterval(seq(-0.95,0.95,0.1), rwbpal)


## White to green
wgpal <- c(hsv(0.35, (0:10)/12, 1-(0:10)/40))
## for 0:1 skill scores
whitegreenskill <- styleInterval(seq(0.05, 0.95, 0.1), wgpal)
## and and 0:100 percentages
whitegreenpct <- styleInterval(seq(5, 95, 10), wgpal)

## constant-intensity rainbow for days of the year
rainbowdoy <- styleInterval(1:365, rainbow_hcl(366, c=50, l=100))



## Create and style the datatable

## Note that different variables need to go into different datatables,
## because otherwise we can't style the numeric columns because their
## ranges don't match.

## (Figure out how to deal with the normal case of no 1 var / no levels later)

for (V in levels(dframe$variable)){
    for (L in levels(dframe$level)){

        v <- paste0(V, sub("p","",L))

        outfile <- paste0(outname, ".", v, ".html")
        outdir  <- paste0(outname, ".", v, "_files")
        caption <- paste0(gsub(".", " ", basename(outname), fixed=TRUE), " ", V, " ", L)

        sdf <- dframe[dframe$variable == V & dframe$level == L,]


        ## Background-color styling for categorical columns

        html <- datatable(sdf,
                          options = list(paging=FALSE),
                          rownames=FALSE,
                          caption=caption,
                          filter="top"
                          )

        ## scenario has special custom coloring
        
        if("scenario" %in% colnames(sdf)){
            html <- html %>%
                formatStyle("scenario",
                            backgroundColor=styleEqual(
                                c("obs","hist","rcp45","rcp85"),
                                c("#FFFFFF","#E4E4FF","#FFFFE4","#FFE4E4")
                                ))
        }


        html <- catcolor("dataset",  "Pastel1")
        html <- catcolor("location", "rainbow_hcl")
        html <- catcolor("month", "rainbow_hcl")
        html <- catcolor("RCM",      "Pastel1",  just="right")
        html <- catcolor("GCM",      "Pastel1",  just="left")
        html <- catcolor("lat",      "RdBu",     just="center")
        html <- catcolor("lon",      "Spectral", just="center", reverse=TRUE)
        html <- catcolor("variable", "Pastel2")
        html <- catcolor("level",    "Greys",    just="left")
        

        ## Background-color or background-bar styling for numeric columns
        
        ## get range of numeric columns in data table
        vrange <- lapply(sdf[,sapply(sdf, is.numeric)], range, na.rm=TRUE)

        
        ## Gray background bars from min to max
        barmetrics <- c("rlevlo", "rlevmid", "rlevhi",       # gev
                        "dmaxval", "dminval",                # tmmd
                        "mean", "skew", "min","max",         # mpdf
                        "q10", "q25", "q50", "q75", "q90",   # mpdf
                        "pislope",                           # intsp
                        "ktau", "q50over", "modeP", "modeT"  # joint
                        )
        
        ## Gray background bars from 0 to max
        devmetrics <- c("freqmad", "intmad", "totmad",       # pfit
                        "dmad",                              # tmmd
                        "stdev"                              # mpdf
                        )

        ## Set min of range for deviation metrics to 0
        for (m in devmetrics){
            if(!is.null(vrange[[m]])) { vrange[[m]][1] <-0 }
        }

        ## apply background bar style
        for(m in c(barmetrics, devmetrics)){
            if(m %in% colnames(sdf)){
                html <- html %>% formatStyle(m, background=styleColorBar(vrange[[m]], "#CCCCCC"))
            }
        }


        ## background color styled metrics

        bgmetrics <- list(
            freqcor = redbluecorr,
            intcor  = redbluecorr,
            totcor  = redbluecorr,
            pdfsk   = whitegreenskill,
            tailsk  = whitegreenskill,
            dpctpos = whitegreenpct,
            dminday = rainbowdoy,
            dmaxday = rainbowdoy
            )

        for(m in names(bgmetrics)){
            if(m %in% colnames(sdf)){
                html <- html %>% formatStyle(m, backgroundColor=bgmetrics[[m]])
            }
        }

        ## NOTE!  NEED TO TEST gev / pfit / tmmd!
        ## ONLY mpdf and bean are currently tested!

        ## Save to file

        htmltools::save_html(html, file=outfile, lib=outdir, background="white")
    }
}
