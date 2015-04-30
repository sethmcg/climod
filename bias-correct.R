library(ncdf4)
library(KernSmooth)
library(devtools)
load_all()

#for(varname in c("prec","tmax","tmin")){
#    for(model in c("gcm","rcm")){
#
#        print("#################")
#        print(paste(varname, model))

varname <- "prec"
model   <- "gcm"

norm <- ifelse(varname == "prec", "power", "zscore")

infile.o <- paste0("test/",varname,".",model,".obs.nc")
infile.c <- paste0("test/",varname,".",model,".cur.nc")
infile.f <- paste0("test/",varname,".",model,".mid.nc")

#>#path <- "../pampas/data/raw/daily/"
#>#ocf <- c("obs","cur","fut","mid")
#>#infile <- as.list(paste0(path, varname, ".", model, ".", ocf, ".nc"))
#>#names(infile) <- ocf        

outfile.c <- paste0("test/",varname,".",model,".cur.bc.nc")
outfile.f <- paste0("test/",varname,".",model,".mid.bc.nc")

file.copy(infile.c, outfile.c, overwrite=TRUE)
file.copy(infile.f, outfile.f, overwrite=TRUE)

fin.o <- nc_open(infile.o)
fin.c <- nc_open(infile.c)
fin.f <- nc_open(infile.f)
#>#fin <- lapply(infile, nc_open)
        
fout.c <- nc_open(outfile.c, write=TRUE)
fout.f <- nc_open(outfile.f, write=TRUE)


## get netcdf data in friendly form
nc.o <- nc_ingest(fin.o)
nc.c <- nc_ingest(fin.c)
nc.f <- nc_ingest(fin.f)
#>#nc <- lapply(fin, nc_ingest)

        
## set up arrays for storing results
save.c <- array(dim=dim(nc.c[[varname]]))
save.f <- array(dim=dim(nc.f[[varname]]))


## generate climatology moving window index arrays
## 360 1-ish-day moving windows w/ 30-day outer pool
cwin.o <- cslice(nc.o$time, outer=30, num=360)
cwin.c <- cslice(nc.c$time, outer=30, num=360)
cwin.f <- cslice(nc.f$time, outer=30, num=360)


#>#time <- lapply(nc, '[[', "time")
#>#cwin <- lapply(time, climwindow, outer=30, num=360)
    
## iterate over space
xname <- names(dimnames(nc.o[[varname]]))[1]
yname <- names(dimnames(nc.o[[varname]]))[2]

nx <- length(nc.o[[xname]])
ny <- length(nc.o[[yname]])
    
for(x in 1:nx){
    for(y in 1:ny){

        print(paste0("x:",x,", y:",y))

        ## get data from nc object 
        data.o <- nc.o[[varname]][x,y,]
        data.c <- nc.c[[varname]][x,y,]
        data.f <- nc.f[[varname]][x,y,]
        #>#data <- lapply(nc,function(d){d[[varname]][x,y,]})

        ## fix up precip
        if(varname == "prec"){
            ## floor all datasets at zero
            data.o <- pmax(data.o, 0)
            data.c <- pmax(data.c, 0)
            data.f <- pmax(data.f, 0)
            #>#data <- lapply(data,pmax,0)

            ## calculate wet/dry fraction
            nt <- length(data.o)
            pwet <- sum(data.o > 0) / nt

            ## find threshold equalizing model wet/dry with obs
            threshold <- sort(data.c)[(1-pwet)*nt]

            ## set values below threshold to zero
            data.c[data.c < threshold] <- 0
            data.f[data.f < threshold] <- 0

            ## set zeros to NA for carrying through calculations
            if(any(is.na(c(data.c, data.f, data.o)))){
                warning("Some values already NA when setting zero -> NA")
            }
            data.o[data.o == 0] <- NA
            data.c[data.c == 0] <- NA
            data.f[data.f == 0] <- NA
        }

        ## if all data for one input dataset is NA, result is NA
        if(all(is.na(data.o)) || all(is.na(data.c)) || all(is.na(data.f))){
            next
        }
        
        ## window data using outer window
        wind.o <- slice(data.o, cwin.o, outer=TRUE)
        wind.c <- slice(data.c, cwin.c, outer=TRUE)
        wind.f <- slice(data.f, cwin.f, outer=TRUE)
        
        
        ## normalize each window separately
        nwd.o   <- lapply(wind.o, normalize, norm=norm)
        nwd.c   <- lapply(wind.c, normalize, norm=norm)
        nwd.f   <- lapply(wind.f, normalize, norm=norm)


        ## construct distribution mappings
        dmaps <- mapply(distmap, nwd.c, nwd.o, SIMPLIFY=FALSE)
        

        ## get normalized inner window for bias correction
        fixme.c <- subslice(nwd.c, cwin.c)
        fixme.f <- subslice(nwd.f, cwin.f)

        ## bias-correct innner window
        fixed.c <- mapply(predict, dmaps, fixme.c)
        fixed.f <- mapply(predict, dmaps, fixme.f)


        ## copy atts from nwd to fixed
        fixed.c <- mapply(copyatts, nwd.c, fixed.c, SIMPLIFY=FALSE)
        fixed.f <- mapply(copyatts, nwd.c, fixed.f, SIMPLIFY=FALSE)
                

        if(varname == "prec"){
          ## check & fix negative values
          neg.c <- sum(unlist(fixed.c) < 0, na.rm=TRUE)
          if(neg.c > 0){
            warning(paste(neg.c,"negative values in current after bias correction"))
            fixed.c <- lapply(fixed.c, pmax, 0)
          }
          neg.f <- sum(unlist(fixed.f) < 0, na.rm=TRUE)
          if(neg.f > 0){
            warning(paste(neg.f,"negative values in future after bias correction"))
            fixed.f <- lapply(fixed.f, pmax, 0)
          }

          
          ## convert NA back to zero
          fixed.c <- lapply(fixed.c, function(x){x[is.na(x)]<-0;return(x)})
          fixed.f <- lapply(fixed.f, function(x){x[is.na(x)]<-0;return(x)})

          
        }

        ## denormalize bias-corrected data
        bc.c <- mapply(denormalize, fixed.c, nwd.o)
        bc.f <- mapply(denormalize, fixed.f, nwd.o, nwd.c)

                
        ## collate BC inner windows back into timeseries
        result.c <- unslice(bc.c, cwin.c)
        result.f <- unslice(bc.f, cwin.f)

        ## save results in array
        save.c[x,y,] <- result.c
        save.f[x,y,] <- result.f
    }
}
    


## write to file --> nc_replace function
## also function for adding history notes

ncvar_put(fout.c, varname, save.c)
ncvar_put(fout.f, varname, save.f)

ncatt_put(fout.c, varname, "bias_correction", "KDDM")
history.c = ncatt_get(fout.c,0,"history")$value
history.c = paste0(date(),
    ": bias-corrected using R package 'climod'\n",
    "    see bias_correction attribute for details\n",
    history.c)
ncatt_put(fout.c,0,"history",history.c)
bcnote = paste(
    "Bias-corrected using KDDM (kernel density distribution mapping)\n",
    "calculated over a 30-day moving window across years\n",
    "observed data from",infile.o,"\n",
    "current data from",infile.c,"\n",
    "future data from",infile.f)
ncatt_put(fout.c,0,"bias_correction",bcnote)

    
        

ncatt_put(fout.f, varname, "bias_correction", "KDDM")
history.f = ncatt_get(fout.f,0,"history")$value
history.f = paste0(date(),
    ": bias-corrected using R package 'climod'\n",
    "    see bias_correction attribute for details\n",
    history.f)
ncatt_put(fout.f,0,"history",history.f)
bcnote = paste(
    "Bias-corrected using kernel density distribution mapping\n",
    "calculated over a 30-day moving window across years\n",
    "observed data from",infile.o,"\n",
    "current data from",infile.c,"\n",
    "future data from",infile.f)
ncatt_put(fout.f,0,"bias_correction",bcnote)
    

nc_close(fout.c)
nc_close(fout.f)

#    }
#}
