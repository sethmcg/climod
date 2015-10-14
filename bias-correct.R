library(ncdf4)
library(KernSmooth)
library(devtools)
load_all()

## width of moving window
mwinwidth = 30

ocf <- c("obs","cur","fut")

for(varname in c("prec","tmax","tmin")){
    
    norm <- ifelse(varname == "prec", "power", "zscore")
    
    print("#################")
    print(paste(varname))
    
    
    indir <- "tests/raw"
    infiles <- paste0(indir, "/", varname, ".", ocf, ".nc")
    names(infiles) <- ocf
            
    
    nc <- lapply(infiles, nc_ingest)
    
    
    ## get netcdf data in friendly form
    
    time <- lapply(nc,"[[","time")
    
    indata <- lapply(nc,"[[",varname)
    
    
    ### set up storage data structures
    
    outdata <- lapply(indata,"+",NA)
    
    
    ### generate climatology moving window index arrays
    ### 360 1-ish-day moving windows w/ 30-day outer pool
    cwin <- lapply(time, cslice, outer=mwinwidth, num=360)
        
    
    ### Assert all 3 datasets have same spatial dimensions and are ordered
    ### lon, lat, time...
    
    nxy <- dim(indata$obs)[1:2]
    nx <- nxy[1]
    ny <- nxy[2]
    
    
    for(x in 1:nx){
        for(y in 1:ny){
    
            print(paste0("x:",x,", y:",y))
    
            ## extract data for this gridcell
            data <- lapply(indata, function(a){a[x,y,]})
                        
    
            ##  For stability, it's best to dedrizzle precip data here        
            if(varname == "prec"){

                data <- dedrizzle(data)
                
                ### MOVE TO UNZERO FUNCTION

                ## Make this warning switchable (default off)
                if(any(is.na(unlist(data)))){
                    warning("Some values already NA")
                }
                data <- lapply(data, function(a){a[a == 0] <- NA; a})
    
                ### END UNZERO FUNCTION
    
            }

            
            ## If all data for one input dataset is NA, result is NA.
            ## No warnings or errors are needed here; all NA is
            ## expected over oceans, outside domain, etc.
            
            if(any(sapply(data, function(a){all(is.na(a))}))){
                for(i in names(result)){
                    outdata[[i]][x,y,] <- NA
                }                
                next
            }
            
            
            ## window data using outer window
            wind <- mapply(slice, data, cwin, MoreArgs=list(outer=TRUE), SIMPLIFY=FALSE)
            
            ## normalize each window separately
            nwd <- rapply(wind, normalize, norm=norm, how="replace")
    
            ## construct distribution mappings
            dmaps <- mapply(distmap, nwd$cur, nwd$obs, SIMPLIFY=FALSE)
            
            ## get normalized inner window for bias correction
            fixme <- mapply(subslice, nwd, cwin, SIMPLIFY=FALSE)
    
            ## (drop obs from fixme)
            fixme <- fixme[-(names(fixme)=="obs")]
            
            ## bias-correct innner window
            fixed <- lapply(fixme, function(x){mapply(predict, dmaps, x)})
    
            ## copy over normalization attributes onto fixed
            fixed$cur <- mapply(copyatts, nwd$cur, fixed$cur, SIMPLIFY=FALSE)
            fixed$fut <- mapply(copyatts, nwd$fut, fixed$fut, SIMPLIFY=FALSE)
    
            if(varname == "prec"){
                ## MOVE TO REZERO FUNCTION
                
                ## fixed values may be negative; set to zero.
                ## (Don't warn; this is normal)
                fixed <- rapply(fixed, function(x){pmax(x,0)}, how="replace")
              
                ## convert NA back to zero
                fixed <- rapply(fixed, function(x){x[is.na(x)]<-0;return(x)}, how="replace")
    
                ## END REZERO FUNCTION          
            }
    
            ## denormalize bias-corrected data
            bc <- list()
            bc$obs <- subslice(wind$obs, cwin$obs)
            bc$cur <- mapply(denormalize, fixed$cur, match=nwd$obs)
            bc$fut <- mapply(denormalize, fixed$fut, match=nwd$obs, adjust=nwd$cur)
                    
            ## collate BC inner windows back into timeseries
            result <- mapply(unslice, bc, cwin)
            
            
            ## save results in array
            for(i in names(result)){
                outdata[[i]][x,y,] <- result[[i]]
            }
        }
    }
        
    ###################
    
    
    
    outdir <- "tests"
    outfiles <- paste0(outdir, "/", varname, ".", ocf, ".bc.nc")
    names(outfiles) <- ocf
    
    
    for(i in c("cur","fut")){
    
        ### BEGIN NC_REPLACE FUNCTION
        
        file.copy(infiles[i], outfiles[i], overwrite=TRUE, copy.mode=FALSE)
        
        fout <- nc_open(outfiles[i], write=TRUE)
        
        ncvar_put(fout, varname, outdata[[i]])
    
        ncatt_put(fout, varname, "bias_correction", "KDDM")
    
    
        ## BEGIN APPEND HISTORY FUNCTION
    
        historyatt <- ncatt_get(fout,0,"history")$value
    
        historyatt = paste0(date(),
            ": bias-corrected using R package 'climod'\n",
            "    see bias_correction attribute for details\n",
            historyatt)   
    
        ncatt_put(fout, 0, "history", historyatt)    
    
        ## END APPEND HISTORY FUNCTION
    
        
        bcnote = paste(
            "Bias-corrected using kernel density distribution mapping\n",
            "calculated over a", mwinwidth, "day moving window across years\n",
            "observed data from",infiles["obs"],"\n",
            "current data from",infiles["cur"],"\n",
            "future data from",infiles["fut"])
        
    
        ncatt_put(fout, 0, "bias_correction", bcnote)
        
        nc_close(fout)
    
        ## END NC REPLACE FUNCTION
    }    
}
