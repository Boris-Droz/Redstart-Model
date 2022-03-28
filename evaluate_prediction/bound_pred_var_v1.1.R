####################################################################################################################################
###					
###  Summary of boundary predicting variables values for each project zone
### ==================================================================###########
###
####################################################################################################################################
####################################################################################################################################
## HISTORIC:
## --------
##  v1.1, February 2021 - Boris DROZ Grafb 
####################################################################################################################################
## DESCRIPTION
###############
# Input: 
#  - folders with divers raster used as current predicting variables
#  - extent (zone) used for the calibration of the model
## - do in batch for diff project zone
####################################################################################################################################
## Load library
##################
## Set library path
# .libPaths("C:/Users/Public/Documents/R/win-library/3.1"
# stat analysis
library(raster)
library(stringr)
library(maptools)

# Set working_dir which contain 
# input folder contening the different folder project with predicting variable raster
# for each of them ...
# like ~/Cdf/pred/current
## =========#############
working_dir <-"~/Actual/rafb_proj_model/analysis/" 

# creat manualy an output folder
outpath <- paste(working_dir,"/Cdf_2021/eval_output/new.rang_eval", sep="")

# see print(list.proj --> line 51 on the code) to select the project
sel.mod <- c(1,3,4,5,7,8,10,11,12)

###################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################
##############################################
setwd(working_dir)

# load all zone of mod
list.proj <- list.dirs(path=working_dir, recursive = FALSE)
list.name <- list.dirs(path=working_dir, full.names = FALSE, recursive = FALSE)
print(list.name)

#  built a data table
out.table <-NULL

for (k in sel.mod) # put 2:length(fns) in case of arcgis version (see above)
  {
  
  inpath <- paste(list.proj[k],"/pred/pred_cur", sep="")

    # pred. var file (batch all file in input) 
    # fns <- list.dirs(inpath,full.names = TRUE); #print(fns) ## for arcgis file in folder
    fns <- list.files(inpath,pattern=".tif$",full.names = TRUE)

	png(paste(outpath,"/",list.name[k],"_var.pred.png",sep=""),width = 560, height = 560) 
    
	plot(stack(fns))
	dev.off()
    
    cali.area <- rgdal::readOGR(paste(list.proj[k],"/zone/zone_mod_LV95.shp",sep=""))
    
    out.dat <- NULL
    out.name <-NULL
    for (n in 1:length(fns)) # put 2:length(fns) in case of arcgis version (see above)
        {
          p.dat <- raster(fns[n]); #plot(p.dat)
          
          #p.dat <- mask(p.dat,cali.area)
          
          r.min <- cellStats(p.dat, stat='min')
          r.max <- cellStats(p.dat, stat='max')
          r.mean <- cellStats(p.dat, stat='mean')
          r.sd <- cellStats(p.dat, stat='sd')
          
          out.dat <-cbind(out.dat,c(min=r.min,mean =r.mean,max=r.max,sd=r.sd))
          out.name <-c(out.name,names(p.dat))
        }
    out.name -> colnames(out.dat)
    
    out.table <- rbind(out.table,data.frame(zone=rep(list.name[k], 4),out.dat))
  }    

# write summary of the data coverage
write.table(out.table, file=paste(outpath,"/boundaries_pred.var_zone.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

