###################################################
## Merge and compare mns and mnt for further calculation
## work under R 4.0.3
###################################################
# v2.3 - Boris Droz - March 2012

###################################################################################################
## Descriptif
##############
# 
# put all input data into two folder labeled "mnt" and a "mns" 
# in the input folder of your projet
# set the parameter and run the script.....
##
## Reference 
#############
##
##
#####################################################################################################################################
#####################################################################################################################################
##### PARAMETER #################
####################################################################

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
n.proj <- "LV95" # name projection
resol <- 2 # resolution needed if smaller just average it
f.type <- ".tif$" # file extension of the mnt and mns
folder <- "YES" # YES -> several commune are merged = several folder
                # NO -> one comune = one folder
working_dir <- "~/Actual/rafb_proj_model/analysis/JuraVD_2021/"

area.mode <- "zone_mod_LV95.shp" # area used for the prediction
w.area <- "CSEAU_VD.shp" # water area (optionel) if "NULL" not use!

#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################
# find out working path
########################
setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")
source("script/function/geo2file_geosv1.R")

check.lib (c("raster","sp","magrittr","rgdal","rgeos",
             "rgrass7","RNetCDF","splitstackshape","maptools","XML"))
use_sp()
#####################################################################################################################################
# Folder and file
#################
inpath <- paste(working_dir,"/input/",sep="")
creat.subDir(working_dir,"pred/pred_cur")
outpath <- paste(working_dir,"/pred/pred_cur",sep="")
setwd(inpath)
#####################################################################################################################################
############################################################################
beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

area.proj <- crs(paste("+init=epsg:",proj,sep="")) # proection of the project

# project area
###############
# zone of interest
shp.zone <-readShapeSpatial(paste(working_dir,"/zone/",area.mode,sep=""),delete_null_obj=TRUE)
proj4string(shp.zone) <- area.proj

## 1) merge mns
################

if (folder=="YES"){
   fol <- list.dirs(paste(inpath,"/mns/",sep=""))
  
   fns <-NULL
  for (i in 2:length(fol)){
      fns <- c(fns,list.files(fol[i],pattern=f.type,full.names = TRUE) )
      }
  
  }else{
  # Creat list of raster file and open it
  fns <- list.files(paste(inpath,"/mns/",sep="")
                    ,pattern=f.type,full.names = TRUE)
  }

## merge all demain the zone of interest
p<-0
for (c  in 1:length(fns))
 {
  cat("> File",c, "is merged ...", "\n",append = FALSE)
  
    ra <- raster(fns[c])
  
    if ( extent(shp.zone)[1] <= extent(ra)[2] &
         extent(shp.zone)[2] >= extent(ra)[1] &
         extent(shp.zone)[3] <= extent(ra)[4] &
         extent(shp.zone)[4] >= extent(ra)[3] )
        {
         p<- p+1
          if (p==1) {
                dem <- ra 
          }else{  dem <- merge(dem,ra)}    
    }else{}
  }  

# check resolution
if (res(dem)[1] == resol) 
  { 
    cat("> Resolution okay ...", "\n",append = FALSE)
    }else{ if ( resol/res(dem)[1] > 1) 
    { # AGGREGATE THE DATA IF RESOLUTION IS SMALLER
      cat("> Data aggregate ...", "\n",append = FALSE)
      dem <- aggregate(dem,fact=resol/res(dem)[1],mean)
      
    }else{ # RESAMPLING IF RESOLUTION IS BIGGER OR SLIGHLY DIFFERENT
      cat("> !! CAUTION Data is resample...", "\n",append = FALSE)
      empty.raster <- raster(res=c(resol,resol),
                             crs=area.proj,ext=extent(proj.area))
      dem <- resample(dem,empty.raster,method="ngb")
    } 
  }
    
writeRaster(dem, filename = paste(getwd(),"/mns_",n.proj,"_2m.tif",sep=""),
              datatype="FLT8S", overwrite=TRUE)  
   
crs(dem) <- area.proj

dem <- crop(dem,shp.zone) # used only selected area
dem <- raster::mask(dem,shp.zone)

# rewrite mns
writeRaster(dem, filename = paste(getwd(),"/mns_",n.proj,"_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)
###############################################################################
## 2) merge mnt
################
# Creat list of raster file and open it
if (folder=="YES"){
  fol <- list.dirs(paste(inpath,"/mnt/",sep=""))
  
  fns <-NULL
  for (i in 2:length(fol)){
    fns <- c(fns,list.files(fol[i],pattern=f.type,full.names = TRUE) )
  }
  
}else{
  # Creat list of raster file and open it
  fns <- list.files(paste(inpath,"/mnt/",sep="")
                    ,pattern=f.type,full.names = TRUE)
}

## merge all demin the zone of interest
p<-0
for (c  in 1:length(fns))
  {
  cat("> File",c, "is merged ...", "\n",append = FALSE)
  
    ra <- raster(fns[c])
    
    if ( extent(shp.zone)[1] <= extent(ra)[2] &
         extent(shp.zone)[2] >= extent(ra)[1] &
         extent(shp.zone)[3] <= extent(ra)[4] &
         extent(shp.zone)[4] >= extent(ra)[3] )
    {
      p<- p+1
      if (p==1) {
        dem <- ra 
      }else{  dem <- merge(dem,ra)}    
    }else{}
  }  

# check resolution
if (res(dem)[1] == resol) 
  { 
    cat("> Resolution okay ...", "\n",append = FALSE)
  }else{ if ( resol/res(dem)[1] > 1) 
          { # AGGREGATE THE DATA IF RESOLUTION IS SMALLER
          cat("> Data aggregate ...", "\n",append = FALSE)
          dem <- aggregate(dem,fact=resol/res(dem)[1],mean)
      
        }else{ # RESAMPLING IF RESOLUTION IS BIGGER OR SLIGHLY DIFFERENT
            cat("> !! CAUTION Data is resample...", "\n",append = FALSE)
            empty.raster <- raster(res=c(resol,resol),
                                   crs=area.proj,ext=extent(proj.area))
            dem <- resample(dem,empty.raster,method="ngb")
    } 
  }

# write mnt
writeRaster(dem, filename = paste(working_dir,"/input/mnt_",n.proj,"_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# creat buffer area for computing all predictor
b.proj.area <-  buffer(shp.zone,width=200)

# rasterrzie and write the shape as project area
r.empty <- raster(res=res(dem),crs=area.proj,ext=extent(b.proj.area))
b.proj.area <- rasterize(b.proj.area,r.empty)

if(w.area=="NULL"){ }else{
  # remove water area to remove border effect
  shp.water <- readShapeSpatial(paste(inpath,w.area,sep=""),delete_null_obj=TRUE)
  proj4string(shp.water) <- area.proj
  r.water <- rasterize(shp.water,r.empty)
  f.water <- reclassify(r.water, c(NA,NA,1))
  f.water <- reclassify(f.water, c(1,Inf,NA))
  # creat water filter
  b.proj.area <- f.water * b.proj.area
  }

names(b.proj.area) <- "buff.proj"
crs(b.proj.area) <- area.proj

writeRaster(b.proj.area, filename = paste(working_dir,"/zone/buff.proj_",n.proj,"_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

crs(dem) <- crs(shp.zone)

dem <- crop(dem,shp.zone) # used only selected area
dem <- raster::mask(dem,shp.zone)

# rewrite mnt
writeRaster(dem, filename = paste(working_dir,"/input/mnt_",n.proj,"_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

########################################################################
## 3) calculate dem 
###################
mns <- raster(paste(getwd(),"/mns_",n.proj,"_2m.tif",sep="")) 
mnt <- raster(paste(getwd(),"/mnt_",n.proj,"_2m.tif",sep="") )

## Cropping all rasters to the lowest extent #####
ex <- intersect(extent(mns), extent(mnt) )
mns <- crop(mns,ex)
mnt <- crop(mnt,ex)

## if raster still do not have the same exact extent
empty.raster <- raster(res=res(mnt),ext=extent(mnt), crs= crs(mnt )) 

if ( extent(mns) == extent(mnt) ) 
    { }else{
      mns <- resample(mns,empty.raster,method="ngb")
    }

dem <- mns -mnt

crs(dem) <- crs(shp.zone)
  
writeRaster(dem, filename = paste(working_dir,"/input/dem_",n.proj,"_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# creat area zone of projection
zone <- reclassify(dem,c(-Inf,Inf,1))

names(zone) <- "area.proj_LV95_2m"

writeRaster(zone, filename = paste(working_dir,"/zone/zone.mod_",n.proj,"_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(dem);plot(dem)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

#####################################################################################################################################
################  END SCRIPT ###############################################
###########################################################################



