#####
#### Generation mfroadw
####=================================================
####
## work under R 4.0.3
## v1.1 B.Droz February 2021

###################################################################################################
## Descriptif
##############
# ##  creat th variable mfroadw from traffic count

## First, check classification attribute and shp field in the shp input file
## Then set the parameter and run the script
## 
#####################################################################################################################################
##### PARAMETER #################
####################################################################

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
n.proj <- "LV95" # name projection
working_dir <-"~/Actual/rafb_proj_model/analysis/JuraVD_2021/" 
#select file
f.traffic <- "MN95_DGMR_TPR_TJM.shp" # file 
h.traffic <- "TJM_15" # selected field of road
  
#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################}
setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","rgdal","maptools", "gstat", "grainchanger" ))

use_sp()

# set file location 
inpath <- paste(working_dir,"/input",sep="")
creat.subDir(working_dir,"/pred/pred_cur")
outpath <- paste(working_dir,"/pred/pred_cur/",sep="")
setwd(inpath)
#############################################################################
beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#####################################################################################################################################
## OPEN Data
################
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )

road <- rgdal::readOGR(f.traffic)
crs(road) <- crs(paste("+init=epsg:",proj,sep=""))
#road <- spTransform(road,crs(b.proj.area))

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(road, byid=TRUE)==FALSE) >0) 
  { road <- gBuffer(road, byid=TRUE, width=0)}else{ }

# crop to zone of pred
road <- crop(road,b.proj.area)

# rasterized and replot as point
r.empty <- raster(res=res(b.proj.area), ext=extent(b.proj.area), 
                  crs= crs(road) ) 
# select data
road <- road[,names(road)==h.traffic]
names(road) <- "traffic"

road$traffic <- as.numeric(road$traffic)

r.road <- rasterize(road,r.empty,field='traffic')

p.road <- rasterToPoints(r.road, spatial=TRUE)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
gs <- gstat(formula=layer~1, locations=p.road, set=list(idp = 2))
r.road <- interpolate(r.empty, gs)
# r.road <- raster::(r.road,b.proj.area)

surf <- raster(paste(working_dir,"zone/surf_terr.tif",sep=""))
w.road <- winmove(r.road, d=100, type = c("circle"), win_fun=sum)/surf

w.road <- crop(w.road,proj.area)

## if raster still do not have the same exact extent
if ( extent(w.road) == extent(proj.area) ) 
    { }else{
      empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
      w.road <- resample(w.road,empty.raster,method="ngb")
    }

w.road <- mask(w.road,proj.area)

names(w.road) <- "mfroadw"

writeRaster(w.road, filename = paste(outpath,"mfroadw.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)
# --- check --
print(w.road);plot(w.road)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

#####################################################################################################################################
################  END SCRIPT ###############################################
###########################################################################

