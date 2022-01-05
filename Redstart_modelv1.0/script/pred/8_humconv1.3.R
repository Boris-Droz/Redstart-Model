#####
#### Generation humcon
####=================================================
####
## work under R 4.0.3
## v1.3 B.Droz February 2021

###################################################################################################
## Descriptif
##############
# ##  creat the variable humcon from habitant housing  (STATPOP) and worker (STATENT)
#  count geolocalized data point
## 
## First, check classification attribute and shp field in the shp input file
## Then set the parameter and run the script
## 
#####################################################################################################################################
##### PARAMETER #################
####################################################################

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
n.proj <- "LV95" # name projection
#f.type <- ".asc$" # type of file of the mnt
working_dir <-"~/Actual/rafb_proj_model/analysis/JuraVD_2021/" 
#select file -- work for shp 
f.hab <- "STATPOP_VD_2016.shp" # file
h.hab <-"NB_HABIT" # header of the data

f.work<- "STATENT_VD_2016.shp" #file
h.work <- "emptot"# header of the data
#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################
#ptm <- proc.time()# ignite timer

setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","rgdal","maptools", "rgeos","gstat", "grainchanger" ))

# use_sp()

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

hab <- rgdal::readOGR(f.hab)
work <- rgdal::readOGR(f.work)
  
# set the coordinate system
proj4string(hab) <- CRS(paste("+init=epsg:", proj,sep="") )
proj4string(work) <- CRS(paste("+init=epsg:", proj,sep="") )

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(hab, byid=TRUE)==FALSE) >0) 
  { hab <- gBuffer(hab, byid=TRUE, width=0)}else{ }

if (sum(gIsValid(work, byid=TRUE)==FALSE) >0) 
  { work <- gBuffer(work, byid=TRUE, width=0)}else{ }

# crop to zone of pred
hab <- crop(hab,b.proj.area)
work <-crop(work,b.proj.area)

## unified field NUM
########################
## creat new class
hab@data <- cbind(hab@data, NUM =rep(1,time=nrow(hab@data) ))
hab$NUM <- as.numeric(hab[,names(hab)==h.hab]@data[,1])
  
## creat new class
work@data <- cbind(work@data, NUM =rep(1,time=nrow(work@data) ))
work$NUM <- as.numeric(work[,names(work)==h.work]@data[,1])

## creat comb set
comb <- bind(hab,work)
comb$NUM<-as.numeric(comb$NUM)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
r.empty <- raster(res=res(b.proj.area), ext=extent(b.proj.area), 
                  crs= crs(comb) ) 

gs <- gstat(formula=NUM~1, locations=comb, set=list(idp = 2))
r.dur <- interpolate(r.empty, gs)
r.dur <- raster::mask(r.dur,b.proj.area)

surf <- raster(paste(working_dir,"zone/surf_terr.tif",sep=""))
w.dur <- winmove(r.dur, d=100, type = c("circle"), win_fun=sum)/surf
                
w.dur <- crop(w.dur,proj.area)

## if raster still do not have the same exact extent
if ( extent(w.dur) == extent(proj.area) ) 
    { }else{
      empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
      w.dur <- resample(w.dur,empty.raster,method="ngb")
    }

w.dur <- raster::mask(w.dur,proj.area)

names(w.dur) <- "humcon"

writeRaster(w.dur, filename = paste(outpath,"humcon.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(w.dur);plot(w.dur)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

#proc.time() - ptm # check time

#####################################################################################################################################
################  END SCRIPT ###############################################
###########################################################################

