#####
#### Generation ddur
####=================================================
####
## work under R 4.0.3
## v1.2 B.Droz February 2021

###################################################################################################
## Descriptif
##############
#  define ddur as the sum of all seatment area
## similar classification as Droz et al. Ornis Fenn. 2015, 92, (3), 112-122.
## ddur = typcouv == dur from couverture_du_sol and batiment + road

## First, check classification attribute and shp field in the shp input file
## adjust the header label and categories on parameter and run the script
##
#####################################################################################################################################
##### PARAMETER #################
####################################################################

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
n.proj <- "LV95" # name projection
#f.type <- ".asc$" # type of file of the mnt
working_dir <-"~/Actual/rafb_proj_model/analysis/Cdf_2021/" 

c
f.landuse<- "mo6_couverture_du_sol.shp" # file name
h.lab.landuse <- "typcou" # header of the landcover name
                          # if select "NO" no selection is made
cat.landuse <- c("dure") # class of land cover

f.bat <- "mo22_batiments.shp"
h.lab.bat <- "typcou"    # if select "NO" no selection is made
cat.bat <- c("ordinaire")
  
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

check.lib (c("raster","sp","magrittr","rgdal","rgeos","RNetCDF",
             "splitstackshape","maptools","XML","grainchanger","gstat"))

# set file location 
inpath <- paste(working_dir,"input/",sep="")
creat.subDir(working_dir,"/pred/pred_cur")
outpath <- paste(working_dir,"/pred/pred_cur/",sep="")
setwd(inpath)
###################################################################################
beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#####################################################################################################################################
## OPEN Data
################
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )

bat <- rgdal::readOGR(f.bat)
c.sol <- rgdal::readOGR(f.landuse)

# set the coordinate system
proj4string(bat) <- CRS(paste("+init=epsg:", proj,sep="") )
proj4string(c.sol) <- CRS(paste("+init=epsg:", proj,sep="") )

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(bat, byid=TRUE)==FALSE) >0) 
  { bat <- gBuffer(bat, byid=TRUE, width=0)}else{ }

if (sum(gIsValid(c.sol, byid=TRUE)==FALSE) >0) 
  { c.sol <- gBuffer(c.sol, byid=TRUE, width=0)}else{ }

# crop to zone of pred
bat <- crop(bat,b.proj.area)
c.sol <-crop(c.sol,b.proj.area)

## CREAT RASTER SET ddur
########################
# gave number to each typcou
# 1:dure, 2:verte, 3:bois, 4:eau
## select the dur land use
if (h.lab.landuse=="NO") { dur <-c.sol}else{
  # dur <-c.sol[c.sol$typcou=="dure",]
  pos <- as.logical(c.sol[,names(c.sol)==h.lab.landuse]@data==cat.landuse)
  dur <- c.sol[pos,]
  dur@data <- cbind(dur@data, rep(1,time=nrow(dur@data) ))
  names(dur) <- c(names(c.sol),"val") }

# remove souterrain from batiment
if (h.lab.bat=="NO") { hous <-bat}else{
  #hous <- bat[bat$typcou=="ordinaire",]
  pos <- as.logical(bat[,names(bat)==h.lab.bat]@data==cat.bat)  
  hous <- bat[pos,]
  hous@data <- cbind(hous@data, rep(1,time=nrow(hous@data) ))
  names(hous) <- c(names(bat),"val") }

# rasterize
r.dur <- rasterize(dur,b.proj.area,field="val",background=0)
r.bat <- rasterize(hous,b.proj.area,field="val",background=0)

r.dur <- r.dur+r.bat
r.dur <- reclassify(r.dur,c(1,Inf,1))*b.proj.area

#surf <- pi*100^2
surf <- winmove(b.proj.area, d=100, type = c("circle"), win_fun=sum, na.rm = TRUE)
w.dur <- winmove(r.dur, d=100, type = c("circle"), win_fun=sum, na.rm = TRUE)/ surf
                #(surf/(res(r.dur)[1]^2)) 

# keep surf for the other pred var calculus
writeRaster(surf, filename = paste(working_dir,"zone/surf_terr.tif",sep=""),
            +             datatype="FLT8S", overwrite=TRUE)

w.dur <- crop(w.dur,proj.area)

## if raster still do not have the same exact extent
if ( extent(w.dur) == extent(proj.area) ) 
    { }else{
      empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
      w.dur <- resample(w.dur,empty.raster,method="ngb")
    }

w.dur <- raster::mask(w.dur,proj.area)

names(w.dur) <- "ddur"

writeRaster(w.dur, filename = paste(outpath,"ddur.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(w.dur);plot(w.dur)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

###########################################################################
############################################################################
###########################################################################
###########################################################################

