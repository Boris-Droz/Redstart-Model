#####
#### Generation lengthfas -- length of building wall
####=================================================
####
## work under R 4.0.3
## v1.1 B.Droz February 2021
###################################################################################################
## Descriptif
##############
##  description in 
## " https://geocompr.robinlovelace.net/geometric-operations.html"

## First, check classification attribute and shp field in the shp input file
## adjust it if necessary in line 65 in the script
## Then set the parameter and run the script
## 
#####################################################################################################################################
##### PARAMETER #################
####################################################################

working_dir <-"~/Actual/rafb_proj_model/analysis/JuraVD_2021/" 

# class of land cover
f.bat <- "mo22_batiments.shp" # file names
h.lab.bat <- "typcou"    # header of the landcover name
                          # if select "NO" no selection is made
cat.bat <- c("ordinaire") # classes of land cover

proj <- 2056

#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################
setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","magrittr","rgdal","maptools", "grainchanger",
             "sf","dplyr", "rgeos"))

# set file location 
inpath <- paste(working_dir,"/input",sep="")
creat.subDir(working_dir,"/pred/pred_cur")
outpath <- paste(working_dir,"/pred/pred_cur/",sep="")
setwd(inpath)
###################################################################################
# ptm <- proc.time()

beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#####################################################################################################################################
## OPEN Data
################
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )

bat <- rgdal::readOGR(f.bat)

# set the coordinate system
proj4string(bat) <- CRS(paste("+init=epsg:", proj,sep="") )

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(bat, byid=TRUE)==FALSE) >0) 
  { bat <- gBuffer(bat, byid=TRUE, width=0)}else{ }

# crop to zone of pred
bat <- crop(bat,b.proj.area)

# remove souterrain from batiment
# remove souterrain from batiment
if (h.lab.bat=="NO") { bat <-bat}else{
  #hous <- bat[bat$typcou=="ordinaire",]
  pos <- as.logical(bat[,names(bat)==h.lab.bat]@data==cat.bat)  
  bat <- bat[pos,] }

## creat buffer around each bat
b.bat <- buffer(bat,width=0.5,dissolve=FALSE)

d.bat <- rgeos::gDifference (b.bat,bat)

# creat dataframe with correct row name
p.df <- data.frame( ID=1:length(d.bat), row.names = 1)

d.bat <- SpatialPolygonsDataFrame(d.bat,p.df,match.ID = FALSE)

# rasterize
empty.raster <- raster(res=c(1,1),ext=extent(b.proj.area), crs= crs(b.proj.area ))
r.bat <- rasterize(d.bat,empty.raster,field="ID",background=0)
r.bat <- reclassify(r.bat,c(1,Inf,1))

w.bat <- aggregate(r.bat,fact=2/1,sum) # end res / actu
w.bat <- w.bat*b.proj.area

# do surface mowing window
surf <- raster(paste(working_dir,"zone/surf_terr.tif",sep=""))
w.bat <- winmove(w.bat, d=100, type = c("circle"), win_fun=sum)/surf

w.bat <- crop(w.bat,proj.area)

## if raster still do not have the same exact extent
if ( extent(w.bat) == extent(proj.area) ) 
    { }else{
      empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
      w.bat <- resample(w.bat,empty.raster,method="ngb")
    }

w.bat <- raster::mask(w.bat,proj.area)

names(w.bat) <- "lengthfas"

writeRaster(w.bat, filename = paste(outpath,"lengthfas.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(w.bat);plot(w.bat)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

#proc.time() - ptm
###########################################################################
############################################################################
###########################################################################
###########################################################################

