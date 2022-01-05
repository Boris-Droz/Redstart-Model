#####
#### Generation dveghaut
####=================================================
####
## work under R 4.0.3
## v1.1 B.Droz February 2021

###################################################################################################
## Descriptif
##############
#  define dveghaut as all meadow, pasture and field present in the area
# similar classification as Droz et al. Ornis Fenn. 2015, 92, (3), 112-122. 
## desnat == pré-champ 
## from couverture_du_sol
##
## First, check classification attribute and shp field in the shp input file
## adjust the header label and categories on parameter and run the script

#####################################################################################################################################
##### PARAMETER #################
####################################################################

working_dir <-"~/Actual/rafb_proj_model/analysis/JuraVD_2021/" 

#select...  file, header and lancover class
f.landuse<- "mo6_couverture_du_sol.shp" # file name
h.lab.landuse <- "desnat" # header of the landcover name
                          # if select "NO" no selection is made
cat.landuse <- c("prÃ©-champ") # classes of land cover
  
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

check.lib (c("raster","sp","magrittr","rgdal","rgrass7","RNetCDF",
             "splitstackshape","maptools","XML","grainchanger","gstat"))

use_sp()

# set file location 
inpath <- paste(working_dir,"/input",sep="")
creat.subDir(working_dir,"/pred/pred_cur")
outpath <- paste(working_dir,"/pred/pred_cur/",sep="")
setwd(inpath)
################################################################################
beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#####################################################################################################################################
## OPEN Data
################
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )

c.sol <- rgdal::readOGR(f.landuse)

# set the coordinate system
proj4string(c.sol) <- CRS(paste("+init=epsg:", proj,sep="") )

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(c.sol, byid=TRUE)==FALSE) >0) 
  { c.sol <- gBuffer(c.sol, byid=TRUE, width=0)}else{ }

# crop to zone of pred
c.sol <-crop(c.sol,b.proj.area)

## CREAT RASTER SET dveghaut
########################
## typcou == verte desnat == pré-champ 
if (h.lab.landuse=="NO") { veg <-c.sol}else{
  pos <- as.logical(c.sol[,names(c.sol)==h.lab.landuse]@data==cat.landuse)
  if (length(cat.landuse) >1){
    for (i in 2:length(cat.landuse))
        {pos<- pos| as.logical(c.sol[,names(c.sol)==h.lab.landuse]@data==cat.landuse[i] ) }
    }else{}
  veg <- c.sol[pos,]
  veg@data <- cbind(veg@data, rep(1,time=nrow(veg@data) ))
  names(veg) <- c(names(c.sol),"val") }

# rasterize
r.veg <- rasterize(veg,b.proj.area,field="val",background=0)
r.veg <- r.veg *b.proj.area

surf <- raster(paste(working_dir,"zone/surf_terr.tif",sep=""))
w.veg <- winmove(r.veg, d=100, type = c("circle"), win_fun=sum, na.rm = TRUE)/ surf

w.veg <- crop(w.veg,proj.area)

## if raster still do not have the same exact extent
if ( extent(w.veg) == extent(proj.area) ) 
    { }else{
      empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
      w.veg <- resample(w.veg,empty.raster,method="ngb")
    }

w.veg <- raster::mask(w.veg,proj.area)

names(w.veg) <- "dveghaut"

writeRaster(w.veg, filename = paste(outpath,"dveghaut.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(w.veg);plot(w.veg)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

###########################################################################
############################################################################
###########################################################################
###########################################################################

