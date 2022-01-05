####################################################################################################################################
###
###     --- Reclassified green area based on ndvi calibrate model  ------
#####============================================================

# v1.0 - Boris Droz - April 2021
####################################################################
## Descriptif
##############
#  combine land cover data and orthophoto detection...
## classification follow : Droz et al. Ornis Fenn. 2015, 92, (3), 112-122.
## 
## Orthophoto classification based on ndvi calibrate in Cdfs with rgbir
## ORTHOPHOTO from https://www.swisstopo.admin.ch/fr/geodata/images/ortho/swissimage-rs.html#d_tails_techniques
## classied as band 1.NIR 2 Red, 3green 4 Blue !! in 16bit
## manually calibrate on a sample of orthophoto from same fly coverage
################################################################################
## Input : bench of geotif in vignette 
##         should be in a "ortho" folder in the input folder of your project
## Output: similar data as input with lower resolution
##          creat a "rgbir" folder in the input folder of your project

##   ----- FIRST ---
## !!! check classification attribute and shp field in the shp input file
## !!  adjust it if necessary in above 163 in the script
## !!Then set the parameter and run the script
##
###############################################################################
####################################################################################################################################
## SCRIPT PARAMETER
########################
# Set folder 
working_dir <- "~/Actual/rafb_proj_model/analysis/Locle_2021/"

#FINAL RESOLUTION
end.res <- 2   # in m if GRS80 degre if wgs84
proj <- 2056 # input projection new swiss projection 21781 # old swiiss

#select file and categories for vegetation and terrenue
f.landuse<- "mo6_couverture_du_sol.shp" 
h.lab.veg <- "desnat"
cat.veg <- c("pâturage","jardin","tourbière",
             "pâturage boisé dense","tourbière boisée")

h.lab.terre <- "desnat"
cat.terre <- c("sans végétati")

####################################################################################################################################
## START SCRIPT 
#####################
## Load library
##################}
setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","grainchanger","rgeos"))

# set file location 
inpath <- paste(working_dir,"/input/ortho",sep="")
creat.subDir(working_dir,"/input/rgbir")

outpath <- paste(working_dir,"/va/rgbir",sep="")
setwd(inpath)
################################################################################
## --- PARRALLEL SET UP----
############################
beginCluster() # ACTIVATE THIS MULTI CORE CALCULATION 

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 
###################################################################################
## OPEN Data
################
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )

# open land use
c.sol <- rgdal::readOGR(paste(working_dir,"/input/",f.landuse,sep=""))
proj4string(c.sol) <- CRS(paste("+init=epsg:", proj,sep="") )

#################################################################################
### 1) Lower the resoltuion of all picture by align with raster project
#######################################################################

rast.list <- list.files(inpath,pattern=".tif$",full.names = TRUE) #; print(rast.list) # get path for each raster file
f.name.list <- list.files(inpath,pattern=".tif$",full.names = FALSE)

for (z in 1 :length(rast.list))
  {
    cat("> Vignet ",z,"of",length(rast.list)," is aggregated now...", "\n",append = FALSE)
    
    pred.var <- stack(rast.list[z]) #; print(pred.var) ## used stack instead of raster for multi layer raster
    proj4string(pred.var) <- CRS(paste("+init=epsg:", proj,sep="") )# gave crs values
    
    # check if inside area of interest
    if ( extent(b.proj.area)[1] <= extent(pred.var)[2] &
         extent(b.proj.area)[2] >= extent(pred.var)[1] &
         extent(b.proj.area)[3] <= extent(pred.var)[4] &
         extent(b.proj.area)[4] >= extent(pred.var)[3] )
        {
        
        pred.var <- crop(pred.var, extent(b.proj.area))# refined to the proj raster
        pred.var <- aggregate(pred.var,fact=end.res/res(pred.var)[1],mean) #agregate 
        pred.var <- projectRaster(from=pred.var, to=b.proj.area, method="bilinear")# to match with the project grid
          
    	 writeRaster(pred.var, filename = paste(working_dir,"/input/rgbir/",f.name.list[z],sep=""),
                datatype="FLT8S", overwrite=TRUE)  
        }else{}
    }

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#################################################################################
### 2) Merge all vignette of the area and set origine to project buffer
#######################################################################

# get list of file
fns <- list.files(paste(working_dir,"input/rgbir/", sep="") , pattern=".tif$",full.names = TRUE) 

## merge all demain the zone of interest
p<-0

temp.v <-NULL
for (c  in 1: length(fns))
{
  cat("> File",c,"of", length(fns), "is merged ...", "\n",append = FALSE)
  r <- stack(fns[c]) # open raster 4 bands
  
  # check if inside area of interest
  if ( extent(b.proj.area)[1] <= extent(r)[2] &
       extent(b.proj.area)[2] >= extent(r)[1] &
       extent(b.proj.area)[3] <= extent(r)[4] &
       extent(b.proj.area)[4] >= extent(r)[3] )
  { 
    p<- p+1
    if (p==1) {
      dem <- r
    }else{  dem <- mosaic(dem,r,fun=mean)}    
  }else{}
}

dem <- crop(dem,extent(b.proj.area))

## if raster still do not have the same exact extent
empty.raster <- raster(res=res(b.proj.area),ext=extent(b.proj.area), 
                       crs= crs(b.proj.area )) 

if ( extent(dem) == extent(b.proj.area) ) 
  { }else{
    dem <- resample(dem,empty.raster,method="ngb")
  }

dem <- raster::mask(dem,b.proj.area)

writeRaster(dem, filename = paste(working_dir,"/input/rgbir_merged.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

## Calcul NDVI
# NDVI = (NIR - Red) / (NIR + Red)
r.ndvi <- (dem[[1]]-dem[[2]])/(dem[[1]]+dem[[2]])

#writeRaster(r.ndvi, filename = paste(inpath,"/ndvi.tif",sep=""),
 #           datatype="FLT8S", overwrite=TRUE)

#################################################################################
### 3) select area for bareground data and green 
#######################################################################
outpath <- paste(working_dir,"/pred/pred_cur/",sep="")

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(c.sol, byid=TRUE)==FALSE) >0) 
{ c.sol <- gBuffer(c.sol, byid=TRUE, width=0)}else{ }

# crop to zone of pred
c.sol <-crop(c.sol,b.proj.area)

## CREAT RASTER SET dvegrase + dterrenu
#######################################
## typcou == verte desnat == pr�-champ 
if (h.lab.veg=="NO") { veg <-c.sol}else{
  pos <- as.logical(c.sol[,names(c.sol)==h.lab.veg]@data==cat.veg)
  if (length(cat.veg) >1){
    for (i in 2:length(cat.veg))
    {pos<- pos| as.logical(c.sol[,names(c.sol)==h.lab.veg]@data==cat.veg[i] ) }
  }else{}
  veg <- c.sol[pos,]
  veg@data <- cbind(veg@data, rep(1,time=nrow(veg@data) ))
  names(veg) <- c(names(c.sol),"val") }

# select terrenue area
if (h.lab.terre=="NO") { terrenu <-c.sol}else{
  pos <- as.logical(c.sol[,names(c.sol)==h.lab.terre]@data==cat.terre)
  if (length(cat.terre) >1){
    for (i in 2:length(cat.terre))
    {pos<- pos| as.logical(c.sol[,names(c.sol)==h.lab.terre]@data==cat.terre[i] ) }
  }else{}
  terrenu <- c.sol[pos,]
  terrenu@data <- cbind(terrenu@data, rep(1,time=nrow(terrenu@data) ))
  names(terrenu) <- c(names(c.sol),"val") }

# rasterize
r.veg <- rasterize(veg,b.proj.area,field="val",background=0)
r.veg <- mask(r.veg,b.proj.area)

#writeRaster(r.veg, filename = paste(inpath,"test.r.veg.tif",sep=""),
 #           datatype="FLT8S", overwrite=TRUE)

r.terrenu <- rasterize(terrenu,b.proj.area,field="val",background=0)
r.terrenu <- mask(r.terrenu,b.proj.area)

### use NDVI calibrate model to reclassified green area into bareground or green
################################################################################
# predict green area

load(paste(working_dir,'script/ndvi_model/glm_ndvi_model',sep='')) # load cali model

names(r.ndvi) <- "NDVI"
r.veg.ndvi <- predict(r.ndvi, glm.ndvi, type="response",na.rm = TRUE)
  
# reclassified  with max TSS
m <- as.numeric (c(0, glm.stepwise.max.tss.tresh, 0, glm.stepwise.max.tss.tresh, 1, 1))
r.veg.ndvi <- reclassify(r.veg.ndvi, m)

#writeRaster(r.veg.ndvi, filename = paste(inpath,"/test.r.vegjardin.tif",sep=""),
 #           datatype="FLT8S", overwrite=TRUE)

r.veg.ndvi <- r.veg.ndvi*r.veg

r.veg.terrenu <- r.veg - r.veg.ndvi

r.terrenu <- r.terrenu + r.veg.terrenu

r.veg <- r.veg.ndvi

######## do mowing window
#########################
## mowing window 
surf <- pi*100^2
dvegras <- winmove(r.veg, d=100, type = c("circle"), win_fun=sum)/
  (surf/(res(r.veg)[1]^2))

dterrenu <- winmove(r.terrenu, d=100, type = c("circle"), win_fun=sum)/
  (surf/(res(r.terrenu)[1]^2))

# resize to project
dvegras <- crop(dvegras,proj.area)
dterrenu <- crop(dterrenu,proj.area)

## if raster still do not have the same exact extent
if ( extent(dvegras) == extent(proj.area) ) 
  { }else{
    empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
    dvegras <- resample(dvegras,empty.raster,method="ngb")
  }

if ( extent(dterrenu) == extent(proj.area) ) 
  { }else{
    empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
    dterrenu <- resample(dterrenu,empty.raster,method="ngb")
  }

dvegras <- mask(dvegras,proj.area)
dterrenu <- mask(dterrenu,proj.area)

names(dvegras) <- "dvegras"
names(dterrenu) <- "dterrenu"

writeRaster(dvegras, filename = paste(outpath,"dvegras.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)
writeRaster(dterrenu, filename = paste(outpath,"dterrenu.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
par(mfrow=c(1,2)) 

print(dvegras);plot(dvegras)
print(dterrenu);plot(dterrenu,new=TRUE)

dev.off() 

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

##################################################################################################################
##################################################################################################################
##
## ------------------ END ------------------
##
##################################################################################################################
##################################################################################################################
