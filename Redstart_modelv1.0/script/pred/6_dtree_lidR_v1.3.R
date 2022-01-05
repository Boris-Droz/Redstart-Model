#####
#### Generation dtree -- canopy coverage
### 
####=================================================
####
## work under R 4.0.3
## v1.3 B.Droz February 2021
##
###################################################################################################
## Descriptif
##############
#  Procedure from :"https://jean-romain.github.io/lidRbook/index.html"
## Roussel et al. Remote Sens. Environ. 2020, 251, 112061. 
###  https://doi.org/10.1016/j.rse.2020.112061

## Classifcation of lidar data follow:
## https://www.swisstopo.admin.ch/fr/geodata
## /height/surface3d.html#technische_details

#############################################################################
## Need to creat a distingue folder named las with the las data in input folder
## Then set the parameter and run the script
## 
#####################################################################################################################################
##### PARAMETER #################
####################################################################

f.type <- ".las$" # type of file of the lidar
working_dir <-"~/Actual/rafb_proj_model/analysis/JuraVD_2021/" 

#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
## parameter for canopy detection
#################################

pr.res <- 0.5 # Point-to-raster resolutions --> resolution of the lidar
wind <- 5 # window routine procedure
subcir <- 0.2 #sucircle tweak

#####################################################################################
#####################################################################################
## Load library
##################

setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","rgdal","lidR",
             "maptools", "grainchanger","rgeos" ))
#use_sp()

# set file location 
inpath <- paste(working_dir,"/input/",sep="")
creat.subDir(inpath,"canopycover")
creat.subDir(working_dir,"/pred/pred_cur")
outpath <- paste(working_dir,"pred/pred_cur/",sep="")
setwd(inpath)
###################################################################################
beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#####################################################################################################################################
## OPEN Data
################
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )

#ptm <- proc.time()# ignite timer

#######################################
#  COMPUTE LAS DATA
#######################################
# Creat list of raster file and open it
fns <- list.files(paste(inpath,"/las/",sep="")
                  ,pattern=f.type,full.names = TRUE)

f.name <-list.files(paste(inpath,"/las/",sep="")
                      ,pattern=f.type,full.names = FALSE)

## merge all las
for (c  in 1: length(fns))
{
  cat("> Canopycoverage number",c," is compute...", "\n",append = FALSE)
  
  las <- readLAS(fns[c],select = "xyzr", 
                  filter = "-drop_z_below 5 -keep_class 3 -keep_NDVI 0.5 1.0")
  
  crs(las) <- crs(b.proj.area)
  
  name <- unlist( strsplit(f.name[c],"[.]") )[1]

  # check if inside area of interest
  if ( extent(b.proj.area)[1] <= extent(las)[2] &
       extent(b.proj.area)[2] >= extent(las)[1] &
       extent(b.proj.area)[3] <= extent(las)[4] &
       extent(b.proj.area)[4] >= extent(las)[3] )
    {
    
  #######################################################################
  ############# DELIMIT CROWNS --- START
  ## #  Procedure from :"https://jean-romain.github.io/lidRbook/index.html"
  ########################################################################
  # Point-to-raster 2 resolutions
  chm_p2r <- grid_canopy(las, pr.res, p2r(subcircle = subcir))

  ## Then the same tree detection routine with a constant windows size of w in m is applied to each CHM:
  ttops_chm_p2r <- find_trees(chm_p2r, lmf(w=wind))

  ## Segmentation of the point-cloud
  algo_p2r <- dalponte2016(chm_p2r,ttops_chm_p2r)

  # segment point cloud
  las_p2r <- segment_trees(las, algo_p2r) 

  # det crowns
  crowns_p2r <- delineate_crowns(las_p2r)
          
  # for each crown ad value of one
  crowns_p2r@data <- cbind(crowns_p2r@data,val=rep(1,nrow(crowns_p2r@data)))

  # raseterize
  r.crowns_p2r <- rasterize(crowns_p2r,b.proj.area,field="val",background=0)

  #r.crowns_p2r <- mask(r.crowns_p2r,proj.area)
  r.crowns_p2r <- reclassify(r.crowns_p2r,c(1,Inf,1))

  # write tempory raster
  writeRaster(r.crowns_p2r, filename = paste(inpath,"canopycover/cc_",name,".tif",sep=""),
              datatype="FLT8S", overwrite=TRUE)

    } else { }
  }
##############################################################################
################
# Creat list of raster file and open it
fns <- list.files(paste(inpath,"canopycover/",sep="")
                  ,pattern=".tif$",full.names = TRUE)

## merge all canopy
r.set <- calc(stack(fns),sum)
r.set <- raster::mask(r.set,b.proj.area)

#r.set <- reclassify(r.set, c(0,Inf,1)) 

writeRaster(r.set, filename = paste(inpath,"cc",pr.res,"_LV95_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

##############################################################################
## mowing window 
surf <- raster(paste(working_dir,"zone/surf_terr.tif",sep=""))
w.crowns_p2r <- winmove(r.set, d=100, type = c("circle"), win_fun=sum)/surf

w.crowns_p2r <- crop(w.crowns_p2r, proj.area)

## if raster still do not have the same exact extent
if ( extent(w.crowns_p2r) == extent(proj.area) ) 
    { }else{
      empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
      w.crowns_p2r <- resample(w.crowns_p2r,empty.raster,method="ngb")
    }

w.crowns_p2r <- raster::mask(w.crowns_p2r, proj.area)

writeRaster(w.crowns_p2r, filename = paste(outpath,"dtree.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(w.crowns_p2r);plot(w.crowns_p2r)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

#proc.time() - ptm # check time

###############################################################################
#########################################
##### END COFFE TIME