################################################################################################
################################################################################################
###  					  ###					
###    CREAT SHAPE ZONE CONSERVATION AND SHAPE PROMOTION OF DIVERCITY!!!
###    ================================================================
#
################################################################################################
################################################################################################
## work under R 4.0.3
## v1.0  March 2021
###############################################################################################
## DESCRIPTION
## ===========
## Make simple nicer smooth shape file for conservation 
##  see packege help on https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html
###############################################################################################
## SET SCRIPT PARAMETER
## ====================

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
species=c("RaFB") # model name
working_dir <- "~/Actual/rafb_proj_model/analysis/Planchette_2021/"

surf.min <- 31400 ## surface min in m2 for considering a area (less than a territory)
fil.hole <- 20000 # m2 hole are filled 
d.terri <- 50 # distance in m to connect two features (less than a territory)
smooth.numb <- 100
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
source("script/function/clean.rasttopolyv2.R")


check.lib (c("raster","rgdal","maptools","geosphere", 
             "spdep", "sp", "cpr", "smoothr", "sf","units", "rgeos", "rangemap"))

# set path
creat.subDir(working_dir,"proj/proj_plan.action")
outpath <- paste(working_dir,"proj/proj_plan.action",sep="")

################################################################################
# load raster
#############
cur.proj <- raster(paste(working_dir,
                         "/proj/proj_cur/ENSEMBLE_bin_RaFB_E4_cur.tif", sep="") )

cons.proj <- raster(paste(working_dir,
                         "/proj/proj_cons0.4/ENSEMBLE_bin_RaFB_E4_cons0.4.tif", sep=""))

## creat a raster with area with value equal to 4
cur.proj <- reclassify(cur.proj, c(-Inf,3,NA, 3,Inf,1))
cons.proj <- reclassify(cons.proj, c(-Inf,3,NA, 3,Inf,1))

crs(cur.proj) <- CRS(paste("+init=epsg:", proj,sep="") )
crs(cons.proj) <- CRS(paste("+init=epsg:", proj,sep="") )

s.cur.proj <- clean.rasttopoly(cur.proj,
                               surf.min=surf.min,fil.hole=fil.hole, # both in m2
                               method = "ksmooth", nb.smooth=smooth.numb )

s.cons.proj <- clean.rasttopoly(cons.proj,
                               surf.min=surf.min,fil.hole=fil.hole, # both in m2
                               method = "ksmooth", nb.smooth=smooth.numb)

# write zone
area.proj <- strsplit(working_dir,"/")[[1]][length(strsplit(working_dir,"/")[[1]])]
shapefile(s.cur.proj, paste(outpath,"/",species,"_cons_",area.proj,".shp",sep=""),
          overwrite=TRUE)
shapefile(s.cons.proj, paste(outpath,"/",species,"_promo_",area.proj,".shp",sep=""),
          overwrite=TRUE)
