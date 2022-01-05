###############################################################################                                                                           #
#=============================================================================#
#    FUNCTION TO PREDICT MODELS THAT ARE FIT USING THE run4models FUNCTION    #
#=============================================================================#
## HISTORIC:
## --------                                                                             #                                                                          #
#   WRITTEN BY RAFAEL W?EST, DEC. 2011, WSL BIRMENSDORF, SWITZERLAND           #
#   rafael.wueest@wsl.ch													  #
#	UPDATED: Christophe F. Randin, May 2012, Uni. Basel, Switezrland 
# v3.0 used for Droz et al. Landscape Urban Plann. 2019
# v4.0 B. Droz February 2021
#  v5.1 --> allow zonal selection based on predictor boundarie used for model calibration
## optimized to extent prediction projection and for parrallel work 
#                                                                             #
###############################################################################
####################################################################################################################################
## DESCRIPTION
###----------
##                    ####################################
##                    ####################################
##                    VERSION FOR PC running under R 3.5.0
###                   ####################################
##                    ####################################
## NEED FOR RUNNING:
##------------------
## Install Java from here http://www.filehippo.com/download_jre_64/
## work with Java Runtime Environment 64-bit 8.0-build-271
##
## Install Maxent v3.3.3k from https://www.cs.princeton.edu/~schapire/maxent/
## dezip the file in your personnal package library
## like "R/win-library/3.5/dismo/java/maxent.jar"
##############################################################################
##
## SET SCRIPT PARAMETER
## ====================
##############################################################################

species=c("RaFB") # model name
# name of scenario -- 
## defined the name of the folders used if all specified check folder automatically
scenario="all" # pred_cons0.4" #"cur" "all"
working_dir <- "~/Actual/rafb_proj_model/analysis/JuraVD_2021/"

# restrict or not the area based on predictor boundarie used for model calibration 
area.pred <- "NO" 

##############################################################################
##############################################################################
#### SCRIPT START HERE
###########################
## Load additional function
###########################
setwd(working_dir)

source("script/function/packagesv1.0.R")
source("script/function/Top.models.proj.v5.1.RaFB.R")
source("script/function/ensemble_projv5.1.Rafb.R")
source("script/function/creat_subDirv1.R")
source("script/function/restri.to_cali_areav1.R")

##################
## Load R. library
##################
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
Sys.setenv(NOAWT=TRUE)

check.lib (c("rms","gam","gbm","rJava","dismo","raster","rgdal","maptools") )
maxent()
########################################################################
beginCluster() # ACTIVATE THIS MULTI CORE CALCULATION

#file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 
#########################################################################
## Define projected scenario
############################
if (scenario=="all") {
    scen <- list.dirs(path=paste(working_dir,"pred",sep=""), recursive = FALSE)
    }else {scen <-paste(working_dir,"pred/",scenario, sep="") }

for (i in 1:length(scen))
  {
  #########################################################################################################################################################
  # Load predictors
  ###################
  path.pred <- scen[i]
  #pred.list <- list.files(path.pred,pattern=".tif$",full.names = TRUE)
  
  # Names of selected predictors
  names.pred <- c("ddur","dtree","dterrenu","dveghaut","dvegras","fsumrad_2m","humcon","lengthfas","mfroadw")
  ddur <-raster(paste(path.pred,"/ddur.tif",sep="") )
  dtree <-raster(paste(path.pred,"/dtree.tif",sep=""))
  dterrenu <-raster(paste(path.pred,"/dterrenu.tif",sep=""))
  dveghaut <-raster(paste(path.pred,"/dveghaut.tif",sep=""))
  dvegras <-raster(paste(path.pred,"/dvegras.tif",sep=""))
  fsumrad_2m <-raster(paste(path.pred,"/fsumrad_2m.tif",sep=""))
  humcon <-raster(paste(path.pred,"/humcon.tif",sep=""))
  lengthfas <-raster(paste(path.pred,"/lengthfas.tif",sep=""))
  mfroadw <-raster(paste(path.pred,"/mfroadw.tif",sep=""))
  # Create predictor stack
  env.stack.proj <- stack(ddur,dtree,dterrenu,dveghaut,dvegras,fsumrad_2m,humcon,lengthfas,mfroadw)
  
  # Create predictor stack
  #env.stack.proj <- stack(pred.list)
  #names.pred <- names(env.stack.proj)

  # set location of calibrate model and projection
  ################################################
  inpath <- paste(working_dir,"script/geographic_proj_model/",sep="")
  envstack <- env.stack.proj
  proj.name <- strsplit(scen[i],"_")
  proj.name <- proj.name[[1]][length(proj.name[[1]])]
  creat.subDir(working_dir,paste("proj/proj_",proj.name,sep=""))
  outpath.proj <- paste(working_dir,"/proj/proj_",proj.name,"/", sep="")
  ######################################################################
  ## Test crop  extent
  #e <- extent(c(2546000, 2548000, 1211000, 1212000))
  #r.test <- raster(e,resolution=2)

  #env.stack.proj <- crop(env.stack.proj,r.test)
  #envstack=env.stack.proj
  ################################################################################
  ########### MODEL RUN  #########################
  ###############################################################################
  # ptm <- proc.time()# ignite timer

  ### restrict area of prediction
  if (area.pred =="YES"){
  
    f.path <- paste(working_dir,"script/geographic_proj_model/boundaries_calibrate_pred.txt",sep="")
  
    env.stack.proj <- restri.to_cali_area (f.path,env.stack.proj)
  
  }else{}

  # Run code
  top.models.proj(inpath,GLM.mod=TRUE,GLM.ModAv=TRUE,
                GLM.Step=TRUE,GAM.Step=TRUE,GBM.mod=FALSE, 
                ME.mod=TRUE,
                species=c("RaFB"),envstack=env.stack.proj,
	              outpath.proj,scenario=proj.name)

  ensemble.spat.proj(inpath=outpath.proj,
                  species="RaFB",scenario=proj.name,
                  glm.modav.mod=TRUE,n.rep.glm.modav=1,
                  glm.stepwise.mod=TRUE,n.rep.glm.stepwise=1,	
                  gam.stepwise.mod=TRUE,n.rep.gam.stepwise=1,
                  gbm.mod=FALSE,n.rep.gbm=10,
                  rf.mod=FALSE,n.rep.rf=10,
                  me.mod=TRUE,n.rep.me=10,
                  ens.code="E4",
                  outpath=outpath.proj)
  
  file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

  }

endCluster() # END OF MULTICORE CALCULATION

#proc.time() - ptm # check time

####################  END #####################################################
