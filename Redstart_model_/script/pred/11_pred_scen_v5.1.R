################################################################################################
################################################################################################
###  					  ###					
###    CREATE PREDICTOR UNDER CONSERVATION AND THREAT SCENARIO 
###    =======================================================
##
##     Boris DROZ & Christophe RANDIN, University of Basel, BIB
##    Projet Rafb
################################################################################################
################################################################################################
## work under R 4.0.3
## v4.4, Novembre 2013 
## v5.1  March 2021
###############################################################################################
## DESCRIPTION
## ===========
## Allow landscape modification only in a restricted area.
## 
###############################################################################################
## SET SCRIPT PARAMETER
## ====================

species=c("RaFB") # model name
scenario=c("cons") #,"dens") # name of scenario -- defined the name of folder 
                          # two possibilities cons and dens (threat)

working_dir <- "~/Actual/rafb_proj_model/analysis/JuraVD_2021/"

## Parameter scenario
######################################
## proportion increment (supplementary ratio compare to initial)
i = 0.1
## number of increment
p = 4
## Value of Optimal tree proportion - defined during model calibration
opt.dtree = 0.4

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

check.lib (c("raster","rgdal","maptools"))

# set file location -- start
inpath <- paste(working_dir,"pred/pred_cur",sep="")
setwd(inpath)
################################################################################################
######################################
################################################################################################
# Load starting predictors
###################
#ptm <- proc.time()

names.pred <- c("ddur","dtree","dterrenu","dvegras","dveghaut")
ddur <-raster("ddur.tif")
dtree <-raster("dtree.tif")
dterrenu <-raster("dterrenu.tif")
dveghaut <-raster("dveghaut.tif")
dvegras <-raster("dvegras.tif")

# Load plan d'ammenagement - simplified version
## --> AT100 Affectations primaires
###################################################################
dplan.cs <- raster(paste(working_dir,"/zone/zone_pa.tif",sep="") )
dplan.cs -> dplan.ds # similar amenagement plan -- could be changed
# dplan.cs <- raster("solocc_1998cs")
# dplan.ds <- raster("solocc_1998ds")
################################################################################################
################################################################
## Zone de restriction -> PA
rclmat <- matrix(c(-Inf,0,0, 0,Inf,1), ncol=3, byrow=TRUE)
fplan.cs <- reclassify(dplan.cs, rclmat)
fplan.ds <- reclassify(dplan.ds, rclmat)

#####################################
# Restriction des variable a la zone
#####################################
## Variante conservation 
##----------------------
if (any(scenario=="cons") ){
    cs.veghaut <- dveghaut*fplan.cs
    
    cs.dterrenu <- dterrenu*fplan.cs
    
    cs.dvegras <- dvegras*fplan.cs
    
    dtot0 <- cs.dterrenu + cs.veghaut + cs.dvegras
    
    dtree0 <- dtree
}else {}

if (any(scenario=="dens") ){
  ## Variante densification
  ##-----------------------
  mod.veghaut <- dveghaut*fplan.ds
  unmod.veghaut <- dveghaut - dveghaut*fplan.ds
  
  mod.dterrenu <- dterrenu*fplan.ds
  unmod.dterrenu <- dterrenu - dterrenu*fplan.ds
  
  mod.dvegras <- dvegras*fplan.ds
  unmod.dvegras <- dvegras - dvegras*fplan.ds
  
}else {}

## Conservation de valeur initial (mean)
##---------------------------------------
mdur.i <- mean(as.array(ddur),na.rm=TRUE)
mtree.i <- mean(as.array(dtree),na.rm=TRUE)
mterrenu.i <-mean(as.array(dterrenu),na.rm=TRUE)
mveghaut.i <-mean(as.array(dveghaut),na.rm=TRUE)
mvegras.i <- mean(as.array(dvegras),na.rm=TRUE)		
################################################################################################
## Creat new predictor
################################################################################################
## initilise matrice for calculation reel augmentation sur la zone

mat.land.mod <- matrix(0,(p+1),5,dimnames=list(c((0:p)*i), c(names.pred)))

mat.land.mod[1,] <- c(mdur.i,mtree.i,mterrenu.i,mvegras.i,mveghaut.i)
	
for (n in 1:p)
    {  
      if (any(scenario=="cons") ){
        # creat output folder for scenario predictor
        creat.subDir(working_dir,paste("pred/pred_cons",n*i,sep=""))
        outpath <- paste(working_dir,"pred/pred_cons",n*i,"/", sep="")
        ## Conservaton scenario
        ##-----------------------
        r1 <- matrix(c(0, Inf, 1,-Inf, 0, 0), ncol=3, byrow=TRUE)
        r2 <- matrix(c(0, Inf, 0,-Inf, 0, 1), ncol=3, byrow=TRUE)
        
        # No modification on the possible area of modification
        dtree <- dtree*reclassify((dtree-dtot0), r1) + (dtree+i)*reclassify((dtree-dtot0), r2)          
        
        # Account for the optimal value of tree (opt.dtree) - not allow more
        mor.opt <- matrix(c(opt.dtree, Inf, 1,-Inf, opt.dtree, 0), ncol=3, byrow=TRUE)
        les.opt <- matrix(c(opt.dtree, Inf, 0,-Inf, opt.dtree, 1), ncol=3, byrow=TRUE)
        
        dtree <- dtree*reclassify(dtree,les.opt) + dtree0*reclassify(dtree, mor.opt) 
          
        dtree0 <- dtree 
        
        ## write output for conservation scenario  
        writeRaster(dtree, filename=paste(outpath,"dtree.tif",sep=""),
                    datatype="FLT8S", overwrite=TRUE)
        
        # copy all other predictor
        l.pred <- c("ddur","dterrenu","dvegras","dveghaut","fsumrad_2m",
                   "humcon","lengthfas","mfroadw" )
        for (j in 1:length(l.pred))
            { 
            file.copy(paste(inpath,"/",l.pred[j],".tif",sep=""), outpath )
            }
        
      }else{ }
  
    if (any(scenario=="dens") ){
        # creat output folder for scenario predictor
        creat.subDir(working_dir,paste("pred/pred_dens",n*i,sep=""))
        outpath <- paste(working_dir,"pred/pred_dens",n*i,"/", sep="")
        
        ## Densification scenario
        ##-----------------------
	      ## definition pour chaque variable quand on peut densifier
        r4 <- matrix(c(-Inf,i/3,0, i/3,1,1, 1,Inf,0), ncol=3, byrow=TRUE)
        r6 <- matrix(c(-Inf,0,0, 0,1,1, 1,Inf,0), ncol=3, byrow=TRUE)

        rterrenu <- reclassify(mod.dterrenu, r4) 
        rveghaut <- reclassify(mod.veghaut, r4)
        rvegras <- reclassify(mod.dvegras, r4)
	      rdur <- reclassify(ddur, r6)	
        
        ## dens.possible when terrenu, veghaut ou vegrase disponible
        ## but zone de densif comprend aussi ddur
        r7 <- matrix(c(-Inf,1,0, 1,4,1), ncol=3, byrow=TRUE)
        r8 <- matrix(c(-Inf,0,1, 0,1,0), ncol=3, byrow=TRUE)
        
        dens.pos <- reclassify((rterrenu+rveghaut+rvegras+rdur),r7)
        no.dens <- reclassify(dens.pos,r8)
                  
        ddur <- ddur*no.dens + (ddur+i)*dens.pos 
        
        # repartition equitable entre les disponibiliter des 3 variables
        # reclassifier les Inf (pas de division zero possible)
        rinc <- matrix(c(-Inf,0,0, 0,(i/3),(i/3), (i/3),(i/2),(i/2), 
                        (i/2),i,i, i,Inf,0), ncol=3, byrow=TRUE) 
        
        inc <- reclassify((i/(rterrenu + rveghaut + rvegras)),rinc)
        
	      ## Mod possible only where variable present suffisament
        mod.dterrenu <- mod.dterrenu*no.dens + (mod.dterrenu-inc*rterrenu)*dens.pos
        
        mod.veghaut <- mod.veghaut*no.dens + (mod.veghaut-inc*rveghaut)*dens.pos
        
        mod.dvegras <- mod.dvegras*no.dens + (mod.dvegras-inc*rvegras)*dens.pos
        
        ## reclassfier pour valeur plus grande  que 1 dans ddur 
        ## et plus petite que 0 dans autres variables
        r9 <- matrix(c(1,Inf,1), ncol=3, byrow=TRUE)
        r10 <- matrix(c(-Inf,0,0), ncol=3, byrow=TRUE)
        
        ddur <- reclassify(ddur,r9)
        
        mod.dterrenu <- reclassify(mod.dterrenu,r10)
        
        mod.veghaut <- reclassify(mod.veghaut,r10)
        
        mod.dvegras <- reclassify(mod.dvegras,r10)
        
        ## concatened mod avec zone unmodified              
        dterrenu2 <- unmod.dterrenu + mod.dterrenu
        
        dveghaut2 <- unmod.veghaut + mod.veghaut
        
        dvegras2 <- unmod.dvegras + mod.dvegras 
        
        ## write Output for densification scenario  
        writeRaster(ddur, filename=paste(outpath,"ddur.tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE)
        writeRaster(dterrenu2, filename=paste(outpath,"dterrenu.tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE)              
        writeRaster(dveghaut2, filename=paste(outpath,"dveghaut.tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE) 
        writeRaster(dvegras2, filename=paste(outpath,"dvegras.tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE) 
        
        # copy all other predictor
        l.pred <- c("dtree","fsumrad_2m",
                    "humcon","lengthfas","mfroadw" )
        
        for (j in 1:length(l.pred))
        { 
          file.copy(paste(inpath,"/",l.pred[j],".tif",sep=""), outpath )
        }
  
    }else{
      dterrenu2 <- dterrenu
      dveghaut2 <-dveghaut
      dvegras2 <- dvegras
  }
    
  ## Calculate reel augmentation
  mdur <- mean(as.array(ddur),na.rm=TRUE)
  mtree <- mean(as.array(dtree),na.rm=TRUE)
  mterrenu <-mean(as.array(dterrenu2),na.rm=TRUE)
  mveghaut <-mean(as.array(dveghaut2),na.rm=TRUE)
  mvegras <- mean(as.array(dvegras2),na.rm=TRUE)
        
  mat.land.mod[(n+1),] <- c(mdur,mtree,mterrenu,mvegras,mveghaut)     
          
  }

## Write table - reel augmentation
write.table(mat.land.mod,file=paste(working_dir,"pred/land.mod.txt",sep=""),
            sep="\t",append=F,row.names=T,col.names=T,quote=F)

#proc.time() - ptm