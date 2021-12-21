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
## optimized to extent prediction projection and for parrallel work 
#                                                                             #
###############################################################################
####################################################################################################################################
## DESCRIPTION
###----------
## VERSION FOR PC running under R 3.5.0
## NEED FOR RUNNING:
##------------------
## Install Java from here http://www.filehippo.com/download_jre_64/
## work with Java Runtime Environment 64-bit 8.0-build-271
##
## Install Maxent v3.3.3k from https://www.cs.princeton.edu/~schapire/maxent/
## dezip the file in your personnal package library
## like "R/win-library/3.5/dismo/java/maxent.jar"
###################################################################################
# Function starts here

ensemble.spat.proj <- function(inpath=inpath,
                               species="RaFB",scenario="Cur",
                               glm.modav.mod=TRUE,n.rep.glm.modav=1,
                               glm.stepwise.mod=TRUE,n.rep.glm.stepwise=1,	
                               gam.stepwise.mod=TRUE,n.rep.gam.stepwise=1,
                               gbm.mod=TRUE,n.rep.gbm=10,
                               rf.mod=TRUE,n.rep.rf=10,
                               me.mod=TRUE,n.rep.me=10,
                               ens.code="E4",
                               outpath=outpath)
{
  cat("Starting ENSEMBLE projections for scenario ",scenario," and type ",ens.code,"...", "\n",append = F)
  
  # GBM
  if (gbm.mod==TRUE)
  {
    
    cat("-> Processing GBM projections...", "\n",append = FALSE)
    
    for (i in 1:n.rep.gbm)
    {
      cat("--> Replication ",i, "\n",append = FALSE)
      
      if (i==1 && n.rep.gbm==1)
      {
        
      } else {
        
        assign(paste("gbm.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_gbm_bin_rep_",i,"_",scenario,".tif",sep=""))) 
      }
      
      if (i ==1)
      {
        sum.gbm.bin.proj <- get(paste("gbm.bin.proj.",i,sep=""))
      } else {
        sum.gbm.bin.proj <- sum.gbm.bin.proj + get(paste("gbm.bin.proj.",i,sep=""))
      }
      
      # rm(get(paste("gbm.bin.proj.",i,sep="")))		
    }
    
    n.maj.gbm <- ceiling(n.rep.gbm / 2)
    sum.gbm.bin.proj[sum.gbm.bin.proj<n.maj.gbm] <- 0
    sum.gbm.bin.proj[sum.gbm.bin.proj>=n.maj.gbm] <- 1
    
    
    ensemble.bin.all <- sum.gbm.bin.proj
    
    # End of loop GBM					
  }
  ###################################################################################################################################################
  # ME
  if (me.mod==TRUE)
  {
    
    cat("-> Processing ME projections...", "\n",append = FALSE)
    
    for (i in 1:n.rep.me)
    {
      cat("--> Replication ",i, "\n",append = FALSE)
      
      if (i==1 && n.rep.gbm==1)
      {
        
      } else {
        
        assign(paste("me.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_me_bin_rep_",i,"_",scenario,".tif",sep=""))) 
      }
      
      if (i ==1)
      {
        sum.me.bin.proj <- get(paste("me.bin.proj.",i,sep=""))
      } else {
        sum.me.bin.proj <- sum.me.bin.proj + get(paste("me.bin.proj.",i,sep=""))
      }
      
      # rm(get(paste("me.bin.proj.",i,sep="")))		
    }
    
    n.maj.me <- ceiling(n.rep.me / 2)
    sum.me.bin.proj[sum.me.bin.proj<n.maj.me] <- 0
    sum.me.bin.proj[sum.me.bin.proj>=n.maj.me] <- 1
    
    if (gbm.mod==FALSE)
    {
      ensemble.bin.all <- sum.me.bin.proj
      
    } else	
    {
      ensemble.bin.all <- ensemble.bin.all + sum.me.bin.proj
    }
    # End of loop ME			
  }
  ###################################################################################################################################################
  # RF
  if (rf.mod==TRUE)
  {
    
    cat("-> Processing RF projections...", "\n",append = FALSE)
    
    for (i in 1:n.rep.rf)
    {
      cat("--> Replication ",i, "\n",append = FALSE)
      
      if (i==1 && n.rep.gbm==1)
      {
        
      } else {
        
        assign(paste("rf.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_rf_bin_rep_",i,"_",scenario,".tif",sep=""))) 
      }
      
      if (i ==1)
      {
        sum.rf.bin.proj <- get(paste("rf.bin.proj.",i,sep=""))
      } else {
        sum.rf.bin.proj <- sum.rf.bin.proj + get(paste("rf.bin.proj.",i,sep=""))
      }
      
      # rm(get(paste("rf.bin.proj.",i,sep="")))		
    }
    
    n.maj.rf <- ceiling(n.rep.rf / 2)
    sum.rf.bin.proj[sum.rf.bin.proj<n.maj.rf] <- 0
    sum.rf.bin.proj[sum.rf.bin.proj>=n.maj.rf] <- 1
    
    if (gbm.mod==FALSE && me.mod==FALSE)
    {
      ensemble.bin.all <-  sum.rf.bin.proj
      
    } else	
    {
      ensemble.bin.all <- ensemble.bin.all + sum.rf.bin.proj
    }
    
    
    # End of loop RF			
  }
  ###################################################################################################################################################
  # GLM ModAv
  if (glm.modav.mod==TRUE)
  {
    
    cat("-> Processing GLM ModAv projections...", "\n",append = FALSE)
    
    for (i in 1:n.rep.glm.modav)
    {
      cat("--> Replication ",i, "\n",append = FALSE)
      
      if (i==1 && n.rep.glm.modav==1)
      {
        assign(paste("glm.modav.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_glm_bin_modav_",scenario,".tif",sep=""))) 
      } else {
        
        assign(paste("glm.modav.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_glm_bin_modav_rep_",i,"_",scenario,".tif",sep=""))) 
      }
      
      if (i ==1)
      {
        sum.glm.modav.bin.proj <- get(paste("glm.modav.bin.proj.",i,sep=""))
      } else {
        sum.glm.modav.bin.proj <- sum.glm.modav.bin.proj + get(paste("glm.modav.bin.proj.",i,sep=""))
      }
      
      # rm(get(paste("glm.modav.bin.proj.",i,sep="")))		
    }
    
    n.maj.glm.modav <- ceiling(n.rep.glm.modav / 2)
    sum.glm.modav.bin.proj[sum.glm.modav.bin.proj<n.maj.glm.modav] <- 0
    sum.glm.modav.bin.proj[sum.glm.modav.bin.proj>=n.maj.glm.modav] <- 1
    
    if (gbm.mod==FALSE && me.mod==FALSE && rf.mod==FALSE)
    {
      ensemble.bin.all <-  sum.glm.modav.bin.proj
      
    } else	
    {
      ensemble.bin.all <- ensemble.bin.all + sum.glm.modav.bin.proj
    }
    # End of loop GLM ModAv		
  }
  ###################################################################################################################################################
  # GLM Stepwise
  if (glm.stepwise.mod==TRUE)
  {
    
    cat("-> Processing GLM Stepwise projections...", "\n",append = FALSE)
    
    for (i in 1:n.rep.glm.stepwise)
    {
      cat("--> Replication ",i, "\n",append = FALSE)
      
      if (i==1 && n.rep.glm.stepwise==1)
      {
        assign(paste("glm.stepwise.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_glm_bin_stepwise_",scenario,".tif",sep=""))) 
      } else {
        
        assign(paste("glm.stepwise.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_glm_bin_stepwise_rep_",i,"_",scenario,".tif",sep=""))) 
      }
      
      if (i ==1)
      {
        sum.glm.stepwise.bin.proj <- get(paste("glm.stepwise.bin.proj.",i,sep=""))
      } else {
        sum.glm.stepwise.bin.proj <- sum.glm.stepwise.bin.proj + get(paste("glm.stepwise.bin.proj.",i,sep=""))
      }
      
      # rm(get(paste("glm.stepwise.bin.proj.",i,sep="")))		
    }
    
    n.maj.glm.stepwise <- ceiling(n.rep.glm.stepwise / 2)
    sum.glm.stepwise.bin.proj[sum.glm.stepwise.bin.proj<n.maj.glm.stepwise] <- 0
    sum.glm.stepwise.bin.proj[sum.glm.stepwise.bin.proj>=n.maj.glm.stepwise] <- 1
    
    if (gbm.mod==FALSE && me.mod==FALSE && rf.mod==FALSE && glm.modav.mod==TRUE)
    {
      ensemble.bin.all <-  sum.glm.stepwise.bin.proj
      
    } else	
    {
      ensemble.bin.all <- ensemble.bin.all + sum.glm.stepwise.bin.proj
    }
    #End of loop GLM Stepwise		
  }
  ###################################################################################################################################################
  # GAM Stepwise
  if (gam.stepwise.mod==TRUE)
  {
    
    cat("-> Processing GAM Stepwise projections...", "\n",append = FALSE)
    
    for (i in 1:n.rep.gam.stepwise)
    {
      cat("--> Replication ",i, "\n",append = FALSE)
      
      if (i==1 && n.rep.gam.stepwise==1)
      {
        assign(paste("gam.stepwise.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_gam_bin_stepwise_",scenario,".tif",sep=""))) 
      } else {
        
        assign(paste("gam.stepwise.bin.proj.",i,sep=""),
               raster(paste(inpath,"proj_RaFB_gam_bin_stepwise_rep_",i,"_",scenario,".tif",sep=""))) 
      }
      
      if (i ==1)
      {
        sum.gam.stepwise.bin.proj <- get(paste("gam.stepwise.bin.proj.",i,sep=""))
      } else {
        sum.gam.stepwise.bin.proj <- sum.gam.stepwise.bin.proj + get(paste("gam.stepwise.bin.proj.",i,sep=""))
      }
      
      # rm(get(paste("gam.stepwise.bin.proj.",i,sep="")))		
    }
    
    n.maj.gam.stepwise <- ceiling(n.rep.gam.stepwise / 2)
    sum.gam.stepwise.bin.proj[sum.gam.stepwise.bin.proj<n.maj.gam.stepwise] <- 0
    sum.gam.stepwise.bin.proj[sum.gam.stepwise.bin.proj>=n.maj.gam.stepwise] <- 1
    
    if (gbm.mod==FALSE && me.mod==FALSE && rf.mod==FALSE && glm.modav.mod==TRUE && glm.stepwise.mod==FALSE)
    {
      ensemble.bin.all <-  sum.gam.stepwise.bin.proj
      
    } else	
    {
      ensemble.bin.all <- ensemble.bin.all + sum.gam.stepwise.bin.proj
    }
    # End of loop GAM stepwise		
  }
  
  v.unique <- unique(ensemble.bin.all)
  max.ens <- max(as.vector(na.omit(v.unique)))		
  
  writeRaster(ensemble.bin.all,paste(outpath,"ENSEMBLE_bin_",species,"_",ens.code,"_",scenario,".tif",sep=""),
              datatype="FLT8S", overwrite=TRUE)
  
  ens.types <- max.ens
  
  #	while(ens.types > 0)
  #	{
  #		tmp.ensemble <- ensemble.bin.all
  #		tmp.ensemble[tmp.ensemble < ens.types] <- 0
  #		tmp.ensemble[tmp.ensemble >= ens.types] <- 1
  #		
  #		cat("--> Exporing ENSEMBLE ",ens.types, "\n",append = F)
  #			
  #		write.asc(tmp.ensemble,paste(outpath,"ENSEMBLE_",ens.types,"_",species,"_",ens.code,"_",scenario,".asc",sep=""))
  #        
  #        ens.types <- (ens.types-1)
  #        rm(tmp.ensemble)
  #		
  #	# End of loop ENSEMBLE
  #	}	
  
  # End of function	
}	
