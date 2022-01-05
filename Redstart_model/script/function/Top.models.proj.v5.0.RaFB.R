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
## optimized to extent prediction projection and for parrallel work                                                                         #
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

# function and parameters start here
top.models.proj <-function(inpath,
                           GLM.mod=TRUE,GLM.ModAv=TRUE,GLM.Step=TRUE,
                           GAM.Step=TRUE,GBM.mod=TRUE,ME.mod=TRUE,
                           species=c("RaFB"),
                           envstack=env.stack.proj,
                           outpath.proj,
                           scenario="Cur")
  # function starts here
{
  ##############
  # species loop
  ##############
  
  for(i in 1:length(species))
  {
    cat(paste('\n\n\n\n################################\n### ',species[i],'\n################################',sep=''),"\n")
    ##########################################     		      
    #-------------------
    # predict GLM models
    #-------------------
    if(GLM.mod==TRUE)
    {
      #####################################################
      if (GLM.ModAv==TRUE)
      {
        
        cat('> GLM projections ',"\n")
        #library needed
        require(raster) 
        require(rms) #instead of Design
        require(gam)
        require(gbm)
        require(rJava)
        require(dismo)
        require(raster)
        require(rgdal)
        require(maptools)
        
        ## laod data
        load(paste(inpath,species[i],'_glm_modav_model',sep=''))
        
        form <- glm.modav.fmula.proj # text formula
        
        eval(parse(text = paste("custom.overlay <- function(",
                                paste(names(envstack), collapse = ", "),
                                ") {",
                                form,
                                "}" ))) # eval text to formula
        
        glm.modav.prob.proj <- overlay(envstack, fun = custom.overlay)
        
        #original     ### 
        ## glm.modav.prob.proj <- eval(parse(text = glm.modav.fmula.proj) )
        
        m <- as.numeric (c(0, glm.modav.max.tss.tresh, 0, glm.modav.max.tss.tresh, 1, 1))
        rclmat.glm.modav <- matrix(m, ncol=3, byrow=TRUE)
        glm.modav.bin.proj <- reclassify(glm.modav.prob.proj, rclmat.glm.modav)
        
        #	writeRaster(glm.modav.prob.proj, filename=paste(outpath.proj,"proj_",species[i],"_glm_prob_modav_",scenario,".tif",sep=""), 
        #		datatype="FLT8S", overwrite=TRUE)		
        
        writeRaster(glm.modav.bin.proj, filename=paste(outpath.proj,"proj_",species[i],"_glm_bin_modav_",scenario,".tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE)		
        
        rm(glm.modav.prob.proj,glm.modav.bin.proj)	
        
      }
      
      if (GLM.Step==TRUE)
      {
        cat('> GLM projections ',"\n")
        
        library(raster) 
        ## laod data
        load(paste(inpath,species[i],'_glm_stepwise_model',sep=''))
        
        glm.stepwise.prob.proj<-predict(envstack,glm.tmp.step,type='response')
        
        m <- as.numeric(c(0, glm.stepwise.max.tss.tresh, 0, glm.stepwise.max.tss.tresh, 1, 1))
        rclmat.glm.stepwise <- matrix(m, ncol=3, byrow=TRUE)
        glm.stepwise.bin.proj <- reclassify(glm.stepwise.prob.proj, rclmat.glm.stepwise)
        
        #	writeRaster(glm.stepwise.prob.proj, filename=paste(outpath.proj,"proj_",species[i],"_glm_prob_stepwise_",scenario,".tif",sep=""), 
        #		datatype="FLT8S", overwrite=TRUE)		
        
        writeRaster(glm.stepwise.bin.proj, filename=paste(outpath.proj,"proj_",species[i],"_glm_bin_stepwise_",scenario,".tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE)		
        
        rm(glm.stepwise.prob.proj,glm.stepwise.bin.proj)	
        
      }
    }
    
    ##########################################         
    #-----------------
    # predict GAM models
    #-------------------
    if (GAM.Step==TRUE)
    {
      cat('> GAM projections ',"\n")
      
      library(raster) 
      ## laod data
      load(paste(inpath,species[i],'_gam_stepwise_model',sep=''))
      
      gam.stepwise.prob.proj<-predict(envstack,gam.tmp.step,type='response')
      
      m <- as.numeric(c(0, gam.stepwise.max.tss.tresh, 0, gam.stepwise.max.tss.tresh, 1, 1))
      rclmat.gam.stepwise <- matrix(m, ncol=3, byrow=TRUE)
      gam.stepwise.bin.proj <- reclassify(gam.stepwise.prob.proj, rclmat.gam.stepwise)
      
      #	writeRaster(gam.stepwise.prob.proj, filename=paste(outpath.proj,"proj_",species[i],"_gam_prob_stepwise_",scenario,".tif",sep=""), 
      #		datatype="FLT8S", overwrite=TRUE)		
      
      writeRaster(gam.stepwise.bin.proj, filename=paste(outpath.proj,"proj_",species[i],"_gam_bin_stepwise_",scenario,".tif",sep=""), 
                  datatype="FLT8S", overwrite=TRUE)		
      
      rm(gam.stepwise.prob.proj,gam.stepwise.bin.proj)	
      
    }
    #-----------------
    # predict GBM models
    #-------------------
    
    if(GBM.mod==TRUE)
    {            
      library(raster)
      cat('> GBM projections',"\n")
      
      #cmd <- paste("ls ", paste(inpath,sep=''),paste("*",species[i],"_gbm_model_rep*",sep=''), sep="")
      #list.rep.gbm <- system(cmd, intern=TRUE)
      
      list.rep.gbm <- list.files(inpath,
                                 pattern=paste(species[i],"_gbm_model_rep",sep=""),
                                 full.names = TRUE)
      
      for (r in 1:length(list.rep.gbm))
      {
        cat('--> GBM rep ',r,"\n")
        
        load(paste(list.rep.gbm[r],sep=''))
        
        gbm.prob.proj <- predict(envstack,gbm.tmp.cal,type='response',n.trees=tmp.best.itr.gbm)
        m <- as.numeric(c(0, gbm.max.tss.tresh, 0, gbm.max.tss.tresh, 1, 1))
        rclmat.gbm <- matrix(m, ncol=3, byrow=TRUE)
        gbm.bin.proj <- reclassify(gbm.prob.proj, rclmat.gbm)
        
        #writeRaster(gbm.prob.proj, filename=paste(outpath.proj,"proj_",species[i],"_gbm_prob_rep_",r,"_",scenario,".asc",sep=""), 
        #           datatype="FLT8S", overwrite=TRUE)		
        
        writeRaster(gbm.bin.proj, filename=paste(outpath.proj,"proj_",species[i],"_gbm_bin_rep_",r,"_",scenario,".tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE)		
        
        rm(gbm.prob.proj,gbm.bin.proj)	  	 				
      }
      
    }
    #---------------------
    # predict maxent models
    #----------------------
    if(ME.mod==TRUE)
    {      
      
      cat('> MaxEnt projections ',"\n")
      
      #cmd <- paste("ls ", paste(inpath,sep=''),paste("*",species[i],"_me_model_rep*",sep=''), sep="")
      #list.rep.me <- system(cmd, intern=TRUE)
      list.rep.me <- list.files(inpath,
                                pattern=paste(species[i],"_me_model_rep",sep=""),
                                full.names = TRUE)
      
      for (r in 1:length(list.rep.me))
      {
        cat('--> MaxEnt rep ',r,"\n")
        
        load(paste(list.rep.me[r],sep=''))
        
        #jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
        #Sys.setenv(NOAWT=TRUE)
        #library(dismo)
        # maxent()
        me.prob.proj <- predict(me,envstack,progress='text',args=c('logistic'))
        
        m <- as.numeric(c(0, me.max.tss.tresh, 0, me.max.tss.tresh, 1, 1))
        rclmat.me <- matrix(m, ncol=3, byrow=TRUE)
        me.bin.proj <- reclassify(me.prob.proj, rclmat.me)
        
        #		writeRaster(me.prob.proj, filename=paste(outpath.proj,"proj_",species[i],"_me_prob_rep_",r,"_",scenario,".tif",sep=""), 
        #		datatype="FLT8S", overwrite=TRUE)		
        
        writeRaster(me.bin.proj, filename=paste(outpath.proj,"proj_",species[i],"_me_bin_rep_",r,"_",scenario,".tif",sep=""), 
                    datatype="FLT8S", overwrite=TRUE)		
        
        rm(me.prob.proj,me.bin.proj)	  	 				
      }
      
    }
    
  }# End of loop species   
} # End of loop function