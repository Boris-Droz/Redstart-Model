####################################################################################################################################
###					
###  Calculate metric to evaluate the model prediction performance in a new range
### ============================================================================
######################################################################################################################################
####################################################################################################################################
## HISTORIC:
## --------
##   v1.0, August 2021 - Boris DROZ 
## 
## Description
##############
#   MESS index to define how similar point in space is to a 
##   reference set of point (Elith 2010) and Mateo 2015
#  Boyce index to evaluate the capacity to predict new range (Hirzel 2006)

## References
#############
### Elith, J.; Kearney, M.; Phillips, S., The art of modelling range-shifting
# species. Methods in Ecology and Evolution 2010, 1, (4), 330-342.
### Mateo, R. G.; Broennimann, O.; Petitpierre, B.; Munoz, J.; van Rooy, J.; 
# Laenen, B.; Guisan, A.; Vanderpoorten, A., What is the potential of spread 
# in invasive bryophytes? Ecography 2015, 38, (5), 480-487.
### Hirzel, A. H.; Le Lay, G.; Helfer, V.; Randin, C.; Guisan, A., Evaluating 
# the ability of habitat suitability models to predict species presences. 
# Ecol. Model. 2006, 199, (2), 142-152.

####################################################################################################################################
## PARAMETER
###############
# set working directory 
working_dir <-"~/Actual/rafb_proj_model/analysis/" 
# set path to function folder of one of your Rafb_model_project
pathtofunction_fold <- paste(working_dir,"Cdf_2021/script/function",sep="")

# set ouput file location
creat.subDir(working_dir,"/Cdf_2021/eval_output/new.rang_eval")
outpath <- paste(working_dir,"/Cdf_2021/eval_output/new.rang_eval/",sep="")
creat.subDir(working_dir,"/Cdf_2021/eval_output/new.rang_eval/MESS_raster")
outpath.mess <- paste(working_dir,"/Cdf_2021/eval_output/new.rang_eval/MESS_raster",sep="")

species <- "RaFB"

# see print(list.proj) to selection
#sel.mod <- c(4,5,6,8,9,10,11,12)
sel.mod <- 7

proj <- 2056 # input CRS code projection new swiss projection 21781 # old swiiss
dat.cal <- "Cdf_2021/input/dat.obs/terri_2021.txt"
dat.bgrd <- "Cdf_2021/input/dat.obs/bckgrd_2021.txt"

#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################}
setwd(working_dir)

# load all zone of mod
list.proj <- list.dirs(path=working_dir, recursive = FALSE)
list.name <- list.dirs(path=working_dir, full.names = FALSE, recursive = FALSE)
print(list.name)

source(paste(pathtofunction_fold,"/creat_subDirv1.R", sep="") )
source(paste(pathtofunction_fold,"/packagesv1.0.R", sep="") )

check.lib (c("raster","sp","magrittr","rgdal","rgeos","RNetCDF","dismo",
             "splitstackshape","maptools","XML","grainchanger","gstat", "sf",
             "ecospat"))

###################################################################################
file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

#####################################################################################################################################
## OPEN Data OBS territory 
dat.obs <- read.table(OBS.path, header=TRUE)

# calculate mediane of territory
##creat table median
x<-boxplot(as.numeric(dat.obs$X_COORD)~dat.obs$ID_UN)
y<-boxplot(as.numeric(dat.obs$Y_COORD)~dat.obs$ID_UN)

M<-data.frame(x$stats[3,],y$stats[3,],x$names)
names(M)<-c("X","Y","Terri")

# consider all area around each median point.
OBS.p <- M
coordinates(OBS.p) <- ~X+Y
proj4string(OBS.p)<-CRS("+init=epsg:2056") 

OBS.area <- buffer(OBS.p, width=100, dissolve=FALSE)
#####################################################################
## open data cal
# open data cal and background
d.cal <- read.table(dat.cal,header=TRUE)
d.bgrd <-read.table(dat.bgrd,header=TRUE)
refpt <- rbind(d.cal[,c(1:2,4:12)],d.bgrd)

######################################################################
## INDEX
#############
  out.mess <- NULL
  out.t.f <-NULL
  out.boy <-NULL
  
  for (n in sel.mod) # loop for zone
  {
    cat("> Zone",list.name[n],"is compute...", "\n",append = FALSE)
    
    # load zone data
    fns <- list.files(path=paste(list.proj[n],"/pred/pred_cur", sep=""),
                      pattern = ".tif$",full.names = TRUE)
    
    r.pred <- stack(fns)
    
    # open binary and probability ensemble + zone of pred
    ens <- raster(paste(list.proj[n],"/proj/proj_cur/ENSEMBLE_bin_RaFB_E4_cur.tif", sep=""))
    ens.p <- raster(paste(list.proj[n],"/proj/proj_cur/ENSEMBLE_prob_RaFB_E4_cur.tif", sep=""))
    zone <- raster(paste(list.proj[n],"/zone/zone.mod_LV95_2m.tif", sep="") )
     
    #################### 
    # MESS  Index
    ###################
    mess <- mess(x=r.pred, v=refpt[,3:11], full=FALSE)
    
    png(paste(outpath.mess,"/",list.name[n],"_MESS.Index.png",sep=""),width = 480, height = 480) 
    
      plot(mess,main=list.name[n])
    
    dev.off()
    
    writeRaster(mess, filename = paste(outpath.mess,"/",list.name[n],"_MESS.Index.tif",sep=""),
                datatype="FLT8S", overwrite=TRUE)
    
  #} # end loop to comput only messs
  #  cellStats(mess,mean)
    
    # eval for zone and for 0 to 1 pred
    # thres.mess <- c(">=0","0to-10","-10to-20","<20")
    thres.mess <- c(">=0","0to-20","<20")
    mess0 <- reclassify(mess, c(-Inf,0, NA, 0,Inf,1 ), include.lowest=TRUE)
    mess020 <- reclassify(mess, c(-Inf,-20, NA, -20,0,1, 0,Inf,NA ), include.lowest=TRUE)
    #mess010 <- reclassify(mess, c(-Inf,-10, NA, -10,0,1, 0,Inf,NA ), include.lowest=TRUE)
   # mess1020 <- reclassify(mess, c(-Inf,-20, NA, -20,-10,1, -10,Inf,NA ), include.lowest=TRUE)
    mess20 <- reclassify(mess, c(-Inf,-20, 1, -20,Inf,NA ), include.lowest=TRUE)
    
    #stack.mess <- stack(mess0,mess010,mess1020,mess20)
    stack.mess <- stack(mess0,mess020,mess20)
    
    # defined obs for the area if some obs exist
    t.obs.area <- crop(OBS.area,zone)
  if (is.null(t.obs.area)){ }else
    {
    	  t.obs.area@data <-cbind(t.obs.area@data, rep(1,time=nrow(t.obs.area@data) ))
    	  names(t.obs.area) <- c(names(OBS.area),"val")
    	  t.obs.area <- rasterize(t.obs.area,zone,field="val",background=0)

   	    # extract xy values for obs cell
    	  ext <- raster::extract(zone, OBS.area, cellnumbers=TRUE)
    	  v <- unlist(ext)# get all cell number
    	  xy.obs <- xyFromCell(zone, v) # get coordinate for all obs
   
    #################
    ## Boyce index
    #################
    # for all data obs get boyce index
    cat("> .....boyce index pred ensemble 0","\n",append = FALSE)
    boyce.index <- ecospat.boyce (fit= ens.p,obs= xy.obs)
    d.boyce <-boyce.index$Spearman.cor  
    
    # for selected mess index category
    for (k in 1:3)
      {
      cat("> .....boyce index pred ensemble",k, "\n",append = FALSE)
        t.ens.p <- ens.p*stack.mess[[k]]
          
        # extract xy values for obs cell
        ext <- raster::extract(zone, OBS.area, cellnumbers=TRUE)
        v <- unlist(ext)# get all cell number
        xy.obs <- xyFromCell(zone, v) # get coordinate for all obs
        
        if (cellStats(t.ens.p,sum)==0) {
          d.boyce <-c(d.boyce, 0)
          }else{
          boyce.index <-ecospat.boyce (fit= t.ens.p,obs= xy.obs)
          d.boyce <-c(d.boyce, boyce.index$Spearman.cor)
          }
      } # end loop cat 1-4
    } # end loop if data
      
    dat <- NULL
    d.tf <- NULL
      for (l in 0:4) # loop ensemble
          {
          cat("> .....MESS index pred ensemble",l, "\n",append = FALSE)
        
        # reclassified zone
          ens.temp <- reclassify(ens, c(-Inf,l-1, NA, 
                                    l-1,l,1,
                                    l,Inf,NA ))
          
          if (is.null(t.obs.area)){ 
            ###############
            # TRUE-FALSE
            ###############
            if (l==0){# true neg
                n.t.neg <- cellStats(ens.temp,stat='sum')
              
                # false neg
                n.f.neg <- 0
              
              }else{# true pos
                n.t.pos <- 0
              
                # true neg
                n.f.pos <- cellStats(ens.temp,stat='sum')
              
                d.tf <-rbind(d.tf,data.frame(n.t.neg=n.t.neg,n.f.neg=n.f.neg,
                                           n.t.pos=n.t.pos,n.f.pos=n.f.pos))
              }
           }else{
          ###############
          # MESS INDEX
          ###############
          # calculus mess for obs data
          t.obs <- ens.temp*stack.mess *zone*t.obs.area
          
          nb.obs <- cellStats(t.obs,stat='sum')
          p.obs <- nb.obs/sum(nb.obs)
          
          # calculus mess for entire zone pred
          t.zpred <- ens.temp*stack.mess *zone
          
          nb.zpred <- cellStats(t.zpred,stat='sum')
          p.zpred <- nb.zpred/sum(nb.zpred)
            
          dat <-cbind(dat,data.frame=cbind(nb.obs, p.obs, nb.zpred,p.zpred))
          
          n.dat <- paste(c("n.obs_","p.obs_","n.zpred_","p.zpred_"),rep(l,4), sep="")
          
          pos <- (length(colnames(dat))-3):length(colnames(dat))
          colnames(dat)[pos] <- n.dat
          
          ###############
          # TRUE-FALSE
          ###############
          if (l==0){# true neg
                  t.neg <- ens.temp- (ens.temp * t.obs.area)
                  n.t.neg <- cellStats(t.neg,stat='sum')
                  
                  # false neg
                  f.neg <- ens.temp * t.obs.area
                  n.f.neg <- cellStats(f.neg,stat='sum')
                  
              }else{# true pos
                   t.pos <- ens.temp * t.obs.area
                   n.t.pos <- cellStats(t.pos,stat='sum')
                   
                   # true neg
                  f.pos <- ens.temp -(ens.temp * t.obs.area)
                  n.f.pos <- cellStats(f.pos,stat='sum')
                  
                  d.tf <-rbind(d.tf,data.frame(n.t.neg=n.t.neg,n.f.neg=n.f.neg,
                                               n.t.pos=n.t.pos,n.f.pos=n.f.pos))
              }
          } # end loop if data
          
      } # end loop ensemble
    
    sum.d.tf <- sum(d.tf[1,1],d.tf[1,2],d.tf[,3],d.tf[,4])
    p.d.tf <- d.tf/sum.d.tf
    
    colnames(p.d.tf) <-c("p.t.neg", "p.f.neg", "p.t.pos", "p.f.pos")
    
    d.tf <- cbind(d.tf,p.d.tf)
     
    if (is.null(t.obs.area)){ }else{       
      out.mess <- rbind(out.mess, data.frame(zone=rep(list.name[n],3), 
                            mess= paste(rep("mess_",3),thres.mess,sep=""),
                            dat))
    
      out.boy <- rbind(out.boy, data.frame(zone=list.name[n],t(d.boyce) ) )
    }
    
    out.t.f <- rbind(out.t.f, data.frame(zone= rep(list.name[n],4),
                                         ens= paste(rep("ens_",4),1:4,sep=""),
                                         d.tf))
      
          }# end loop for zone
  
  colnames(out.boy) <- c("zone","all.obs",paste(rep("mess_",3),thres.mess,sep="") )
      
  # write dataset-----
  write.table(out.mess, file=paste(outpath,"/mess_index.txt",sep=""),sep="\t",
              append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE) 
  
  write.table(out.t.f, file=paste(outpath,"/true-false.txt",sep=""),sep="\t",
              append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
  
  write.table(out.boy, file=paste(outpath,"/boyce_index.txt",sep=""),sep="\t",
              append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
  
  
  