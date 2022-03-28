####################################################################################################################################
###					
###  EVALUATION OF THE FALSE NEG DATA POINT
### =============================================
######################################################################################################################################
####################################################################################################################################
## HISTORIC:
## --------
##   v1.0, August 2021 - Boris DROZ 
## 
## Description
##############
### Evaluate for the false negative area which predicting variables is not satisfied.
## Test all conditon until the 4 model technic are satisfied (optimum territory)
####################################################################################################################################
## PARAMETER
###############
# set working directory
#######################
working_dir <-"~/Actual/rafb_proj_model/analysis/"
# set path to function folder of one of your Rafb_model_project
pathtofunction_fold <- paste(working_dir,"Cdf_2021/script/function",sep="")
# set data obs location
OBS.path <- "~/Actual/rafb_proj_model/analysis/Cdf_2021/input/dat.obs/rafb_med.obs_allzone.txt"
# set ouput file location
###########################
# folder for temporary prediction
creat.subDir(working_dir,"/Cdf_2021/eval_output/new.rang_eval/f.sen.data")
outpath <- paste(working_dir,"/Cdf_2021/eval_output/new.rang_eval/f.sen.data/",sep="")
# folder for final result
creat.subDir(working_dir,"/Cdf_2021/eval_output/new.rang_eval/sen.raster")
outpath.sens <- paste(working_dir,"/Cdf_2021/eval_output/new.rang_eval/sen.raster/",sep="")

# optimum pred var file table
f.opt <-"~/Actual/rafb_proj_model/analysis/Cdf_2021/eval_output/proj_model_2021_eval/Optim_RaFB.txt"

#MODEL PARAMETER
################
species=c("RaFB") # model name
m.tech <- c("GLM.ModAv","GLM.Step","GAM.Step","ME") # model technique
nb.step <- 100

# see print(list.proj --> line 61) to selection -> order of all projects
sel.mod <- c(8,10,11,12)

proj <- 2056 # input projection new swiss projection 21781 # old swiiss

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
source(paste(pathtofunction_fold,"/packagesv1.0.R", sep=""))
source(paste(pathtofunction_fold,"/Top.models.proj.v5.0.RaFB.R", sep=""))
source(paste(pathtofunction_fold,"/ensemble_projv5.0.Rafb.R", sep=""))
##################
## Load R. library
##################
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
Sys.setenv(NOAWT=TRUE)

check.lib (c("rms","gam","gbm","rJava","dismo","raster","rgdal","maptools") )
maxent()

check.lib (c("raster","sp","magrittr","rgdal","rgeos","RNetCDF","dismo",
             "splitstackshape","maptools","XML","grainchanger","gstat", "sf",
             "ecospat"))

###################################################################################
file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 
#####################################################################################################################################
## OPEN Data 
################
# path.dat <- "~/Actual/rafb_proj_model/analysis/Cdf_2021/input/dat.obs/rafb_med.obs_allzone.txt"
OBS.p <- read.table(OBS.path, header=TRUE)

# optimum pred var
p.v.opt <- read.table(f.opt,header=TRUE)
v.opt <- apply(p.v.opt[c(1:3,5),],2,mean)

# consider all area around each median point.
coordinates(OBS.p) <- ~X+Y
proj4string(OBS.p)<-CRS("+init=epsg:2056") 

OBS.area <- buffer(OBS.p, width=100, dissolve=FALSE)

###################################################################
  for (n in sel.mod) # loop for zone
  {
    cat("> Zone",list.name[n],"is compute...", "\n",append = FALSE)
    cat("> ###################################", "\n",append = FALSE)
    
    # load zone data
    fns <- list.files(path=paste(list.proj[n],"/pred/pred_cur", sep=""),
                      pattern = ".tif$",full.names = TRUE)
    
    r.pred <- stack(fns)
    
    ens <- raster(paste(list.proj[n],"/proj/proj_cur/ENSEMBLE_bin_RaFB_E4_cur.tif", sep=""))
    zone <- raster(paste(list.proj[n],"/zone/zone.mod_LV95_2m.tif", sep="") )
    
    # ignite variable
    # defined obs for the area
    t.obs.area <- OBS.area
    t.obs.area@data <-cbind(t.obs.area@data, rep(1,time=nrow(t.obs.area@data) ))
    names(t.obs.area) <- c(names(OBS.area),"val")
    r.obs.area <- rasterize(t.obs.area,zone,field="val",background=NA)
    # reclassified zone
    #ens.temp <- r.obs.area *ens

    # if need to restart do for current 
    ## ens.area <- raster(file.choose() ) ## all.neg.sen...
	## ens.temp <- raster(file.choose() ) ## temporary ensemble bin

    ens.area <- reclassify(ens, c(-Inf,-1, NA, # area final of meth
                                       -1,3,0,
                                       3,Inf,NA))
    
    ens.temp <- reclassify(ens, c(-Inf,-1, NA, 
                                       -1,3,1,
                                       3,Inf,NA))
    i=1
    opt.meth <-c(1:9,12:19,23:29,34:39,45:49,56:59,67:69,78:79,89,
			123:129,134:139,145:149,156:159,167:169,178:179,189,
			234:239,245:249,256:259,267:269,278:279,289,
			345:349,356:359,356:359,367:369,378:379,389,
			456:459,467:469,478:479,489,
			567:569,578:579,589, 678:679,689, 789,
			1234:1239,1245:1249,1256:1259,1267:1269,1278:1279,1289,
			1345:1349,1356:1359,1367:1368,1378:1379,1389,
			1456:1459,1467:1469,1478:1479,1489,
			1567:1569,1578:1579,1589, 1678:1679,1689, 1789,
			2345:2349,2356:2359,2367:2369,2378:2379,2389,
			2456:2459,2467:2469,2478:2479,2489,
			2567:2569,2678:2679,2789,
			3456:3459,3467:3469,3478:3479,3489,
			3567:3569,3578:3579,3589, 3678:3679,3689, 3789,
			4567:4569,4578:4579,4589, 4678:4679,4789,
			5678:5679,5789)
    
    split.num <- function(v.num,tab.order){
              text <- as.character(v.num)
              vec.text <-NULL
              for (z in 1:nchar(text)) {
                vec.text <-c(vec.text, substring(text, z, z) )}
                vec.text <- as.numeric(vec.text)
                rank <- tab.order[vec.text]
                return(rank)
                }
    d.rank <- c(3, 5, 1, 4, 8, 7, 9, 2, 6) # to convert alpha order to imp order
    r.opt.rank <- lapply(opt.meth,split.num, d.rank) #[,2])
                        
    while (cellStats(ens.temp,sum)>0 & i<=length(opt.meth)) ## loop
        {
        cat("> ###################################", "\n",append = FALSE)
        cat(">... optimized",names(r.pred)[ r.opt.rank[[i]] ],"is compute...", "\n",append = FALSE)
        cat("> ###################################", "\n",append = FALSE)
        
        t.pred <- r.pred*ens.temp # extract pred
        t.pred.i <- t.pred
        
        pos <- unlist(r.opt.rank[i])
        
        for (k in 1:length(pos)){
        t.pred[[pos[k]]] <- reclassify(t.pred[[pos[k]]],
                                    c(-Inf,Inf, v.opt[pos[k]] ) )   # optimized desired data
        if (k==1){
            plu_min <- t.pred.i[[pos[k] ]]-t.pred[[pos[k] ]]
            plu_min <- reclassify(plu_min,c(-Inf,0,-1, 0,Inf, 1))                                  
          
          }else{
            pm <- t.pred.i[[pos[k] ]]-t.pred[[pos[k] ]]
            pm <- reclassify(pm,c(-Inf,0,-1, 0,Inf, 1))
            plu_min <- plu_min *pm
          }
        }
        
        plu_min [is.na(plu_min)] <- 0 # for all area
        
        names(t.pred)<-names(r.pred)# rename pred
        
        # set location of calibrate model and projection
        ################################################
        inpath <- paste(list.proj[n],"/script/geographic_proj_model/",sep="")
        envstack <- t.pred
        proj.name <- "sen_anal"
        outpath.proj <- outpath
        
        ##############################
        ########### MODEL RUN  #######
        ##############################
        # Run code
        top.models.proj(inpath,GLM.mod=TRUE,GLM.ModAv=TRUE,
                        GLM.Step=TRUE,GAM.Step=TRUE,GBM.mod=FALSE, 
                        ME.mod=TRUE,
                        species=c(species),envstack=t.pred,
                        outpath.proj,scenario=proj.name)
        
        ensemble.spat.proj(inpath=outpath.proj,
                           species=species,scenario=proj.name,
                           glm.modav.mod=TRUE,n.rep.glm.modav=1,
                           glm.stepwise.mod=TRUE,n.rep.glm.stepwise=1,	
                           gam.stepwise.mod=TRUE,n.rep.gam.stepwise=1,
                           gbm.mod=FALSE,n.rep.gbm=10,
                           rf.mod=FALSE,n.rep.rf=10,
                           me.mod=TRUE,n.rep.me=10,
                           ens.code="E4",
                           outpath=outpath.proj)

        # reload ensemble 
        ens.temp <- raster(paste(outpath.proj,"ENSEMBLE_bin_",species,"_E4_sen_anal.tif",sep=""))
        
        ens.opt <- reclassify(ens.temp, c(-Inf,3, 0, 
                                          3,4,1, 4,Inf,0 ))
        ens.opt [is.na(ens.opt)] <- 0
        
        ens.temp <- reclassify(ens.temp, c(-Inf,-1, NA, 
                          -1,3,1,
                          3,Inf,NA ))
        
        ens.area <- ens.area + ens.opt * opt.meth[i] * plu_min
        
        # reactualized for next loop
        i=i+1
        
    } # end while loop
    
    if (i==1){}else{
      
      names(ens.area) <- "all.neg.sen"
      
      writeRaster(ens.area, filename = paste(outpath.sens,"all.neg.sen",list.name[n],".tif",sep=""),
                  datatype="FLT8S", overwrite=TRUE)
      
      ## ens.area <- raster(file.choose()) ## 
      
      ens.neg <- r.obs.area *ens.area
    
      names(ens.neg) <- "f.neg.sen"
    
      writeRaster(ens.neg, filename = paste(outpath.sens,"f.neg.sen",list.name[n],".tif",sep=""),
                  datatype="FLT8S", overwrite=TRUE)
    }
    
  # make table of unique area
  ens.area <- abs(ens.area) 
  ens.neg <-abs(ens.neg)
  
  pos <- unique(ens.area)
  
  d.out <-NULL
  for (q in 1:length(pos))
      {
      pred.area <- reclassify(ens.area,c(-Inf,pos[q]-1, 0, 
                                         pos[q]-1,pos[q],1,
                                         pos[q],Inf,0 ))
      pred.area <- cellStats(pred.area, sum)
      
      neg.area <- reclassify(ens.neg,c(-Inf,pos[q]-1, 0, 
                                     pos[q]-1,pos[q],1,
                                     pos[q],Inf,0 ))
      neg.area <- cellStats(neg.area, sum)
      
      d.out <- cbind(d.out,c(pred.area,neg.area))
    }
  p.out <- d.out/apply(d.out,1,sum)
  d.out <- cbind(d.out,p.out)
  
  colnames(d.out) <- c(paste("n.var_",pos,sep=""),paste("p.var_",pos,sep=""))
  row.names(d.out) <- c("pred.area","f.neg.area")
  
  write.table(d.out, file=paste(outpath.sens,list.name[n],"_sen.res.txt",sep=""),sep="\t",
              append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
    
  }# end loop for zone

#############################################################################################
### reload data and make tablefor f.neg only
#############################################

# list of file
fns <- list.files(path=outpath.sens,
                  pattern = ".tif$",full.names = TRUE)

list.name <- list.files(path=outpath.sens,
                        pattern = ".tif$",full.names = FALSE)

for (k in c(11:17,19) )# file containing f.neg.only 
{
  r <- raster(fns[k]) # load raster
  
  ens.area <- abs(r) 
  
  pos <- unique(ens.area)
  d.out <-NULL
  for (q in 1:length(pos))
      {
        pred.area <- reclassify(ens.area,c(-Inf,pos[q]-1, 0, 
                                           pos[q]-1,pos[q],1,
                                           pos[q],Inf,0 ))
        pred.area <- cellStats(pred.area, sum)
        
        d.out <- cbind(d.out, pred.area)
      }
  
  p.out <- d.out/apply(d.out,1,sum)
  d.out <- rbind(d.out,p.out)
  
  colnames(d.out) <- pos
  row.names(d.out) <- c("n.var_","p.var_")
  
  # order
  pos <- order(d.out[2,],decreasing =TRUE)
  d.out <- d.out[,pos]
  
  write.table(d.out, file=paste(outpath.sens,list.name[k],".txt",sep=""),sep="\t",
              append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
}
  
  