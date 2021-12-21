###################################################
## Compute incoming solar insolation with RGrass 7.8.5
## work under R 4.0.3
###################################################
# v2.1 - Boris Droz - February 2012
# v2.2 - integre wgs turbidity raster
###################################################################################################
## BEFORE FIRST RUN ####################
## -------------- #########
## --> RGRASS7 requires GRASS v.7 
##  ----- DOWNLOAD IT AT : -------
##     https://grass.osgeo.org
#########################################
## turbidity raster file should be downloaded from
## http://www.soda-pro.com/help/general-knowledge/linke-turbidity-factor
## and put in the folder ... script/linke_value
#######################################
#######################################
## Descriptif
##############
# Comput averaged daily incoming solar insolation (beam radiation) 
# on the basic of a DEM - surface model (mns minus mnt)
# for the month April - Mai -June (Rafb breeding period)
#
# Set the parameter and run the script

## Reference 
#############
## GRASS start up: https://grasswiki.osgeo.org/wiki/R_statistics/rgrass7
### info function https://grass.osgeo.org/grass79/manuals/r.sun.html
#####################################################################################################################################

#####################################################################################################################################
#####################################################################################################################################
##### PARAMETER #################
####################################################################

working_dir <- "~/Actual/rafb_proj_model/analysis/Cdf_2021/"

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
n.proj <- "LV95" # name projection

# tubidity from Remund J., Wald L., Lefevre M., Ranchin T., Page J., 2003.
sele.linke_value <- "NO" # YES automatic selection NO gave value in l_k
l_k <-c(4.1,4.4,5.1)  # depending on time and location see ref

#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################
# find out working path
########################
setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","magrittr","rgdal","rgeos",
             "rgrass7","RNetCDF","splitstackshape","maptools","XML"))
use_sp()
#####################################################################################################################################
## START AND END
################
## period between April and June was selected because 
## it corresponds to the nesting period of the Common Redstart 
day <- c("1Apr16","30Apr16","1May16","31May16","1Jun16","30Jun16")
tmp <- as.POSIXlt(day, format = "%d%b%y")

# Folder and file
#################
inpath <- paste(working_dir,"/input",sep="")
creat.subDir(working_dir,"/pred/pred_cur")
creat.subDir(inpath,"/monthly_rad")
outpath <- paste(working_dir,"/pred/pred_cur",sep="")
setwd(inpath)
#####################################################################################################################################
############################################################################
###########r####################################################################
beginCluster( ) # ACTIVATE THIS MULTI CORE CALCULATION FOR RASTER

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

## CALCULATE RAD 
###################
# Set GRASS environment and database location 
##############################################
loc <- initGRASS(gisBase ="C:/Program Files/GRASS GIS 7.8", 
                 home=getwd(),  override=TRUE, mapset = "PERMANENT" )

## and baseline mapset with required epsg code
execGRASS("g.proj", flags = "c", epsg = proj)
###############################################

#ptm <- proc.time()# ignite timer

# get file dem
fns.dem <- paste(getwd(),"/dem_",n.proj,"_2m.tif",sep="")

# get linke file value
fns.link <- list.files(paste(working_dir,"/script/linke_value/",sep="")
                       ,pattern="tif$",full.names = TRUE)

# Import raster to GRASS, set region and projection
execGRASS("r.in.gdal", flags=c("o","overwrite"), 
            parameters=list(input= fns.dem, output="tmprast"))
execGRASS("g.region", parameters=list(raster="tmprast") ) 

# creat lat long raster
execGRASS("r.latlong",flags=c("overwrite"),
             parameters=list(input= "tmprast", output="r.lat"))
execGRASS("r.latlong",flags=c("overwrite","l"),
             parameters=list(input= "tmprast", output="r.long"))

# get linke_value
# # turbidity from Remund J., Wald L., Lefevre M., Ranchin T., Page J., 2003.
## http://www.soda-pro.com/help/general-knowledge/linke-turbidity-factor
if (sele.linke_value == "NO") {  }else{
  
  # get mean lat and long
  lat <- as.numeric(sub('.*:', '',summary(readRAST("r.lat"))[]$data[4,]))
  long <- as.numeric(sub('.*:', '',summary(readRAST("r.long"))[]$data[4,]))
  
  #get values from raster
  l_k <- NULL
  for (q in 1:3)
    {
     r.link <- raster(fns.link[q])
     
     l_k <- c(l_k, getValues(r.link,row=rowFromY(r.link, lat) ) [colFromX(r.link, long)]/20 )
    }
  }

## compute input parameter for solar radiation
# slope + aspect
execGRASS("r.slope.aspect",  flags=c("overwrite"),
          parameters=list(elevation= "tmprast",  aspect="aspect.dem", slope="slope.dem") )
 
	p <-0 # count   

    for (h in c(1,3,5))
      {
	p<-p+1

      t.d <- seq(from=tmp$yday[h],to=tmp$yday[h+1],by=1)
      
      for (z in 1:length(t.d))
        {
        cat("> Compute Julian day ",t.d[z],"...", "\n",append = FALSE)
        
        if (z==1) {
          ## compute global solar direct radiance per monthly averaged for the month of interest
          execGRASS("r.sun", flags=c("quiet","overwrite"), 
                    parameters=list(elevation= "tmprast", aspect ="aspect.dem", 
                                    slope="slope.dem", lat="r.lat",long="r.long",
                                    day=t.d[z],
                                    # defaut values from  Area Solar Radiation from arcgis
                                    # diffuse prop. (albedo) and transmittance value 
                                    # Solar constant of 1367 W.m-2 is used
                                    albedo_value=0.3, linke_value= l_k[p], 
						                        beam_rad = "radout"))
                                    #glob_rad="radout")) 
          } else{
            execGRASS("r.sun", flags=c("quiet","overwrite"), 
                      parameters=list(elevation= "tmprast", aspect ="aspect.dem", 
                                      slope="slope.dem",lat="r.lat",long="r.long",
                                      day=t.d[z],
						                          albedo_value=0.3, linke_value= l_k[p],
						                          beam_rad = "radout2"))
                                      #glob_rad="radout2"))
            
            execGRASS("r.mapcalc", flags=c("overwrite"), 
                        parameters= list(expression="radout=radout+radout2") )
          }
        } # end sum of a month
          
        # save the monthly solar
        r <- readRAST("radout") # read the grass format
        write.asciigrid(r, paste(inpath,"/monthly_rad/sumrad_",format(tmp,"%m")[h],".tif",sep=""),
                        attr = 1) # write the raster in tif format
                        
        r.t <- raster(paste(inpath,"/monthly_rad/sumrad_",format(tmp,"%m")[h],".tif",sep=""))
        
        if (h==1) { r.out <- r.t /(tmp$yday[h+1]- tmp$yday[h]+1) # per day  
                }else{
                    r.out <- r.out + r.t /(tmp$yday[h+1]- tmp$yday[h]+1)   
                }
      } # end loop months

names(r.out) <- "fsumrad_2m"
crs(r.out)<- CRS(paste("+init=epsg:", proj,sep="") )
          
writeRaster(r.out, filename = paste(outpath,"/fsumrad_2m.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)   

## CHECK ---
r.out;plot(r.out)

file.remove(dir(tmpDir(), full.names=TRUE)) # REMOVE TEMPORARY RASTER 

endCluster() # END OF MULTICORE CALCULATION

#proc.time() - ptm # check time
#####################################################################################################################################
################  END SCRIPT ###############################################
###########################################################################



