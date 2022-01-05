#####
#### Plan Affectation
####=================================================
####
## work under R 4.0.3
## v1.1 B.Droz April 2021

###################################################################################################
## Descriptif
##############
#  define amenagement where to applied possible scenario
## we restrict the scenario by excluding farmland area, downhill resort, forest.
## 
#####################################################################################################################################
##### PARAMETER #################
####################################################################

proj <- 2056 # input projection new swiss projection 21781 # old swiiss
n.proj <- "LV95" # name projection
#f.type <- ".asc$" # type of file of the mnt
working_dir <-"~/Actual/rafb_proj_model/analysis/Planchette_2021/"

#####################################################################################################################################
#####################################################################################################################################
## ---- SCRIPT START HERE -- do it manually step by step
#####################################################################################################################################
#####################################################################################
#####################################################################################
## Load library
##################}
setwd(working_dir)
source("script/function/creat_subDirv1.R")
source("script/function/packagesv1.0.R")

check.lib (c("raster","sp","magrittr","rgdal","rgeos","RNetCDF",
             "splitstackshape","maptools","XML","grainchanger","gstat"))

# set file location 
inpath <- paste(working_dir,"/input/",sep="")
outpath <- paste(working_dir,"/zone/",sep="")
setwd(inpath)
###################################################################################
#####################################################################################################################################
## OPEN Data
################
PA <- rgdal::readOGR(file.choose()) ## should be a shp file "principal use  "
b.proj.area <- raster(paste(working_dir,"/zone/buff.proj_LV95_2m.tif",sep="") )
proj.area <- raster(paste(working_dir,"/zone/zone.mod_LV95_2m.tif",sep="") )

# set the coordinate system
proj4string(PA) <- CRS(paste("+init=epsg:", proj,sep="") )

# any bad polygon --> resolve "invalid: Ring Self-intersection"
if (sum(gIsValid(PA, byid=TRUE)==FALSE) >0) 
{ PA <- gBuffer(PA, byid=TRUE, width=0)}else{ }

# crop to zone of pred
c.PA <- crop(PA,b.proj.area)

####################
## VD classification
####################
unique(c.PA$TYPE_PRINC)

c.PA <- c.PA[c.PA$TYPE_PRINC!="Zone agricole protégée" &
               c.PA$TYPE_PRINC!="Zone agricole" &
               c.PA$TYPE_PRINC!="DP" &
               c.PA$TYPE_PRINC!="Aire forestière" &
               c.PA$TYPE_PRINC!="Zone de piste de ski" &
             c.PA$TYPE_PRINC!="Zone d'extraction et dépôt de matériaux" &
             c.PA$TYPE_PRINC!="Zone équestre",
             ]

####################
## NE classification
####################
unique(c.PA$desprincip)

c.PA <- c.PA[c.PA$desprincip!="Zone de transport" &
               c.PA$desprincip!="Zone de maintien de l'habitat rural" &
               c.PA$desprincip!="Zone agricole / Aire forestière / Cours d'eau et étendue d'eau / Espace de transport" &
               c.PA$desprincip!="Zone d'extraction de matériaux" &
               c.PA$desprincip!="Zone de traitement des déchets" &
               c.PA$desprincip!="Zone d'aérodrome",
                ]

####################
## FR classification -- geodienst
####################
unique(c.PA$desprincip)


c.PA <- c.PA[c.PA$desprincip!="Zone de transport" &
               c.PA$desprincip!="Zone de maintien de l'habitat rural" &
               c.PA$desprincip!="Zone agricole / Aire forestière / Cours d'eau et étendue d'eau / Espace de transport" &
               c.PA$desprincip!="Zone d'extraction de matériaux" &
               c.PA$desprincip!="Zone de traitement des déchets" &
               c.PA$desprincip!="Zone d'aérodrome",
]
###################################################################################
##
# rasterize
r.pa <- rasterize(c.PA,b.proj.area,field=1,background=0)

w.dur <- crop(r.pa,proj.area)

## if raster still do not have the same exact extent
if ( extent(w.dur) == extent(proj.area) ) 
  { }else{
    empty.raster <- raster(res=res(proj.area),ext=extent(proj.area), crs= crs(proj.area )) 
    w.dur <- resample(w.dur,empty.raster,method="ngb")
  }

w.dur <- raster::mask(w.dur,proj.area)

names(w.dur) <- "PA"

writeRaster(w.dur, filename = paste(outpath,"zone_pa.tif",sep=""),
            datatype="FLT8S", overwrite=TRUE)

# --- check --
print(w.dur);plot(w.dur)

