
### Creat a clean shape from a raster
## tutoreil in 
##  https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html
#  https://stackoverflow.com/questions/52911022/dissolve-boundaries-of-polygon-features-in-r
## 
############################################################################
clean.rasttopoly <- function(r,  # raster 
                             surf.min=surf.min,fil.hole=fil.hole, # both in m2
                             method = "ksmooth" , nb.smooth=200,# smoth parameter
                             d.con=0.001, # distance to connect two areas 
                                          # use 0.001 to connect only if adajacent!
                             concave_distance_lim = 2,# hull polygone
                             d.terri=50)  # distance to conc poly 
                             
{ # check lib
  require(raster)
  require(sf)
  require(units)
  require(smoothr)
  require(rgeos)
  require(rangemap)# hull polygone
  
  # shape the raster
  shp.proj <- rasterToPolygons(r, dissolve=FALSE)
  
  #apply a negligible buffer distance around polygon features to combine features 
  #that are only connected by one raster corner (and therefore counted 
  #incorrectly as distinct features after the above steps)
  #and dissolve
  shp.proj <- buffer(shp.proj, width = d.con, dissolve = TRUE)
  
  #disaggregate multi-part buffered polygon into separate features
  shp.proj <- disaggregate(shp.proj) ## too much time
  
  # fill hole
  shp.proj <- fill_holes(shp.proj, set_units(fil.hole, m^2))
  
  #####################################################################
  ## creat concave polygone with hull polygone 
  ############################################
  spl <- NULL
  for (i in 1:length(shp.proj))
    {
      pol <- st_as_sf(shp.proj[i,]) # convert into st object
    
      xy <- st_coordinates(pol)
    
      xy <- as.data.frame(unique(xy))
    
      xy.sp <- sp::SpatialPointsDataFrame(coords = xy[, 1:2], data = xy,
                                    proj4string = crs(r)  )
    
      spl[[i]] <- hull_polygon(xy.sp, hull_type = "concave", 
                             concave_distance_lim = concave_distance_lim,
                 verbose = TRUE) # rangemap
    }
  
  spl <- do.call(bind, spl) # merge all spatial polygon
  
  # buffer all polygone to merge it by define threshold
  #####################################################
  b.spl <- buffer(spl, width=d.terri, dissolve = FALSE)
  
  sp2 <- NULL
  q<-0
  for (i in 1:length(b.spl))
    {
      pol1 <- b.spl[i,]
      
      n <- 1:length(b.spl)
      n <- n[n!=i]
      
      for (j in n) 
        {
          q<- q+1
          
          pol2 <- b.spl[j,]
          sp2[[q]] <- gIntersection(pol1,pol2)
        }
    }
  
  sp2 <- sp2[-which(sapply(sp2, is.null))] # remove null feature
  
  sp2 <- do.call(bind, sp2) # merge all spatial polygon
  
  sp2 <- bind(spl,sp2)  # merge all intersection with the original
   ############### temp save done...
  
  #aggregate to remove stack layer 
  sp2 <- aggregate(sp2,dissolve=TRUE)
  
  #disaggregate multi-part buffered polygon into separate features
  sp2 <- disaggregate(sp2) ## too much time
  
  # find and delet all area size too small
  sp2$area_sqm <- raster::area(sp2)   
  sp2 <- sp2[sp2$area_sqm>surf.min,]
  
  # fill hole again
  sp2 <- fill_holes(sp2, set_units(fil.hole, m^2))
  
  # smooth the line of the shape -take time
  sp2 <- smooth(sp2, method = method, smoothness = nb.smooth)
  
  return(sp2)
}         


