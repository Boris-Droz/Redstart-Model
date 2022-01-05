######################################################################
##
##  FUNCTION TO RESTRICT AREA OF PREDICTION 
##   BASED ON THE BOUNDARIES VALUES USED 
##   DURING THE CALIBRATION OF THE MODEL
##
######################################################################
######################################################################
##
## v1.1 - B.Droz March 2021
##

restri.to_cali_area <-  function(f.path,env.stack.proj)
{
  require(raster)
  
  cali.pred <- read.table(f.path, header =TRUE)
  
  names.cali.pred <- colnames(cali.pred)
  names.raster <- names(env.stack.proj)
  
  p<-0
  
  for (i in 1:ncol(cali.pred)) 
  {
    pos <- names.cali.pred[i]== names.raster
    
    if (any(pos)){
      
      p<-p+1
      
      f.min <- cali.pred[1,i]
      f.max <- cali.pred[3,i]
      
      ## reclassifed to be in calibration boundary
      fil.min <- reclassify( subset( env.stack.proj,names.raster[pos] ), 
                             c(-Inf,f.min,NA, f.min,Inf,1 ),include.lowest=TRUE )
      
      fil.max <- reclassify( subset( env.stack.proj,names.raster[pos] ), 
                             c(-Inf,f.max,1, f.max,Inf,NA),include.lowest=FALSE )                        
      
      if (p==1) { f.out <- stack(fil.min,fil.max) }else{
        f.out <- stack(f.out,fil.min,fil.max) }
    }else{}
    
  }
  # multiplied all raster stacked
  for (j in 1:nlayers(f.out) ) {
      if (j==1) { r.out <- f.out[[j]] 
      }else{ r.out <- r.out * f.out[[j]] }
    }
  
  env.stack.proj <- r.out * env.stack.proj
  
  names.raster -> names(env.stack.proj)
  
  return(env.stack.proj)
}



