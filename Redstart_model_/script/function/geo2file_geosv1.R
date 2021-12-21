#!/usr/bin/env python 
# author: Bruno COMBAL, JRC, European Commission 
# date: 4/12/2008 
# purpose: reproject a GEOS image (MSG) into geographic
# source: https://publications.jrc.ec.europa.eu/repository/bitstream/JRC52438/combal_noel_msg_final.pdf

## adapte for R code by B.Droz - February 2021
##############################################################

geo2file_geos <- function (lon,lat,rflon=0)
  {
  
  #    """ #    input 
  #    vlon: longitude 
  #    vlat: latitude 
  #    rflon: satellite longitude, default=0degrees 
  #    output xr: x position 
  #    yr: y position 
  #    """ 
  
  #geo2file_geos(lon, lat)     
  # Setup constants     
  re = 6378.160           # equatorial radius     
  h  = 42164.0 - re       # Reference altitude     
  rp = 6356.5838     
  lpsi2  = 1              # spin directi
  
  resol = 3712.0  # pixel/line number for M
  
  deltax = 17.83 / resol     # E-W scanning step   scan amplitude is 17.83 deg. for MSG
  deltay = 17.83 / resol     
  
  #dtor = 1*pi/180 # convert degree in radian
  dtor =1 # r function are in degree....
  xlat = dtor * lat     
  xlon = dtor * lon     
  cosxlat = cos(xlat)     
  sinxlat = sin(xlat)     
  cosxlon = cos(xlon)     
  rom = (re*rp) / sqrt( rp*rp*cosxlat*cosxlat + re*re*sinxlat*sinxlat )         
  y   = sqrt( h * h + rom * rom - 2.0 * h * rom * cosxlat * cosxlon )     
  r1  = y*y + rom*rom     
  r2  = h * h     
  
  if (r1 > r2){return (list(-1, -1) )}
  
  rs   = re + h     
  reph = re     
  rpph = rp     
  coslo = cos(rflon * dtor)     
  sinlo = sin(rflon * dtor)     
  teta  = atan( (rpph/reph) * tan( xlat ) )     
  xt    = reph * cos( teta ) * cosxlon      
  yt    = reph * cos( teta ) * sin( xlon )     
  zt    = rpph * sin( teta )         
  px    =  atan(  (coslo  *(yt  -  rs*sinlo)  -  (xt  -  rs*coslo)*sinlo)
                       /(sinlo*(yt  -  rs*sinlo) + (xt - rs*coslo)*coslo) )     
  py    = atan( zt * ( (tan(px)*sinlo - coslo) / (xt-rs*coslo)) *cos(px) )     
  px    = px / dtor     
  py = py / dtor     
  xr    = px / (deltax*lpsi2)     
  yr    = py / (deltay*lpsi2)
  
  if (xr >= 0)
    { xr =  px / (deltax*lpsi2)  + 0.5
    } else{  xr =  px / (deltax*lpsi2)  - 0.5 }
  
  if (yr >= 0){ yr =  py / (deltax*lpsi2)  + 0.5
    } else{  yr =  py / (deltax*lpsi2)  - 0.5 }
  
  xr = xr + 0.5*resol + 0.5     
  yr = yr + 0.5*resol + 0.5
  
  return (c(xr, yr))
}
##test 
# geo2file_geos (lon=6.8,lat=47.1,rflon=0)
