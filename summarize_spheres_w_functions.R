# Summarize TLS on multiple spheres. 
# 2021 April
# Aaron Maxwell & Luis Andres Guillen

# AIM: This codes analysis TLS using spherical transformation to understand
# point density at different radius from the scan
# ARGUMENTS: Variables that need to be uploaded are the LAS files, radius and constants. 
# OUTPUT: Raster object or tif file.

#Loading libraries ----
library(rlas)
library(dplyr)
library(ggplot2)
library(plotly)
library(geometry)
library(MASS)
library(raster)
library(sampSurf)
library(rgdal)
library(raster)
library(lidR)
library(tmap)
library(sf)

#Declaring variables ----
#Read in normalized LAS data

tlsData <- readLAS("../FS_fuel_prediction/raw_data/Compartment10_laz_files/C10_Plot_10East.laz")
tlsData<-clip_circle(tlsData,0,0,10)
#Define radius of sphere
radii <- c(1, 2, 3, 4, 5)#, 6, 7, 8, 9, 10) 

#Declaring constants ----
height = 2
thetaDiv=.5
phiDiv=.5
grdDist=.25
r_blank <- raster(ncol= 360/thetaDiv, nrow = 180/phiDiv, 
                  xmn=-180, xmx=180,ymn=-90,ymx=90)

# Declaring functions ----

cld_fun<-function(r,tlsData){
  #Aim:Calculate coordinates where each ray would intersect the sphere
  #Arguments:r = radius; tlsData: LAS object
  #Output: cld dataframe. 
  cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z - height))
  cld$a <- cld$x^2 + cld$y^2 + cld$z^2
  cld$t1 <- -sqrt(0^2-(4*cld$a*(-r^2)))/(2*cld$a)
  cld$t2 <- sqrt(0^2-(4*cld$a*(-r^2)))/(2*cld$a)
  cld$x1 <- cld$x*cld$t1
  cld$y1 <- cld$y*cld$t1
  cld$z1 <- cld$z*cld$t1
  cld$x2 <- cld$x*cld$t2
  cld$y2 <- cld$y*cld$t2
  cld$z2 <- cld$z*cld$t2
  cld$slp <- cld$z/(sqrt(cld$x^2+ cld$y^2))
  cld$slp1 <- cld$z1/(sqrt(cld$x1^2+ cld$y1^2))
  cld$slp2 <- cld$z2/(sqrt(cld$x2^2+ cld$y2^2))
  
  #Find correct coordinate based on slope
  cld<-cld %>% mutate(xs = case_when(sign(slp)==sign(slp1)~ x1,sign(slp)!=sign(slp1)~x2),
                      ys = case_when(sign(slp)==sign(slp1)~ y1,sign(slp)!=sign(slp1)~y2),
                      zs = case_when(sign(slp)==sign(slp1)~ z1,sign(slp)!=sign(slp1)~z2))
  return(cld)    
}

sphere_func<-function(cld=cld){
  #Aim: calculate spherical coordinates. 
  #Arguments: cld = dataframe with cartesian distances. 
  #Output: list of objects s3 and p3. 
  
  #Filter points that would hit the ground before hitting the sphere
  just_sphere <- cld %>% filter(cld$zs > -height)
  
  #Filter out points that would strike the sphere before the ground
  just_grd <- cld %>% filter(cld$zs <= -height)
  
  #Calculate Euclidean distance from sensor to point
  just_sphere$dist_act <- sqrt(just_sphere$x^2 + 
                                 just_sphere$y^2 + 
                                 just_sphere$z^2)
  #Calculate Euclidean distance from sensor to sphere/point intersection
  just_sphere$dist_sphere <- sqrt(just_sphere$xs^2 + 
                                    just_sphere$ys^2 +
                                    just_sphere$zs^2)
  
  #Filter out points that actually pass through the sphere
  pass_sphere <- just_sphere %>% filter(dist_act >= dist_sphere)
  
  pass_sphere_g <- just_grd %>% filter(z <= -height)
  
  #Create table of x, y, z values where each point did or would have intersected the ground
  allF <- bind_rows(just_sphere, just_grd)
  
  #Create table of x, y, z values where each point that passed the ground intersected it
  passF <- bind_rows(pass_sphere, pass_sphere_g)
  
  
  #Convert from Cartesian to spherical coordinates
  s <- as.data.frame(cart2sph(x=allF$xs, 
                              y=allF$ys, 
                              z=allF$zs))
  p <- as.data.frame(cart2sph(x=passF$xs, 
                              y=passF$ys, 
                              z=passF$zs))
  
  s2 <- s[,c(1,2)]
  p2 <- p[,c(1,2)]
  
  s2$theta <- s$theta*(180/pi)
  s2$phi <- s$phi*(180/pi)
  p2$theta <- p$theta*(180/pi)
  p2$phi <- p$phi*(180/pi)
  
  s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
  p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
  
  multipoint_list<-list(s3,p3) 
  names(multipoint_list)<-c("s3","p3")
  return(multipoint_list)
} 

pointCount <- function(r, pts){
  #From: https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r
  
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

on_sphere_multi_v2 <- function(tlsData, radii, height, thetaDiv=1, phiDiv=1){
  for(ra in radii){
    if(ra == radii[1]){
      r <- ra
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func()
      
      gridS <- pointCount(r_blank,list_point$s3)
      gridP <- pointCount(r_blank, list_point$p3)

      null1 <- raster::calc(gridS, function(x){x[x<=0]<-NA; return(x)})
      null2 <- null1 > 0

      out <- (1-(gridP/gridS))*null2
    }else if(ra == radii[2]){
      r<-ra
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func()
      
      gridPre <- gridS - gridP
    
      gridS <- pointCount(r_blank, list_point$s3) - gridPre
      gridP <- pointCount(r_blank, list_point$p3) 
    
      mask1 <- out
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)})
      mask3 <- mask2 >= 0
    
      out <- (1-(gridP/gridS))
      out2 <- stack(mask1, out)
    }else{
      r<-ra
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func()
      gridPre <- gridS - gridP
      
      gridS <- pointCount(r_blank, list_point$s3) - gridPre
      gridP <- pointCount(r_blank, list_point$p3) 
      
      mask1 <- out
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)})
      mask3 <- mask2 >= 0
      
      out <- (1-(gridP/gridS))*mask3
      out2 <- stack(out2, out)
    }
  }
  names(out2) <- as.character(radii)
  return(out2)
}

# Applying functions ----

test <- on_sphere_multi_v2(tlsData=tlsData, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)

writeRaster(test, "D:/sphere_vox.tif")