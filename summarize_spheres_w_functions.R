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
height = 2 #height of the scanner 
thetaDiv=.5 #
phiDiv=.5 #
grdDist=.25 #
r_blank <- raster(ncol= 360/thetaDiv, nrow = 180/phiDiv, 
                  xmn=-180, xmx=180,ymn=-90,ymx=90) # blank raster to fill in data

# Declaring functions ----

cld_fun<-function(r,tlsData){
  #Aim: Calculate coordinates where each ray would intersect the sphere
  #Arguments:r = radius; tlsData: LAS object
  #Output: cld dataframe. 
  
  # create a dataframe to calculate distances to the scanner
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
  #Arguments: cld: dataframe with cartesian distances obtained with cld_fun
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
  
  s3 <- st_multipoint(x=as.matrix(s2), dim="XY") # make spatial object
  p3 <- st_multipoint(x=as.matrix(p2), dim="XY") # make spatial object
  
  #create lists and finalize function
  multipoint_list<-list(s3,p3)
  names(multipoint_list)<-c("s3","p3")
  return(multipoint_list)
} 

pointCount <- function(r, pts){
  # Idea source: https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r
  # Aim: Calculate the number of points in a raster cell grid. 
  # Argument: r: radious
  # Output: r2: raster object with points. 
  
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
  # Aim: Transforms a point cloud into a raster showing the points that passed or not pass a certain radious, 
  # depending on the layer of the raster.
  # Arguments: tlsData: cloudpoint, radii: list of radious, height: height of scanner.
  # Output: raster object with different layers, representing each radious. 
  for(ra in radii){ # Begin for radious loop 
    if(ra == radii[1]){ # Logic statement of 1st radius
      r <- ra #read radius
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func() # apply functions to obtain spherical coordinates 
      
      gridS <- pointCount(r_blank,list_point$s3) #create point object of points on the sphere
      gridP <- pointCount(r_blank, list_point$p3)#create point object of points that pass the sphere

      null1 <- raster::calc(gridS, function(x){x[x<=0]<-NA; return(x)}) # make raster??
      null2 <- null1 > 0 #??

      out <- (1-(gridP/gridS))*null2 #create raster for first radius
    }else if(ra == radii[2]){# logical statement 2nd radius
      r<-ra # read radius
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func() # apply functions to obtain spherical coordinates 
      
      gridPre <- gridS - gridP # make temporal raster of difference. 
    
      gridS <- pointCount(r_blank, list_point$s3) - gridPre # create point object of points on the sphere minus the previous sphere
      gridP <- pointCount(r_blank, list_point$p3) #create point object of points that pass the sphere
    
      mask1 <- out # read previous layer. 
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)}) #??
      mask3 <- mask2 >= 0 #??
    
      out <- (1-(gridP/gridS)) # make faster for 2nd radius
      out2 <- stack(mask1, out) # stack rasters
    }else{ # logic stament n_th radius
      r<-ra # read data
      list_point<-cld_fun(r=r,tlsData = tlsData) %>% sphere_func() # apply functions to obtain spherical coordinates
      
      gridPre <- gridS - gridP# make temporal raster of difference.
      
      gridS <- pointCount(r_blank, list_point$s3) - gridPre # create point object of points on the sphere minus the previous sphere
      gridP <- pointCount(r_blank, list_point$p3) #create point object of points that pass the sphere
      
      mask1 <- out # read raster of previous radius.
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)}) #??
      mask3 <- mask2 >= 0 #??
      
      out <- (1-(gridP/gridS))*mask3 # create raster of n radius. 
      out2 <- stack(out2, out) # stack rasters of n radius 
    }
  }#Close radius loop
  names(out2) <- as.character(radii)# name layers
  return(out2) # set function output. 
} 

# Applying functions ----

test <- on_sphere_multi_v2(tlsData=tlsData, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)

writeRaster(test, "D:/sphere_vox.tif") # save data. 