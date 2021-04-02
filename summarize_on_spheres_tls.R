#Load Libraries
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

#From: https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r
pointCount <- function(r, pts){
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

#Read in normalized LAS data
tlsData <- readLAS("D:/lidar/normalize/norm_clf_BLK360 Scan.las")
#Define radius of sphere
radius <- 10 
#Create points in the 3D space
x <- tlsData@data$X
y <- tlsData@data$Y
z <- tlsData@data$Z #Subtract 2 for sensor height

height = 2
thetaDiv=.5
phiDiv=.5

on_sphere <- function(tlsData, radius, height, thetaDiv=1, phiDiv=1){
  r <- radius
  x <- tlsData@data$X
  y <- tlsData@data$Y
  z <- tlsData@data$Z - height
  cld <- as.data.frame(cbind(x=x, y=y, z=z))
  
  #Calculate coordinates where each ray would intersect the sphere
  cld$a <- cld$x^2 + cld$y^2 + cld$z^2
  cld$b <- 0
  cld$c <- -r^2
  cld$t1 <- -cld$b - sqrt(cld$b^2-(4*cld$a*cld$c))/(2*cld$a)
  cld$t2 <- -cld$b + sqrt(cld$b^2-(4*cld$a*cld$c))/(2*cld$a)
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
  for(i in 1:nrow(cld)){
    if(sign(cld$slp[i]) == sign(cld$slp1[i])){
      cld$xs <- cld$x1
      cld$ys <- cld$y2
      cld$zs <- cld$z2
    }else{
      cld$xs <- cld$x2
      cld$ys <- cld$y2
      cld$zs <- cld$z2
    }
  }

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
  just_sphere$xf <- just_sphere$xs 
  just_sphere$yf <- just_sphere$ys 
  just_sphere$zf <- just_sphere$zs 
  just_grd$xf <- just_grd$xs
  just_grd$yf <- just_grd$ys
  just_grd$zf <- just_grd$zs
  allF <- bind_rows(just_sphere, just_grd)

  #Create table of x, y, z values where each point that passed the ground intersected it
  pass_sphere$xf <- pass_sphere$xs 
  pass_sphere$yf <- pass_sphere$ys 
  pass_sphere$zf <- pass_sphere$zs 
  pass_sphere_g$xf <- pass_sphere_g$xs
  pass_sphere_g$yf <- pass_sphere_g$ys
  pass_sphere_g$zf <- pass_sphere_g$zs
  passF <- bind_rows(pass_sphere, pass_sphere_g)


  #Convert from Cartesian to spherical coordinates
  s <- as.data.frame(cart2sph(x=allF$xf, 
                              y=allF$yf, 
                              z=allF$zf))
  p <- as.data.frame(cart2sph(x=passF$xf, 
                              y=passF$yf, 
                              z=passF$zf))
  
  s2 <- s[,c(1,2)]
  p2 <- p[,c(1,2)]
  
  s2$theta <- s$theta*(180/pi)
  s2$phi <- s$phi*(180/pi)
  p2$theta <- p$theta*(180/pi)
  p2$phi <- p$phi*(180/pi)
  
  r_blank <- raster(ncol= 360/thetaDiv, 
                    nrow = 180/phiDiv, 
                    xmn=-180, xmx=180, 
                    ymn=-90, 
                    ymx=90)

  s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
  p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
  
  gridS <- pointCount(r_blank, s3)
  gridP <- pointCount(r_blank, p3)
  
  null1 <- raster::calc(gridS, function(x){x[x<=0]<-NA; return(x)})
  null2 <- null1 > 0
  
  s_prop <- (1-(gridP/gridS))*null2
}

testOut <- on_sphere(tlsData=tlsData, radius=radius, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)
