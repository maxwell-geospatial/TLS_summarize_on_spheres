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
#tlsData <- decimate_points(tlsData1, homogenize(10, 1))

#Define radius of sphere
radii <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) 
#Create points in the 3D space
x <- tlsData@data$X
y <- tlsData@data$Y
z <- tlsData@data$Z #Subtract 2 for sensor height

height = 2
thetaDiv=.5
phiDiv=.5
grdDist=.25


on_sphere_multi <- function(tlsData, radii, height, thetaDiv=1, phiDiv=1){
  for(ra in radii){
    if(ra == radii[1]){
      r <- ra
      cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z - height))
      
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
      cld<-cld %>% mutate(xs = case_when(sign(slp)==sign(slp1)~ x1,sign(slp)!=sign(slp1)~x2),
                          ys = case_when(sign(slp)==sign(slp1)~ y1,sign(slp)!=sign(slp1)~y2),
                          zs = case_when(sign(slp)==sign(slp1)~ z1,sign(slp)!=sign(slp1)~z2))

      
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
      
      out <- (1-(gridP/gridS))*null2 #raster
    }else if(ra == radii[2]){
      r <- ra
      cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z - height))
      
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
      cld<-cld %>% mutate(xs = case_when(sign(slp)==sign(slp1)~ x1,sign(slp)!=sign(slp1)~x2),
                          ys = case_when(sign(slp)==sign(slp1)~ y1,sign(slp)!=sign(slp1)~y2),
                          zs = case_when(sign(slp)==sign(slp1)~ z1,sign(slp)!=sign(slp1)~z2))
      
      
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
      
      r_blank <- raster(ncol= 360/thetaDiv, 
                        nrow = 180/phiDiv, 
                        xmn=-180, xmx=180, 
                        ymn=-90, 
                        ymx=90)
      
      s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
      p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
      
      gridPre <- gridS - gridP #substracting all - passed 1radius
      
      gridS <- pointCount(r_blank, s3) - gridPre
      gridP <- pointCount(r_blank, p3) 
      
      mask1 <- out
      mask2 <- raster::calc(mask1, function(x){x[x==1]<-NA; return(x)})
      mask3 <- mask2 >= 0
      
      out <- (1-(gridP/gridS))
      out2 <- stack(mask1, out)
    }else{
      r <- ra
      cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z - height))
      
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
      cld<-cld %>% mutate(xs = case_when(sign(slp)==sign(slp1)~ x1,sign(slp)!=sign(slp1)~x2),
                          ys = case_when(sign(slp)==sign(slp1)~ y1,sign(slp)!=sign(slp1)~y2),
                          zs = case_when(sign(slp)==sign(slp1)~ z1,sign(slp)!=sign(slp1)~z2))
      
      
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
      
      r_blank <- raster(ncol= 360/thetaDiv, 
                        nrow = 180/phiDiv, 
                        xmn=-180, xmx=180, 
                        ymn=-90, 
                        ymx=90)
      
      s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
      p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
      
      gridPre <- gridS - gridP
      
      gridS <- pointCount(r_blank, s3) - gridPre
      gridP <- pointCount(r_blank, p3) 
      
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



test <- on_sphere_multi(tlsData=tlsData, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)

writeRaster(test, "D:/sphere_vox.tif")


