# Summarize TLS on multiple spheres. 
# 2021 April
# Aaron Maxwell & Luis Andres Guillen

# AIM: This codes analysis TLS using spherical transformation to understand
# point density at different radius from the scan
# ARGUMENTS: Variables that need to be uploaded are the LAS files, radius, and constants. 
# OUTPUT: multi-band raster object or tif file.
#Each cell contains an estimate of the proportion of the 
#returns that reached it that were returned by a feature in the cell. 
#Each cell is defined by an angle of theta and phi. 
#Each layer in the raster stack represents a distance from the sensor. 

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
#Read in shed scan reference
#refScan <- readLAS("D:/lidar/Shedtest1.las")

#mycsf <- csf(TRUE, 0.5, 0.5, time_step = 0.65) #Define settings for ground segmentation algorithm
#las_clf <- classify_ground(refScan, mycsf) #Classify to ground (2) and not ground (1)

#las_norm <- normalize_height(las_clf, tin()) #Normalize Z relative to ground (Z is now height above ground)

#writeLAS(las_norm, "D:/lidar/normalize/shed_norm.las") #Save normalized file to disk

#Read in reference data with all returns (here, we are using a scan of a shed interior)
refData <- readLAS("D:/lidar/Shedtest1.las") 

#Read in normalized scan data
tlsData <- readLAS("D:/lidar/normalize/norm_clf_BLK360 Scan.las") 

#radii or distance from sensor to make calculations
radii <- c(2, 4, 6, 8, 10) 

height = 2 #height of the scanner above ground
thetaDiv=.5 #size of each cell in angle of theta (columns of raster)
phiDiv=1.5 #size of each cell in angle of phi (rows of raster)

#Define a blank raster grid to serve as a templage for creating and aligning all grids.
r_blank <- raster(ncol= 360/thetaDiv, nrow = 180/phiDiv, 
                  xmn=-180, xmx=180,ymn=-90,ymx=90)

# Declaring functions ----

cld_fun<-function(r,tlsData, type="Data"){
  #Aim: Calculate coordinates where each ray would intersect the sphere
  #Arguments:r = radius; tlsData: LAS object
  #Arguments:type: whether reference or data
  #Output: cld dataframe
  if (type=="Data") {  
  # create a dataframe to calculate distances to the scanner. For data, must subtract the scanner height from Z coordinate.
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
  } else if (type=="Ref") {
    # create a dataframe to calculate distances to the scanner. For references, do not need to subtract scanner height from Z coordinate. 
    cld <- as.data.frame(cbind(x=tlsData@data$X, y=tlsData@data$Y, z=tlsData@data$Z))
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
  } else {
      print("Invalid Type Defined.")
    }
  return(cld)    
}

sphere_func<-function(cld=cld, type="Data"){
  #Aim: calculate spherical coordinates
  #Arguments: cld: dataframe with Cartesian distances obtained with cld_fun
  # Arguments:type: whether reference or data
  #Output: list of dataframe objects s3 and p3. 
  if (type=="Data") {
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
    
    #Create table of x, y, z values where each point did or would have intersected the sphere
    allF <- bind_rows(just_sphere, just_grd)
    
    #Create table of x, y, z values where each point that passed the sphere intersected it
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
  } else if (type=="Ref") {
      #Same as above, bur for reference. Do not need to differentiate those that passed the sphere.
      s <- as.data.frame(cart2sph(x=cld$xs, 
                                  y=cld$ys, 
                                  z=cld$zs))
      p <- as.data.frame(cart2sph(x=cld$xs, 
                                  y=cld$ys, 
                                  z=cld$zs))
      
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
  } else {
    print("Invalid Type Defined.")
  }
  return(multipoint_list)
} 

pointCount <- function(r, pts){
  # Idea source: https://gis.stackexchange.com/questions/309407/computing-number-of-points-in-a-raster-grid-cell-in-r
  # Aim: Calculate the number of points in a raster cell grid. 
  # Argument: r: radius
  # Output: r2: raster object with point count in each cell
  
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

on_sphere_multi_v2 <- function(tlsData, refData, radii, height, thetaDiv=1, phiDiv=1){
  # Aim: Transforms a point cloud into a raster showing the points that passed or not pass a certain radius, 
  # depending on the layer of the raster.
  # Arguments: tlsData: cloudpoint, radii: list of radii, height: height of scanner.
  # Output: raster object with different layers, representing each radius. 
  for(ra in radii){ # Begin for radius loop 
    if(ra == radii[1]){ # Logic statement of 1st radius
      r <- ra #read radius
      # apply functions to obtain spherical coordinates for actual data
      list_point<-cld_fun(r=r,tlsData = tlsData, type="Data") %>% sphere_func(type="Data")
      # apply functions to obtain spherical coordinates for reference scan
      ref_points<-cld_fun(r=r, tlsData = refData, type="Ref") %>% sphere_func(type="Ref")
      
      #Obtain grid of point counts per cell for reference data
      gridRef <- pointCount(r_blank,ref_points$s3)
      #Obtain grid of point counts per cell for all data points on sphere
      gridS <- pointCount(r_blank, list_point$s3)
      #Obtain grid of point counts for only data points that passed sphere.
      gridP <- pointCount(r_blank, list_point$p3)
      #Get number of points that returned in voxel by subtracting the pass count from the total count
      gridIn <- gridS - gridP

      #If not returns in the cell in rerference data, make null
      null1 <- raster::calc(gridRef, function(x){x[x<=0]<-NA; return(x)})
      #Make a mask of all non-null cells
      null2 <- null1 >= 0 #??

      out <- gridIn/gridRef # Divide number of returns in voxel by total possible from reference to get a proportion of total returns. 
      out <- raster::calc(out, function(x){x[x<=0]<-0; return(x)}) # Truncate data to range of 0 to 1.
      out <- raster::calc(out, function(x){x[x>=1]<-1; return(x)}) # Truncate data to range of 0 to 1. 
      out <- out*null2 # Multiply result by mask to flag NULL cells with no returns. 
    }else if(ra == radii[2]){# logical statement 2nd radius
      r<-ra # read radius
      
      # apply functions to obtain spherical coordinates for actual data
      list_point<-cld_fun(r=r,tlsData = tlsData, type="Data") %>% sphere_func(type="Data") 
      #Make a copy of result for prior radius
      first <- out
      #If proportion for pior radius is 1, then make null since no returns are avaialble (occlusion)
      first2 <- raster::calc(first, function(x){x[x>=1]<-NA; return(x)})
      #Convert to mask
      mask3 <- first2 >= 0
      
      #Returns available = total - those that were returned in prior cell
      gridPre <- gridRef - gridIn
      #Make copy of prior count
      gridInPre <- gridIn
      
      #Obtain grid of point counts per cell for all data points on sphere
      gridS <- pointCount(r_blank, list_point$s3)
      #Obtain grid of point counts for only data points that passed sphere.
      gridP <- pointCount(r_blank, list_point$p3)
      #Number of returns in voxel = All - Pass - Returns from prior voxel
      gridIn <- gridS - gridP - gridInPre#create point
      
      #Get total number of returns that have returned (Add all returns in current and prior voxels)
      total_re <- gridIn + gridInPre
    
      out <- ((gridIn)/gridPre) # Proportion = count in voxel/total that passed through voxel
      out <- raster::calc(out, function(x){x[x<= 0]<-0; return(x)}) # Truncate to range 0 to 1.
      out <- raster::calc(out, function(x){x[x>= 1]<-1; return(x)}) # Truncate to range 0 to 1. 
      out <- out*mask3 #Mask for occlusion 
      out2 <- stack(first, out) # stack raster results
    }else{ # logic statement n_th radius (All radii after the first two)
      r<-ra # read radius
      # apply functions to obtain spherical coordinates for actual data
      list_point<-cld_fun(r=r,tlsData = tlsData, type="Data") %>% sphere_func(type="Data") 
      
      #Copy preceeding raster stack
      preceed <- out2
      #Copy Prior result
      prior <- out
      #Make a occulusion mask
      prior2 <- raster::calc(prior, function(x){x[x>=1]<-NA; return(x)})
      mask3 <- prior2 >= 0
      
      #Total possible = all pulses - poor returns 
      gridPre <- gridRef - total_re
      
      #Obtain grid of point counts per cell for all data points on sphere
      gridS <- pointCount(r_blank, list_point$s3)
      #Obtain grid of point counts for only data points that passed sphere
      gridP <- pointCount(r_blank, list_point$p3)
      #Total in sphere = all - pass - prior returns
      gridIn <- gridS - gridP - total_re
      
      #Make new copy of total by adding prior total with voxel being processes.
      total_re <- total_re + gridIn
      
      #proportion = total returns from voxel/proportion 
      out <- ((gridIn)/gridPre)
      out <- raster::calc(out, function(x){x[x<=0]<-0; return(x)}) #Truncate to range 0 to 1.
      out <- raster::calc(out, function(x){x[x>=1]<-1; return(x)}) #Truncate to range 0 to 1.
      out <- out*mask3 #Apply maks for occlusion
      out2 <- stack(preceed, out) # stack results
    }
  }#Close radius loop
  names(out2) <- as.character(radii)# name layers
  return(out2) # set function output. 
} 

# Applying functions ----

test <- on_sphere_multi_v2(tlsData=tlsData, refData=refData, radii=radii, height=height, thetaDiv=thetaDiv, phiDiv=phiDiv)

writeRaster(test, "D:/sphere_vox2.tif") # save data. 
