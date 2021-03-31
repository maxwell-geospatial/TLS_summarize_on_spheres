#Read in normalized LAS data
tlsData <- readLAS("D:/lidar/normalize/norm_clf_BLK360 Scan.las")
#Define radius of sphere
radii <- c(2, 4, 6, 8, 10) 
#Create points in the 3D space
x <- tlsData@data$X
y <- tlsData@data$Y
z <- tlsData@data$Z #Subtract 2 for sensor height

height = 2
thetaDiv=.5
phiDiv=.5
grdDist=.25

on_sphere_multi <- function(tlsData, radii, height, thetaDiv=1, phiDiv=1, grdDist=.5){
  for(r in radii){
    if(r==radii[1]){
      if(r < height){
        r <- r
        x <- tlsData@data$X
        y <- tlsData@data$Y
        z <- tlsData@data$Z - height
        cld <- as.data.frame(cbind(x=x, y=y, z=z))
        
        #Calculate coordinates where each ray could intersect the sphere
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
        just_sphere <- cld
        
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
        npass_sphere <- just_sphere %>% filter(dist_act < dist_sphere)
        
        #Create table of x, y, z values where each point did or would have intersected the ground
        just_sphere$xf <- just_sphere$xs 
        just_sphere$yf <- just_sphere$ys 
        just_sphere$zf <- just_sphere$zs 

        
        #Create table of x, y, z values where each point that passed the ground intersected it
        pass_sphere$xf <- pass_sphere$xs 
        pass_sphere$yf <- pass_sphere$ys 
        pass_sphere$zf <- pass_sphere$zs 

        
        #Create table of x, y, z values of returns within sphere
        npass_sphere$xf <- npass_sphere$x 
        npass_sphere$yf <- npass_sphere$y 
        npass_sphere$zf <- npass_sphere$z 

        
        #Convert from Cartesian to spherical coordinates
        s <- as.data.frame(cart2sph(x=just_sphere$xf, 
                                    y=just_sphere$yf, 
                                    z=just_sphere$zf))
        p <- as.data.frame(cart2sph(x=pass_sphere$xf, 
                                    y=pass_sphere$yf, 
                                    z=pass_sphere$zf))
        np <- as.data.frame(cart2sph(x=npass_sphere$xf, 
                                     y=npass_sphere$yf,
                                     z=npass_sphere$zf))
        
        s2 <- s[,c(1,2)]
        p2 <- p[,c(1,2)]
        np2 <- np[,c(1,2)]
        
        s2$theta <- s$theta*(180/pi)
        s2$phi <- s$phi*(180/pi)
        p2$theta <- p$theta*(180/pi)
        p2$phi <- p$phi*(180/pi)
        np2$theta <- np$theta*(180/pi)
        np2$phi <- np$phi*(180/pi)
        
        r_blank <- raster(ncol= 360/thetaDiv, 
                          nrow = 180/phiDiv, 
                          xmn=-180, xmx=180, 
                          ymn=0, 
                          ymx=90)
        
        s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
        p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
        np3 <- st_multipoint(x=as.matrix(np2), dim="XY")
        
        gridS <- pointCount(r_blank, s3)
        gridP <- pointCount(r_blank, p3)
        gridNP <- pointCount(r_blank, np3)
      
        outLst <- as.list(as.list(gridS, gridP))
        
      } else{
        r <- r
        x <- tlsData@data$X
        y <- tlsData@data$Y
        z <- tlsData@data$Z - height
        cld <- as.data.frame(cbind(x=x, y=y, z=z))
        
        #Calculate coordinates where each ray could intersect the sphere
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
        just_grd <- setdiff(cld, just_sphere)
        
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
        npass_sphere <- just_sphere %>% filter(dist_act < dist_sphere)
        
        #Determine coordinates where points would intersect the ground
        just_grd$grd_z <- -height
        just_grd$t <- just_grd$grd_z/just_grd$z
        just_grd$grd_x <- just_grd$t*just_grd$x
        just_grd$grd_y <- just_grd$t*just_grd$y
        just_grd$dist_act <- sqrt(just_grd$x^2 + just_grd$y^2 + just_grd$z^2)
        just_grd$dist_grd <- sqrt(just_grd$grd_x^2 + just_grd$grd_y^2 + (just_grd$grd_z)^2)
        npass_grd <- just_grd %>% filter(dist_act < dist_grd)
        pass_grd <- just_grd %>% filter(dist_act >= dist_grd)
        
        #Create table of x, y, z values where each point did or would have intersected the ground
        just_sphere$xf <- just_sphere$xs 
        just_sphere$yf <- just_sphere$ys 
        just_sphere$zf <- just_sphere$zs 
        just_grd$xf <- just_grd$grd_x
        just_grd$yf <- just_grd$grd_y
        just_grd$zf <- just_grd$grd_z
        
        #Create table of x, y, z values where each point that passed the ground intersected it
        pass_sphere$xf <- pass_sphere$xs 
        pass_sphere$yf <- pass_sphere$ys 
        pass_sphere$zf <- pass_sphere$zs 
        pass_grd$xf <- pass_grd$grd_x
        pass_grd$yf <- pass_grd$grd_y
        pass_grd$zf <- pass_grd$grd_z
        
        #Create table of x, y, z values of returns within sphere
        npass_sphere$xf <- npass_sphere$x 
        npass_sphere$yf <- npass_sphere$y 
        npass_sphere$zf <- npass_sphere$z 
        npass_grd$xf <- npass_grd$x
        npass_grd$yf <- npass_grd$y
        npass_grd$zf <- npass_grd$z
        
        #Convert from Cartesian to spherical coordinates
        s <- as.data.frame(cart2sph(x=just_sphere$xf, 
                                    y=just_sphere$yf, 
                                    z=just_sphere$zf))
        p <- as.data.frame(cart2sph(x=pass_sphere$xf, 
                                    y=pass_sphere$yf, 
                                    z=pass_sphere$zf))
        np <- as.data.frame(cart2sph(x=npass_sphere$xf, 
                                     y=npass_sphere$yf,
                                     z=npass_sphere$zf))
        
        s2 <- s[,c(1,2)]
        p2 <- p[,c(1,2)]
        np2 <- np[,c(1,2)]
        
        s2$theta <- s$theta*(180/pi)
        s2$phi <- s$phi*(180/pi)
        p2$theta <- p$theta*(180/pi)
        p2$phi <- p$phi*(180/pi)
        np2$theta <- np$theta*(180/pi)
        np2$phi <- np$phi*(180/pi)
        
        r_blank <- raster(ncol= 360/thetaDiv, 
                          nrow = 180/phiDiv, 
                          xmn=-180, xmx=180, 
                          ymn=0, 
                          ymx=90)
        
        s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
        p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
        np3 <- st_multipoint(x=as.matrix(np2), dim="XY")
        
        gridS <- pointCount(r_blank, s3)
        gridP <- pointCount(r_blank, p3)
        gridNP <- pointCount(r_blank, np3)
        
        gs <- just_grd[,c("xf", "yf")]
        gp <- pass_grd[,c("xf", "yf")]
        gnp <- npass_grd[,c("xf", "yf")]
        
        gs2 <- st_multipoint(x=as.matrix(gs), dim="XY")
        gp2 <- st_multipoint(x=as.matrix(gp), dim="XY")
        gnp2 <- st_multipoint(x=as.matrix(gnp), dim="XY")
        
        r_blank2 <- raster(ncol= r/grdDist, 
                           nrow = r/grdDist, 
                           xmn=-r, xmx=r, 
                           ymn=-r, 
                           ymx=r)
        
        gridSg <- pointCount(r_blank2, gs2)
        gridPg <- pointCount(r_blank2, gp2)
        gridNPg <- pointCount(r_blank2, gnp2)
        
        center <- st_sfc(st_point(x=c(0, 0), dim="XY"))
        circle <- as(st_buffer(center, r), Class="Spatial")
        gridSg2 <- mask(gridSg, circle)
        gridPg2 <- mask(gridPg, circle)
        gridNPg2 <- mask(gridNPg, circle)
        
        outLst <- as.list(as.list(gridS, gridP, gridSg2, gridSg2)) 
      }
    }else{
      if(r < height){
        r <- r
        x <- tlsData@data$X
        y <- tlsData@data$Y
        z <- tlsData@data$Z - height
        cld <- as.data.frame(cbind(x=x, y=y, z=z))
        
        #Calculate coordinates where each ray could intersect the sphere
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
        just_sphere <- cld
        
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
        npass_sphere <- just_sphere %>% filter(dist_act < dist_sphere)
        
        #Create table of x, y, z values where each point did or would have intersected the ground
        just_sphere$xf <- just_sphere$xs 
        just_sphere$yf <- just_sphere$ys 
        just_sphere$zf <- just_sphere$zs 

        #Create table of x, y, z values where each point that passed the ground intersected it
        pass_sphere$xf <- pass_sphere$xs 
        pass_sphere$yf <- pass_sphere$ys 
        pass_sphere$zf <- pass_sphere$zs 

        #Create table of x, y, z values of returns within sphere
        npass_sphere$xf <- npass_sphere$x 
        npass_sphere$yf <- npass_sphere$y 
        npass_sphere$zf <- npass_sphere$z 

        
        #Convert from Cartesian to spherical coordinates
        s <- as.data.frame(cart2sph(x=just_sphere$xf, 
                                    y=just_sphere$yf, 
                                    z=just_sphere$zf))
        p <- as.data.frame(cart2sph(x=pass_sphere$xf, 
                                    y=pass_sphere$yf, 
                                    z=pass_sphere$zf))
        np <- as.data.frame(cart2sph(x=npass_sphere$xf, 
                                     y=npass_sphere$yf,
                                     z=npass_sphere$zf))
        
        s2 <- s[,c(1,2)]
        p2 <- p[,c(1,2)]
        np2 <- np[,c(1,2)]
        
        s2$theta <- s$theta*(180/pi)
        s2$phi <- s$phi*(180/pi)
        p2$theta <
        np2$theta <- np$theta*(180/pi)
        np2$phi <- np$phi*(180/pi)
        
        r_blank <- raster(ncol= 360/thetaDiv, 
                          nrow = 180/phiDiv, 
                          xmn=-180, xmx=180, 
                          ymn=0, 
                          ymx=90)
        
        s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
        p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
        np3 <- st_multipoint(x=as.matrix(np2), dim="XY")
        
        gridS <- pointCount(r_blank, s3)
        gridP <- pointCount(r_blank, p3)
        gridNP <- pointCount(r_blank, np3)
        
        outLst <- c(outLst, as.list(gridS, gridP))

      } else{
        r <- r
        x <- tlsData@data$X
        y <- tlsData@data$Y
        z <- tlsData@data$Z - height
        cld <- as.data.frame(cbind(x=x, y=y, z=z))
        
        #Calculate coordinates where each ray could intersect the sphere
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
        just_grd <- setdiff(cld, just_sphere)
        
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
        npass_sphere <- just_sphere %>% filter(dist_act < dist_sphere)
        
        #Determine coordinates where points would intersect the ground
        just_grd$grd_z <- -height
        just_grd$t <- just_grd$grd_z/just_grd$z
        just_grd$grd_x <- just_grd$t*just_grd$x
        just_grd$grd_y <- just_grd$t*just_grd$y
        just_grd$dist_act <- sqrt(just_grd$x^2 + just_grd$y^2 + just_grd$z^2)
        just_grd$dist_grd <- sqrt(just_grd$grd_x^2 + just_grd$grd_y^2 + (just_grd$grd_z)^2)
        npass_grd <- just_grd %>% filter(dist_act < dist_grd)
        pass_grd <- just_grd %>% filter(dist_act >= dist_grd)
        
        #Create table of x, y, z values where each point did or would have intersected the ground
        just_sphere$xf <- just_sphere$xs 
        just_sphere$yf <- just_sphere$ys 
        just_sphere$zf <- just_sphere$zs 
        just_grd$xf <- just_grd$grd_x
        just_grd$yf <- just_grd$grd_y
        just_grd$zf <- just_grd$grd_z
        
        #Create table of x, y, z values where each point that passed the ground intersected it
        pass_sphere$xf <- pass_sphere$xs 
        pass_sphere$yf <- pass_sphere$ys 
        pass_sphere$zf <- pass_sphere$zs 
        pass_grd$xf <- pass_grd$grd_x
        pass_grd$yf <- pass_grd$grd_y
        pass_grd$zf <- pass_grd$grd_z
        
        #Create table of x, y, z values of returns within sphere
        npass_sphere$xf <- npass_sphere$x 
        npass_sphere$yf <- npass_sphere$y 
        npass_sphere$zf <- npass_sphere$z 
        npass_grd$xf <- npass_grd$x
        npass_grd$yf <- npass_grd$y
        npass_grd$zf <- npass_grd$z
        
        #Convert from Cartesian to spherical coordinates
        s <- as.data.frame(cart2sph(x=just_sphere$xf, 
                                    y=just_sphere$yf, 
                                    z=just_sphere$zf))
        p <- as.data.frame(cart2sph(x=pass_sphere$xf, 
                                    y=pass_sphere$yf, 
                                    z=pass_sphere$zf))
        np <- as.data.frame(cart2sph(x=npass_sphere$xf, 
                                     y=npass_sphere$yf,
                                     z=npass_sphere$zf))
        
        s2 <- s[,c(1,2)]
        p2 <- p[,c(1,2)]
        np2 <- np[,c(1,2)]
        
        s2$theta <- s$theta*(180/pi)
        s2$phi <- s$phi*(180/pi)
        p2$theta <- p$theta*(180/pi)
        p2$phi <- p$phi*(180/pi)
        np2$theta <- np$theta*(180/pi)
        np2$phi <- np$phi*(180/pi)
        
        r_blank <- raster(ncol= 360/thetaDiv, 
                          nrow = 180/phiDiv, 
                          xmn=-180, xmx=180, 
                          ymn=0, 
                          ymx=90)
        
        s3 <- st_multipoint(x=as.matrix(s2), dim="XY")
        p3 <- st_multipoint(x=as.matrix(p2), dim="XY")
        np3 <- st_multipoint(x=as.matrix(np2), dim="XY")
        
        gridS <- pointCount(r_blank, s3)
        gridP <- pointCount(r_blank, p3)
        gridNP <- pointCount(r_blank, np3)
        
        gs <- just_grd[,c("xf", "yf")]
        gp <- pass_grd[,c("xf", "yf")]
        gnp <- npass_grd[,c("xf", "yf")]
        
        gs2 <- st_multipoint(x=as.matrix(gs), dim="XY")
        gp2 <- st_multipoint(x=as.matrix(gp), dim="XY")
        gnp2 <- st_multipoint(x=as.matrix(gnp), dim="XY")
        
        r_blank2 <- raster(ncol= r/grdDist, 
                           nrow = r/grdDist, 
                           xmn=-r, xmx=r, 
                           ymn=-r, 
                           ymx=r)
        
        gridSg <- pointCount(r_blank2, gs2)
        gridPg <- pointCount(r_blank2, gp2)
        gridNPg <- pointCount(r_blank2, gnp2)
        
        center <- st_sfc(st_point(x=c(0, 0), dim="XY"))
        circle <- as(st_buffer(center, r), Class="Spatial")
        gridSg2 <- mask(gridSg, circle)
        gridPg2 <- mask(gridPg, circle)
        gridNPg2 <- mask(gridNPg, circle)
        
        outLst <- c(outLst, as.list(gridS, gridP, gridSg2, gridSg2))
    
      }
    }
  }
  names(outLst) <- as.character(radii)
  return(outLst)
}

test <- on_sphere_multi(tlsData, radii, height, thetaDiv=1, phiDiv=1, grdDist=.5)
