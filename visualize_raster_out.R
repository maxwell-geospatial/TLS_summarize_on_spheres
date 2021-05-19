zR1 <- test$X2 == 0
zR1[zR1 == 0] <- NA
zR2 <- test$X4 == 0
zR2[zR2 == 0] <- NA
zR3 <- test$X6 == 0
zR3[zR3 == 0] <- NA
zR4 <- test$X8 == 0
zR4[zR4 == 0] <- NA
zR5 <- test$X10 == 0
zR5[zR5 == 0] <- NA

library(tmap)
plotR1 <- tm_shape(test$X2)+
  tm_raster(palette="Greens", colorNA="lightblue")+
  tm_shape(zR1)+
  tm_raster(palette="gray")
plotR2 <- tm_shape(test$X4)+
  tm_raster(palette="Greens", colorNA="lightblue")+
  tm_shape(zR2)+
  tm_raster(palette="gray")
plotR3 <- tm_shape(test$X6)+
  tm_raster(palette="Greens", colorNA="lightblue")+
  tm_shape(zR3)+
  tm_raster(palette="gray")
plotR4 <- tm_shape(test$X8)+
  tm_raster(palette="Greens", colorNA="lightblue")+
  tm_shape(zR4)+
  tm_raster(palette="gray")
plotR5 <- tm_shape(test$X10)+
  tm_raster(palette="Greens", colorNA="lightblue")+
  tm_shape(zR5)+
  tm_raster(palette="gray")
tmap_arrange(plotR1, plotR2, plotR3, plotR4, plotR5)


library(tmap)
plotR1 <- tm_shape(test$X2)+
  tm_raster(palette="Greens", colorNA="red")
plotR2 <- tm_shape(test$X4)+
  tm_raster(palette="Greens", colorNA="red")
plotR3 <- tm_shape(test$X6)+
  tm_raster(palette="Greens", colorNA="red")
plotR4 <- tm_shape(test$X8)+
  tm_raster(palette="Greens", colorNA="red")
plotR5 <- tm_shape(test$X10)+
  tm_raster(palette="Greens", colorNA="red")
tmap_arrange(plotR1, plotR2, plotR3, plotR4, plotR5)

head <- read.lasheader("D:/lidar/Shedtest1.las")
sheddf <- cld_fun(10, refData, type="Ref")
sheddf2 <- sheddf[,16:18]
names(sheddf2) <- c("X", "Y", "Z")
write.las("D:/lidar/ShedSphere.las", head, sheddf2)


head <- read.lasheader("D:/lidar/normalize/norm_clf_BLK360 Scan.las")
tlsdf <- cld_fun(10, tlsData, type="Data")
tlsdf2 <- tlsdf[,16:18]
names(tlsdf2) <- c("X", "Y", "Z")
write.las("D:/lidar/tlsSphere.las", head, tlsdf2)
