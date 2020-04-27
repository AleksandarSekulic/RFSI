library(sp)
library(meteo) # sa rforge
library(plotGoogleMaps)
library(ggplot2)
library(nabor)
library(ranger)
library(rgdal)
library(raster)
library(doParallel)
library(gstat)
library(plyr)
library(caret)
library(CAST)
library(ggpubr)
library(automap)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(rasterVis)
library(RSAGA)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

utm33 <- "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

border <- readOGR("borders/osm/union_of_selected_boundaries_AL2-AL2.shp")
border <- spTransform(border, CRS(utm33))

dir.create("predictions")

v = "TEMP"
year="2008"
# time=seq(as.POSIXct(paste(year,"-01-01", sep="")), as.POSIXct(paste(year,"-12-31", sep="")), by="day")
# days<-gsub("-","",time,fixed=TRUE)
# daysNum <- length(time)
covariates <- c("TEMP", "HRdem","HRdsea","Lat","Lon","HRtwi","INSOL","ctd","MODIS.LST")

setwd(paste(wd, "/temp_data/", sep = ""))
load(file = "stfdf.rda")

load(file = "pc.HRlocsMODIS.f.rda")
load(file = "vgm.st.rda")
load(file = "grids.rda")
load(file = "LSTdate.rda")
load(file = "LST.listday.rda")
load(file = "tscale.rda")
load(file = "pca.rda")
load(file = "lms.pc.rda")

setwd(paste(wd, "/predictions/", sep = ""))
dir.create("strk")
dir.create("idw")
dir.create("rf")
dir.create("rfsi")
dir.create("rfsp")

load("../models/RF.rda")
load("../models/RFSI.rda")
load("../models/RFsp.rda")

g.TEMP <- gstat(id=c("rTEMP2008"), formula=rTEMP2008~1, data=pc.HRlocsMODIS.f, nmin=20, nmax=50, model=vgm.st) #, set=list(cn_max=1e10))

# available MODIS images:
slices <- c(5,13,20,28)
# new locations:
grids.new <- grids["mask"]
grids.xy <- as(grids[c("HRdem", "HRdsea", "Lat", "Lon", "HRtwi", "mask")], "SpatialPixelsDataFrame")
theta <- min(pc.HRlocsMODIS.f$cday)

###### Predictions ######

for(i in slices) {  
  ### STRK ###
  newlocs.xyt <- data.frame(grids.xy)
  day <- LSTdate[i]
  slice <- floor(unclass(as.POSIXct(LSTdate[i]))/86400)[[1]]
  newlocs.xyt$cday <- rep(slice, length(newlocs.xyt[1]))
  newlocs.xyt$cdays <- tscale * newlocs.xyt$cday
  # insolation:
  newlocs.xyt$INSOL <- grids@data[grids.xy@grid.index, paste("INSOL", slice-round(unclass(as.POSIXct("2007-12-31"))/86400, 0)[1], sep="")]
  # MODIS LST:
  LSTname <- strsplit(strsplit(LST.listday[i], "LST/")[[1]][2], ".LST_")[[1]][1]
  # filter missing pixels step 1:
  LSTname1M <- strsplit(strsplit(LST.listday[i-1], "LST/")[[1]][2], ".LST_")[[1]][1]
  LSTname1P <- strsplit(strsplit(LST.listday[i+1], "LST/")[[1]][2], ".LST_")[[1]][1]
  grids.new$tmp <- ifelse(is.na(grids@data[,LSTname]), (grids@data[,LSTname1M]+grids@data[,LSTname1P])/2, grids@data[,LSTname])
  gc()
  # filter missing pixels step 2:
  writeGDAL(grids.new["tmp"], "MODIS_tmp.sdat", "SAGA")
  rsaga.geoprocessor(lib="grid_tools", module=7, param=list(INPUT="MODIS_tmp.sgrd", MASK="mask.sgrd", RESULT="MODIS_f.sgrd"))
  newlocs.xyt$MODIS.LST <- readGDAL("MODIS_f.sdat", silent=TRUE)$band1[grids.xy@grid.index]
  # Extract Principal Components (using the existing "pca" object):
  pc.newlocs <- predict(pca, newlocs.xyt)
  pc.newlocs.xyt <- cbind(data.frame(pc.newlocs), newlocs.xyt[,c("mask","cday","cdays","x","y")]) 
  # convert to grids:
  coordinates(pc.newlocs.xyt) <- c("x","y","cday")
  proj4string(pc.newlocs.xyt) <- CRS(proj4string(grids))
  # TEMP.ok <- predict.gstat(g.TEMP, pc.newlocs.xyt, debug.level=-1)  # takes 5 mins per map! Possible singular matrix problems;
  TEMP.ok <- predict(g.TEMP, pc.newlocs.xyt, debug.level=-1)  # takes 5 mins per map! Possible singular matrix problems;
  TEMP.reg <- predict(lms.pc, pc.newlocs.xyt)
  pc.newlocs.xyt@data[,"strk"] <- TEMP.reg + TEMP.ok$rTEMP2008.pred
  pc.newlocs.xyt@data[,"strk_var"] <- TEMP.ok$rTEMP2008.var
  pc.newlocs.xyt <- data.frame(pc.newlocs.xyt)
  coordinates(pc.newlocs.xyt) <- c("x","y")
  gridded(pc.newlocs.xyt) <- TRUE
  p <- raster( as(pc.newlocs.xyt["strk"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  p <- p * 100
  writeRaster(p, paste("strk/", LSTname, sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  p <- raster( as(pc.newlocs.xyt["strk_var"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  p <- round(1.349 * sqrt(p) * 100)
  writeRaster(p, paste("strk/", LSTname, "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  
  ### IDW ###
  p = 1.8
  n = 11
  data <- stfdf@sp@data
  st <- stfdf[, day]
  st@data <- cbind(st@data, data)
  st.all <- st
  st <- st[!is.na(st$TEMP), ]
  gs <- gstat(formula=TEMP~1, locations=st, nmax = n, set=list(idp = p))
  pc.newlocs.xyt@proj4string <- CRS(utm33)
  pred <- predict(gs, pc.newlocs.xyt)
  pc.newlocs.xyt$idw <- pred$var1.pred
  p <- raster( as(pc.newlocs.xyt["idw"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  p <- p * 100
  writeRaster(p, paste("idw/", LSTname, sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  
  ### RF ###
  newlocs.xyt$ctd <- cos((newlocs.xyt$cday-theta)*pi/180)
  pre <- predict(rf_model, newlocs.xyt)
  pc.newlocs.xyt$rf <- pre$predictions
  p <- raster( as(pc.newlocs.xyt["rf"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  p <- p * 100
  writeRaster(p, paste("rf/", LSTname, sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  pred_int <- 0.5# 0.6827 # one sigma
  quantiles = c((1-pred_int)/2, 0.5, 1-(1-pred_int)/2)
  quant_reg <- predict(rf_model, newlocs.xyt, type = "quantiles", quantiles = quantiles)
  iqr <- (quant_reg$predictions[,3] - quant_reg$predictions[,1]) / 2
  iqr <- round(iqr * 100)
  pc.newlocs.xyt$rf_iqr <- iqr
  p <- raster( as(pc.newlocs.xyt["rf_iqr"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  writeRaster(p, paste("rf/", LSTname, "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  
  ### RFSI ###
  loc <- newlocs.xyt
  coordinates(loc) <- c("x", "y")
  n_obs = 10
  nearest_obs <- near.obs(
    locations = loc,
    observations = st,
    zcol = "TEMP",
    n.obs = n_obs
  )
  df <- cbind(newlocs.xyt, nearest_obs)
  pre <- predict(rfsi_model, df)
  pc.newlocs.xyt$rfsi <- pre$predictions
  p <- raster( as(pc.newlocs.xyt["rfsi"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  p <- p * 100
  writeRaster(p, paste("rfsi/", LSTname, sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  pred_int <- 0.5# 0.6827 # one sigma
  quantiles = c((1-pred_int)/2, 0.5, 1-(1-pred_int)/2)
  quant_reg <- predict(rfsi_model, df, type = "quantiles", quantiles = quantiles)
  iqr <- (quant_reg$predictions[,3] - quant_reg$predictions[,1]) / 2
  iqr <- round(iqr * 100)
  pc.newlocs.xyt$rfsi_iqr <- iqr
  p <- raster( as(pc.newlocs.xyt["rfsi_iqr"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  writeRaster(p, paste("rfsi/", LSTname, "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  
  ### RFsp ###
  
  ### Create BBOX ###
  # xmin=min(gg_po@coords[, "x"])-1000; xmax=max(gg_po@coords[, "x"])+1000;
  # ymin=min(gg_po@coords[, "y"])-1000 ; ymax=max(gg_po@coords[, "y"])+1000 # Serbia (proj bbox is 166021.4431, 833978.5569, 0.0000, 9329005.1825)
  # bbox = extent(xmin, xmax, ymin, ymax)
  # bbox <- as(bbox, "SpatialPolygons")
  # bbox@proj4string = CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  # fishnet <- raster(bbox)
  # res(fishnet) <- c(1000, 1000) # 1 x 1 km
  # fishnet[] <- runif(length(fishnet), -10, 10)
  # fishnet = as(fishnet, "SpatialPixelsDataFrame")
  # # plot(fishnet)
  # 
  # st <- as(fishnet, "SpatialPoints")
  fishnet = as(pc.newlocs.xyt, "SpatialPixelsDataFrame")
  
  ## Buffer distances only:
  grid.distP <- GSIF::buffer.dist(st.all["staid"], fishnet[1], as.factor(1:nrow(st.all)))
  ov_toma <- over(pc.newlocs.xyt, grid.distP)
  df <- cbind(newlocs.xyt, ov_toma)
  pre <- predict(rfsp_model, df)
  pc.newlocs.xyt$rfsp <- pre$predictions
  p <- raster( as(pc.newlocs.xyt["rfsp"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  p <- p * 100
  writeRaster(p, paste("rfsp/", LSTname, sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  pred_int <- 0.5# 0.6827 # one sigma
  quantiles = c((1-pred_int)/2, 0.5, 1-(1-pred_int)/2)
  quant_reg <- predict(rfsp_model, df, type = "quantiles", quantiles = quantiles)
  iqr <- (quant_reg$predictions[,3] - quant_reg$predictions[,1]) / 2
  iqr <- round(iqr * 100)
  pc.newlocs.xyt$rfsp_iqr <- iqr
  p <- raster( as(pc.newlocs.xyt["rfsp_iqr"], "SpatialPixelsDataFrame"))
  p = crop(p, border)
  p = mask(p, border)
  writeRaster(p, paste("rfsp/", LSTname, "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  
}


###### Prediction plots###################################################################

load('../temp_data/stfdf.rda')
# stfdf <- stfdf[!is.na(stfdf[, 1]$TEMP), ]

summary(stfdf[ , LSTdate[slices[1]], "TEMP"])
summary(stfdf[ , LSTdate[slices[2]], "TEMP"])
summary(stfdf[ , LSTdate[slices[3]], "TEMP"])
summary(stfdf[ , LSTdate[slices[4]], "TEMP"])

met <- c("idw", "rf", "rfsi", "rfsp") # "strk", 

# paste("rfsi/", LSTname, "_iqr", sep = "")
LSTname <- c("LST2008_02_02", "LST2008_04_06", "LST2008_06_01", "LST2008_08_04")
time <- LSTdate[slices]

pred <- raster(paste(met[1], "/",
                     LSTname[1], ".tif", sep = "")) / 100
pred <- projectRaster(pred, crs=CRS("+proj=longlat +datum=WGS84"))
pred <- as(pred, "SpatialPixelsDataFrame")
names(pred) <- paste(met[1], 1, sep = "")

for (t in 1:length(time)){
  for (m in 1:length(met)){
    r <- raster(paste(met[m], "/", LSTname[t], ".tif", sep = "")) / 100
    # r <- projectRaster(r, crs=CRS("+proj=longlat +datum=WGS84"))
    r <- as(r, "SpatialPixelsDataFrame")
    names(r) <- met[m]
    pred[[paste(met[m], t, sep = "")]] <- r[[met[m]]]
  }
}

### levelplot ###

for (t in 1:length(time)){
  for (m in 1:length(met)){
    r <- raster(paste(met[m], "/", LSTname[t], ".tif", sep = "")) / 100
    r <- projectRaster(r, crs=CRS("+proj=longlat +datum=WGS84"))
    assign(paste(met[m], t, sep = ""), r)
  }
}

summary(stfdf[, time[1], "TEMP"])
# summary(pred@data$strk1)
summary(pred@data$idw1)
summary(pred@data$rf1)
summary(pred@data$rfsi1)
summary(pred@data$rfsp1)

max(pred@data) # 30.43
max(pred@data$idw4)
max(pred@data$rf4)
max(pred@data$rfsp4)
max(pred@data$rfsi4)
# max(pred@data$strk4)
min(pred@data) # -3.95

# colorRampPalette(c("blue", "deepskyblue3","white", "orange", "red"))(60)
# pal <- colorRampPalette(brewer.pal(9, "Blues"))
# pal <- c(colorRampPalette(c("blue", "deepskyblue3","white"))(3)[-3], colorRampPalette(c("white", "orange", "red"))(13)[-1])
# cut <- c(floor(min(pred@data)), -2, 0, 2, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, ceiling(max(pred@data)))
# pal <- colorRampPalette(c("blue", "deepskyblue3","white", "orange", "red"))(7)
# cut <- c(floor(min(pred@data)), 0, 2, 4, 6, 8, 10, ceiling(max(pred@data)))
# cut2 <- cut
# cut2 <- c(floor(min(pred@data)), -2, 0, 2, 5, 10, 15, 20, 25, ceiling(max(pred@data)))

pal <- c(colorRampPalette(c("blue", "deepskyblue3","white"))(2)[-2], colorRampPalette(c("white", "orange", "red"))(15)[-1])
cut <- c(floor(min(stfdf[, time[1]]$TEMP, na.rm=T)), 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, ceiling(max(stfdf[, time[1]]$TEMP, na.rm=T)))
cut2 <- c(floor(min(stfdf[, time[1]]$TEMP, na.rm=T)), 0, 2, 4, 6, 8, 10, 12, ceiling(max(stfdf[, time[1]]$TEMP, na.rm=T)))

s <- stack(# strk1, strk2, strk3, strk4,
           idw1,# idw2, idw3, idw4,
           rf1,# rf2, rf3, rf4,
           rfsp1,# rfsp2, rfsp3, rfsp4,
           rfsi1)#, rfsi2, rfsi3, rfsi4)

# s <- stack(strk1, idw1, rf1, rfsi1, rfsp1)

met2 <- c("IDW", "RF", "RFsp", "RFSI") # "STRK", 

time_max <- c()
for (i in 1){#:4){
  time_max <- c(time_max, paste("max: ", round(max(stfdf[ , time[i]]$TEMP, na.rm=T), 1), "°C       ", #\n",
                                "min: ", round(min(stfdf[ , time[i]]$TEMP, na.rm=T), 1), "°C\n",
                                sep = ""))
}

# tiff("../plot/predictions.tiff", width = 100, height = 120, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/predictions.jpeg", width = 100, height = 120, units = 'mm', res = 600)
levelplot(s, 
          margin=FALSE,                       
          colorkey=list(
            space='bottom',                   
            labels=list(at=cut2,cex=0.7),
            axis.line=list(col='black'),
            width=0.75
          ),
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(x=list(color="black", draw=TRUE, alternating=c(3), cex=0.6, abbreviate=T, minlength=1, rot=45),
                      y=list(color="black", draw=TRUE, alternating=c(3), cex=0.6, abbreviate=T, minlength=1, rot=45)),            
          col.regions=pal,
          at=cut,
          names.attr=met2,#rep('', nlayers(s)),
          par.strip.text=list(cex=0.75), #, lines=1, col="blue, fontfamily='Serif'),
          # xlab.top=list(as.character(time), cex=0.5, fontface='bold'),
          xlab=list(as.character(time_max), cex=0.7),
          ylab = NULL,
          # ylab=list(rev(met2), cex=0.5, fontface='bold'),
          layout=c(2, 2),
          # layout=c(4, 4),
          # layout=c(4, 5),
          # layout=c(3, 2),
          as.table = TRUE) +
  # layer(sp.points(stfdf@sp, pch=3, col="black", cex=0.1, alpha=0.5)) +
  layer(sp.polygons(border, lwd=0.1)) +
  layer({
    SpatialPolygonsRescale(layout.north.arrow(type = 2),
                           offset = c(18.2,43.8),
                           scale = 1)
  }, packets = 2) +
  layer({
    grid.text(x= 17.7, y=44.3, "N", gp=gpar(cex=1.6), rot=0,
              default.units='native')
  }, packets = 2)+
  layer({
    xs <- seq(14.4, 16.1, by=0.85)
    grid.rect(x=xs[-3], y=42.8,
              width=0.85, height=0.14,
              gp=gpar(fill=c('black', 'transparent')),
              default.units='native')
    # x = xs - 40000
    xs[3] <- xs[3] + 0.3
    grid.text(x= xs - 0.4, y=42.5, c(0, 75, "150km"),# seq(0, 150, by=75),
              gp=gpar(cex=0.6), rot=0,
              default.units='native')
  }, packets = 4)
dev.off()
#
