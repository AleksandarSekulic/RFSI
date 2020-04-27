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

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

# lonmin=0 ;lonmax=4 ; latmin=40 ;latmax=43 # Ctalonia
border <- readOGR("catalonia_border/union_of_selected_boundaries_AL4-AL4.shp")

setwd(paste(wd, "/temp_data/", sep = ""))

v = "prcp"
years = 2016:2018
# time<-seq(as.Date("2016-01-01"), as.Date("2016-01-04"), by="day")
time<-as.Date(c("2016-01-01", "2016-01-02", "2016-01-03", "2016-01-04"))
days<-gsub("-","",time,fixed=TRUE)
daysNum <- length(time)

dir.create("../predictions")

dem <- readGDAL("../dem_twi/dem_cat.tif")
r = raster(dem)
names(dem) = 'dem'
# dem$twi<- readGDAL('../dem_twi/twi_cat.tif')$band1############################################################
gg <- as(dem, "SpatialPixelsDataFrame")
# gg=SpatialPointsDataFrame(coordinates(dem),data=dem@data) 
gg@proj4string = CRS("+proj=longlat +datum=WGS84")
rm(dem) 
gg2 <- spTransform(gg, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
gg_po <- as(gg2, "SpatialPointsDataFrame")
gg_po$staid <- 1:length(gg_po)

###### STRK - prediction ###################################################################

dir.create("../predictions/strk")

load(file = "stfdf_temp.rda")
load(file='../models/STRK_lm.rda')
load(file='../models/STRK_fit_vgm.rda')

sp.nmax = 10
coef = coefficients(lm)
vario = sumMetric_Vgm
temp <- stfdf[ , time]

registerDoParallel(cores=detectCores()-1)
pred <- foreach(i = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  
  i_1=c(rep(1,1),1:(daysNum -1))
  ip1=c(1:daysNum)
  data <- temp[, i_1[i]:ip1[i], "lm_res", drop=F]
  data <-  data[!is.na(data$lm_res), drop=F]
  
  df <- list()
  tmin <- raster(paste("../min/", days[i], ".tif", sep = ""))
  df$tmin <- raster::extract(tmin, gg) / 10
  df <- as.data.frame(df)
  tmax <- raster(paste("../max/", days[i], ".tif", sep = ""))
  df$tmax <- raster::extract(tmax, gg) / 10
  df$tmin <- ifelse(df$tmin > df$tmax, df$tmax-1, df$tmin)
  
  version <- c('06B', '05B', '04B')
  imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V06A.1day.tif", sep = "")))
  if(inherits(imerg, "try-error")) {
    for (v in version) {
      imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V", v, ".1day.tif", sep = "")))
      if(inherits(imerg, "try-error")) {
        next
      } else {
        break
      }
    }
    if(inherits(imerg, "try-error")) {
      return(NULL)
    }
  }
  df$imerg <- raster::extract(imerg, gg)
  
  tlm  = coef[1] + coef[2]*df$imerg+coef[3]*df$tmax + coef[4]*df$tmin
  
  res <- krigeST(as.formula("lm_res~1"),
                 data=data, 
                 newdata=STF(as(gg,"SpatialPoints"),
                             temp@time[i],
                             temp@endTime[i]),
                 # newdata=STF(as(temp@sp[i,],"SpatialPoints"),
                 #             temp@time[ii],  
                 #             temp@endTime[ii]),     
                 modelList=vario,
                 # nmax = sp.nmax,
                 computeVar=T) # @data[,1]
  res$rk_pred <- tlm + res$var1.pred
  res$rk_pred <- ifelse(res$rk_pred<0, 0, res$rk_pred)
  res
}
stopImplicitCluster()

for(j in 1:daysNum){
  pre = round(pred[[j]]$rk_pred * 100)
  iqr <- round(1.349 * sqrt(pred[[j]]$var1.var) * 100)
  gg@data$pred=pre
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/strk/", days[j], sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  gg@data$pred=iqr
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/strk/", days[j], "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  rm(pre, iqr, p)
}

###### IDW - prediction ###################################################################

dir.create("../predictions/idw")

load(file = "stfdf_temp.rda")

pp = 2.2
n = 13

registerDoParallel(cores=detectCores()-1)
foreach(i = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  
  st <- stfdf[, i]
  st <- st[!is.na(st$prcp), "prcp"]
  
  ##############
  
  gs <- gstat(formula=prcp~1, locations=st, nmax = n, set=list(idp = pp))
  pre <- predict(gs, gg)
  pre <- round(pre$var1.pred * 100)
  
  gg@data$pred=pre
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/idw/", days[i], sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  rm(pre, p)
  
}
stopImplicitCluster()

###### RF - prediction ###################################################################

dir.create("../predictions/rf")

load(file = "../models/RF.rda")

registerDoParallel(cores=detectCores()-1)
foreach(i = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  
  df <- list()
  tmin <- raster(paste("../min/", days[i], ".tif", sep = ""))
  df$tmin <- raster::extract(tmin, gg) / 10
  df <- as.data.frame(df)
  tmax <- raster(paste("../max/", days[i], ".tif", sep = ""))
  df$tmax <- raster::extract(tmax, gg) / 10
  df$tmin <- ifelse(df$tmin > df$tmax, df$tmax-1, df$tmin)
  
  version <- c('06B', '05B', '04B')
  imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V06A.1day.tif", sep = "")))
  if(inherits(imerg, "try-error")) {
    for (v in version) {
      imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V", v, ".1day.tif", sep = "")))
      if(inherits(imerg, "try-error")) {
        next
      } else {
        break
      }
    }
    if(inherits(imerg, "try-error")) {
      return(NULL)
    }
  }
  df$imerg <- raster::extract(imerg, gg)
  
  ##############
  pre <- predict(rf_model, df)
  pre <- round(pre$predictions * 100)
  pred_int <- 0.5# 0.6827 # one sigma
  quantiles = c((1-pred_int)/2, 0.5, 1-(1-pred_int)/2)
  quant_reg <- predict(rf_model, df, type = "quantiles", quantiles = quantiles)
  iqr <- (quant_reg$predictions[,3] - quant_reg$predictions[,1]) / 2
  # head(iqr)
  # hist(iqr)
  # mean(iqr); sqrt(rf_avg$prediction.error)
  iqr <- round(iqr * 100)
  
  gg@data$pred=pre
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/rf/", days[i], sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  gg@data$pred=iqr
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/rf/", days[i], "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  rm(pre, iqr, p)
  
}
stopImplicitCluster()

###### RFSI - prediction ###################################################################

dir.create("../predictions/rfsi")

load('stfdf_temp.rda')
load(file = "../models/RFSI.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

temp <- stfdf[ , time]

n_obs = 7

registerDoParallel(cores=detectCores()-1)
foreach(i = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  
  obs_st <- temp[, i, "prcp"]
  obs_st <-  obs_st[!is.na(obs_st$prcp),]
  
  if (length(obs_st) < 2) {
    return(NULL)
  }
  
  df <- list()
  tmin <- raster(paste("../min/", days[i], ".tif", sep = ""))
  df$tmin <- raster::extract(tmin, gg) / 10
  df <- as.data.frame(df)
  tmax <- raster(paste("../max/", days[i], ".tif", sep = ""))
  df$tmax <- raster::extract(tmax, gg) / 10
  df$tmin <- ifelse(df$tmin > df$tmax, df$tmax-1, df$tmin)
  
  version <- c('06B', '05B', '04B')
  imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V06A.1day.tif", sep = "")))
  if(inherits(imerg, "try-error")) {
    for (v in version) {
      imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V", v, ".1day.tif", sep = "")))
      if(inherits(imerg, "try-error")) {
        next
      } else {
        break
      }
    }
    if(inherits(imerg, "try-error")) {
      return(NULL)
    }
  }
  df$imerg <- raster::extract(imerg, gg)
  
  nearest_obs <- near.obs(
    locations = gg2,
    observations = obs_st,
    zcol = "prcp",
    n.obs = n_obs
  )
  
  df <- cbind(df, nearest_obs)
  ##############
  pre <- predict(rfsi_model, df)
  pre <- round(pre$predictions * 100)
  pred_int <- 0.5# 0.6827 # one sigma
  quantiles = c((1-pred_int)/2, 0.5, 1-(1-pred_int)/2)
  quant_reg <- predict(rfsi_model, df, type = "quantiles", quantiles = quantiles)
  iqr <- (quant_reg$predictions[,3] - quant_reg$predictions[,1]) / 2
  # head(iqr)
  # hist(iqr)
  # mean(iqr); sqrt(rf_avg$prediction.error)
  iqr <- round(iqr * 100)
  
  gg@data$pred=pre
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/rfsi/", days[i], sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  gg@data$pred=iqr
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/rfsi/", days[i], "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  rm(pre, iqr, p)
  
}
stopImplicitCluster()

###### RFsp - prediction ###################################################################

dir.create("../predictions/rfsp")

load('stfdf_temp.rda')
load(file = "../models/RFsp.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

temp <- stfdf[ , time]

### Create BBOX ###
xmin=min(gg_po@coords[, "x"])-1000; xmax=max(gg_po@coords[, "x"])+1000;
ymin=min(gg_po@coords[, "y"])-1000 ; ymax=max(gg_po@coords[, "y"])+1000 # Serbia (proj bbox is 166021.4431, 833978.5569, 0.0000, 9329005.1825)
bbox = extent(xmin, xmax, ymin, ymax)
bbox <- as(bbox, "SpatialPolygons")
bbox@proj4string = CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
fishnet <- raster(bbox)
res(fishnet) <- c(1000, 1000) # 1 x 1 km
fishnet[] <- runif(length(fishnet), -10, 10)
fishnet = as(fishnet, "SpatialPixelsDataFrame")
# plot(fishnet)

st <- as(fishnet, "SpatialPoints")

## Buffer distances only:
grid.distP <- GSIF::buffer.dist(stfdf@sp["staid"], fishnet[1], as.factor(1:nrow(stfdf@sp)))
ov_toma <- over(gg_po, grid.distP)

registerDoParallel(cores=detectCores()-1)
foreach(i = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","meteo")) %dopar% {
  
  df <- list()
  tmin <- raster(paste("../min/", days[i], ".tif", sep = ""))
  df$tmin <- raster::extract(tmin, gg) / 10
  df <- as.data.frame(df)
  tmax <- raster(paste("../max/", days[i], ".tif", sep = ""))
  df$tmax <- raster::extract(tmax, gg) / 10
  df$tmin <- ifelse(df$tmin > df$tmax, df$tmax-1, df$tmin)
  
  version <- c('06B', '05B', '04B')
  imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V06A.1day.tif", sep = "")))
  if(inherits(imerg, "try-error")) {
    for (v in version) {
      imerg <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V", v, ".1day.tif", sep = "")))
      if(inherits(imerg, "try-error")) {
        next
      } else {
        break
      }
    }
    if(inherits(imerg, "try-error")) {
      return(NULL)
    }
  }
  df$imerg <- raster::extract(imerg, gg)
  
  df <- cbind(df, ov_toma)
  
  ##############
  pre <- predict(rfsp_model, df)
  pre <- round(pre$predictions * 100)
  pred_int <- 0.5 # 0.6827 # one sigma
  quantiles = c((1-pred_int)/2, 0.5, 1-(1-pred_int)/2)
  quant_reg <- predict(rfsp_model, df, type = "quantiles", quantiles = quantiles)
  iqr <- (quant_reg$predictions[,3] - quant_reg$predictions[,1]) / 2
  # head(iqr)
  # hist(iqr)
  # mean(iqr); sqrt(rf_avg$prediction.error)
  iqr <- round(iqr * 100)
  
  gg@data$pred=pre
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/rfsp/", days[i], sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  gg@data$pred=iqr
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste("../predictions/rfsp/", days[i], "_iqr", sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  rm(pre, iqr, p)
  
}
stopImplicitCluster()

###### Prediction plots###################################################################

load('stfdf_temp.rda')
stfdf <- stfdf[!is.na(stfdf[, 1]$prcp), ]

summary(stfdf[ , time[1], "prcp"])
summary(stfdf[ , time[2], "prcp"])
summary(stfdf[ , time[3], "prcp"])
summary(stfdf[ , time[4], "prcp"])

met <- c("strk", "idw", "rf", "rfsp", "rfsi")

pred <- raster(paste("../predictions/", met[1], "/",
                     days[1], ".tif", sep = "")) / 100
pred <- as(pred, "SpatialPixelsDataFrame")
names(pred) <- paste(met[1], 1, sep = "")

for (t in 1:length(time)){
  for (m in 1:length(met)){
    
    r <- raster(paste("../predictions/", met[m], "/",
                      days[t], ".tif", sep = "")) / 100
    r <- as(r, "SpatialPixelsDataFrame")
    names(r) <- met[m]
    pred[[paste(met[m], t, sep = "")]] <- r[[met[m]]]
  }
}

### levelplot ###

for (t in 1:length(time)){
  for (m in 1:length(met)){
    assign(paste(met[m], t, sep = ""), raster(paste("../predictions/", met[m], "/",
                                                    days[t], ".tif", sep = "")) / 100)
  }
}

max(pred@data)
max(pred@data$rf4)
max(pred@data$rfsp4)
max(pred@data$rfsi4)
max(pred@data$idw4)
max(pred@data$strk4)
min(pred@data)

pal <- colorRampPalette(brewer.pal(9, "Blues"))
cut <- c(floor(min(pred@data)), 0.5, 1, 1.5, 2, 3, 4, 6, 8, 15, ceiling(max(pred@data)))
# cut2 <- cut
cut2 <- c(floor(min(pred@data)), 2, 4, 6, 8, 15, ceiling(max(pred@data)))

s <- stack(strk1, strk2, strk3, strk4,
           idw1, idw2, idw3, idw4,
           rf1, rf2, rf3, rf4,
           rfsp1, rfsp2, rfsp3, rfsp4,
           rfsi1, rfsi2, rfsi3, rfsi4)

met2 <- c("STRK", "IDW", "RF", "RFsp", "RFSI")

time_max <- c()
for (i in 1:4){
  time_max <- c(time_max, paste("max: ", max(stfdf[ , time[i]]$prcp), " mm\n", sep = ""))
}

# tiff("../plot/predictions.tiff", width = 100, height = 150, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/predictions.jpeg", width = 100, height = 150, units = 'mm', res = 600)
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
          scales=list(color="black", draw=TRUE, alternating=c(3), cex=0.5, abbreviate=T, minlength=2, rot=45),            
          col.regions=pal,
          at=cut,
          names.attr=rep('', nlayers(s)),
          xlab.top=list(as.character(time), cex=0.7, fontface='bold'), xlab=list(as.character(time_max), cex=0.7),
          ylab=list(rev(met2), cex=0.7, fontface='bold'),
          layout=c(4, 5),
          as.table = TRUE) +
  layer(sp.points(stfdf@sp, pch=3, col="black", cex=0.1, alpha=0.5)) +
  layer(sp.polygons(border, lwd=0.1)) +
  layer({
    SpatialPolygonsRescale(layout.north.arrow(type = 2),
                           offset = c(2.5,40.5),
                           scale = 1)
  }, packets = 4) +
  layer({
    grid.text(x= 2.2, y=40.9, "N", gp=gpar(cex=0.8), rot=0,
              default.units='native')
  }, packets = 4)+
  layer({
    xs <- seq(1.45, 3.15, by=0.85) # 
    grid.rect(x=xs[-3], y=40.9,
              width=0.85, height=0.07,
              gp=gpar(fill=c('black', 'transparent')),
              default.units='native')
    grid.text(x= xs - 0.4, y=40.7, c(0, 75, "150km"),# seq(0, 150, by=75), # 0 100 200
              gp=gpar(cex=0.5), rot=0,
              default.units='native')
  }, packets = 20)
dev.off()

###### IQR plots###################################################################

load('stfdf_temp.rda')
stfdf <- stfdf[!is.na(stfdf[, 1]$prcp), ]

met <- c("strk", "rf", "rfsp", "rfsi")

pred <- raster(paste("../predictions/", met[1], "/",
                     days[1], "_iqr.tif", sep = "")) / 100
pred <- as(pred, "SpatialPixelsDataFrame")
names(pred) <- paste(met[1], 1, sep = "")

for (t in 1:length(time)){
  for (m in 1:length(met)){
    
    r <- raster(paste("../predictions/", met[m], "/",
                      days[t], "_iqr.tif", sep = "")) / 100
    r <- as(r, "SpatialPixelsDataFrame")
    names(r) <- met[m]
    pred[[paste(met[m], t, sep = "")]] <- r[[met[m]]]
  }
}

### levelplot ###

for (t in 1:length(time)){
  for (m in 1:length(met)){
    assign(paste(met[m], t, sep = ""), raster(paste("../predictions/", met[m], "/",
                                                    days[t], "_iqr.tif", sep = "")) / 100)
  }
}

max(pred@data)
min(pred@data)

pal <- colorRampPalette(brewer.pal(9, "YlOrRd")[1:9])
cut <- c(0, 0.5, 1, 2, 3, 4, 5, 6, 7, 13, ceiling(max(pred@data)))

s <- stack(strk1, strk2, strk3, strk4,
           rf1, rf2, rf3, rf4,
           rfsp1, rfsp2, rfsp3, rfsp4,
           rfsi1, rfsi2, rfsi3, rfsi4)

met2 <- c("STRK", "RF", "RFsp", "RFSI")

# tiff("../plot/iqr.tiff", width = 100, height = 130, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/iqr.jpeg", width = 100, height = 130, units = 'mm', res = 600)
levelplot(s, 
          margin=FALSE,                       
          colorkey=list(
            space='bottom',                   
            labels=list(at=cut[-2],cex=0.7),
            axis.line=list(col='black'),
            width=0.75
          ),    
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(color="black", draw=TRUE, alternating=c(3), cex=0.5, abbreviate=T, minlength=2, rot=45),         
          col.regions=pal,
          at=cut,
          names.attr=rep('', nlayers(s)),
          xlab.top=list(as.character(time), cex=0.7, fontface='bold'), xlab=NULL,
          ylab=list(rev(met2), cex=0.7, fontface='bold'),
          layout=c(4, 4)) +
  layer(sp.points(stfdf@sp, pch=3, col="black", cex=0.1, alpha=0.5)) +
  layer(sp.polygons(border, lwd=0.1)) +
  layer({
    SpatialPolygonsRescale(layout.north.arrow(type = 2),
                           offset = c(2.5,40.5),
                           scale = 1)
  }, packets = 4) +
  layer({
    grid.text(x= 2.2, y=40.9, "N", gp=gpar(cex=0.8), rot=0,
              default.units='native')
  }, packets = 4)+
  layer({
    xs <- seq(1.45, 3.15, by=0.85)
    grid.rect(x=xs[-3], y=40.9,
              width=0.85, height=0.07,
              gp=gpar(fill=c('black', 'transparent')),
              default.units='native')
    grid.text(x= xs - 0.4, y=40.7, c(0, 75, "150km"),# seq(0, 150, by=75),
              gp=gpar(cex=0.5), rot=0,
              default.units='native')
  }, packets = 16)
dev.off()
