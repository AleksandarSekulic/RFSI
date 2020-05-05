### install latest version of R meteo package ###
# install.packages("meteo", repos="http://R-Forge.R-project.org")
library(meteo) # rforge

library(sp)
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
library(ggpubr)
library(grid)
library(hexbin)
library(ggsn)
library(mlr)
library(tuneRanger)

###### functions  ######
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}
scale_x_longitude <- function(xmin=-180, xmax=180, step=1, limits=NA, ...) {
  xbreaks <- seq(xmin,xmax,step)
  xlabels <- unlist(lapply(xbreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*W")), ifelse(x > 0, parse(text=paste0(x,"^o", "*E")),x))))
  if (is.na(limits)){
    return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), ...))
  } else {
    return(scale_x_continuous("Longitude", breaks = xbreaks, labels = xlabels, expand = c(0, 0), limits=limits, ...))
  }
  
}
scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, limits=NA, ...) {
  ybreaks <- seq(ymin,ymax,step)
  ylabels <- unlist(lapply(ybreaks, function(x) ifelse(x < 0, parse(text=paste0(x,"^o", "*S")), ifelse(x > 0, parse(text=paste0(x,"^o", "*N")),x))))
  if (is.na(limits)){
    return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), ...))
  } else {
    return(scale_y_continuous("Latitude", breaks = ybreaks, labels = ylabels, expand = c(0, 0), limits=limits, ...))
  }
  
}

###### initial data ######

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(42)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

utm33 <- "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

border_wgs84 <- readOGR("borders/osm/union_of_selected_boundaries_AL2-AL2.shp")
border <- spTransform(border_wgs84, CRS(utm33))

dir.create("plot")
dir.create("temp_data")
setwd(paste(wd, "/temp_data/", sep = ""))

v = "TEMP"
year="2008"
# time=seq(as.POSIXct(paste(year,"-01-01", sep="")), as.POSIXct(paste(year,"-12-31", sep="")), by="day")
# days<-gsub("-","",time,fixed=TRUE)
# daysNum <- length(time)
covariates <- c("TEMP", "HRdem","HRdsea","Lat","Lon","HRtwi","INSOL","ctd","MODIS.LST")


###### prepare Croatian dataset (2018) ###################################################################

load("HRtemp2008locsMODIS.rda")
utm33 <- "+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

ob <- HRtemp2008locsMODIS[, c(1:5,7)]
sta <- HRtemp2008locsMODIS[, c(2,8:17)]
sta <- sta[!duplicated(sta), ]

stfdf <- meteo2STFDF ( obs      = ob,
                       stations = sta,
                       crs      = CRS(utm33), # CRS("+proj=longlat +datum=WGS84"),
                       obs.ST_ID.time=c(2,1),
                       stations.ST_ID.lon.lat=c(1,10,11)
)
stfdf <- stfdf[, -367]

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

stfdf = rm.dupl(stfdf, zcol = 1, zero.tol = 0.6) # 0 stations removed

nrowsp <- length(stfdf@sp)

numNA <- apply(matrix(stfdf@data[,"TEMP"],
                      nrow=nrowsp,byrow=F), MARGIN=1,
               FUN=function(x) sum(is.na(x)))

# Remove stations out of covariates
rem <- numNA != daysNum
stfdf <-  stfdf[rem,drop=F] # 1 station

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

theta <- min(stfdf$cday, na.rm=T)
stfdf$ctd <- cos((stfdf$cday-theta)*pi/180)
stfdf$cday <- NULL
stfdf@sp@data$staid <- 1:nrow(stfdf@sp)

# save(stfdf, file="stfdf.rda")

###### stations plot ######

load("stfdf.rda")
load("grids.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=longlat +datum=WGS84"))
st_proj <- stfdf@sp

r <- grids["HRdem"]
# r <- spTransform(r, CRS("+proj=longlat +datum=WGS84"))
r <- raster(r)
# writeRaster(r, "../dem.tiff", "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
r <- projectRaster(r, crs=CRS("+proj=longlat +datum=WGS84"))
r = crop(r, border_wgs84)
r = mask(r, border_wgs84)
r <- as(r, "SpatialPixelsDataFrame")

theme = theme_set(theme_minimal())
sta_dem_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(r), aes(x=x, y=y, fill=HRdem), alpha=0.8) +
  scale_fill_gradientn(colours = terrain.colors(100), name = "DEM [m]") +
  geom_polygon(data = border_wgs84, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  geom_point(data = as.data.frame(stfdf), aes(x = lon, y = lat, color = "red", shape = as.factor("staid")), size = 0.5) +
  theme(plot.title = element_text(hjust = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6),
        text = element_text(size = 6),
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=6, face="bold"),
        legend.text=element_text(size=6),
        legend.position = c(1.15,.4),
        plot.margin=unit(c(5.5,60,5.5,5.5),"points")) +
  labs(x = "x", y = "y") +
  scale_colour_manual(name = "Stations",
                      labels = c("Croatia"),
                      values = c("red"="red")) +
  scale_shape_manual(name = "Stations",
                     labels = c("Croatia"),
                     values = c(17)) +
  scalebar(x.min = 13, x.max = 20,
           y.min = 42, y.max = 47,
           st.size = 3.5, location="bottomright", st.dist=0.04, border.size=0.5,
           dist = 75, dist_unit = "km",
           transform = T, model = "WGS84",
           anchor=c(x=15.5, y=43)) +
  north(x.min = 13, x.max = 20,
        y.min = 42, y.max = 47,
        scale = 0.2, symbol = 3, location="bottomright",
        anchor=c(x=19, y=43.5)) +
  scale_x_longitude(xmin=13, xmax=20, step=1) +
  scale_y_latitude(ymin=42, ymax=47, step=1) +
  coord_fixed()
  # scalebar(x.min = 400000, x.max = 800000,
  #          y.min = 4700000, y.max = 5150000,
  #          st.size = 2.5, location="bottomright", st.dist=0.05, border.size=0.5,
  #          dist = 75, dist_unit = "km",
  #          transform = F, model = "WGS84",
  #          anchor=c(x=550000, y=4750000)) +
  # north(x.min = 400000, x.max = 800000,
  #       y.min = 4700000, y.max = 5150000,
  #       scale = 0.2, symbol = 3, location="bottomright",
  #       anchor=c(x=800000, y=4850000)) +
  # coord_fixed() +
  # scale_x_continuous(breaks = seq(300000, 900000, 100000), labels = seq(300000, 900000, 100000)/100000, expand = c(0, 0)) + 
  # scale_y_continuous(breaks = seq(4600000, 5200000, 100000), labels = seq(4600000, 5200000, 100000)/100000, expand = c(0, 0))
# sta_dem_plot

### Fig7 ###
# tiff("../plot/stations.tiff", width = 100, height = 62, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/stations.jpeg", width = 100, height = 62, units = 'mm', res = 600)
sta_dem_plot
dev.off()

###### create 10 space-time folds ###################################################################

load("stfdf.rda")
hist(stfdf@data$TEMP)

### Fig8 ###
# tiff("../plot/histogram.tiff", width = 80, height = 75, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/histogram.jpeg", width = 80, height = 75, units = 'mm', res = 600)
par(cex = 0.7, mar=c(4,4,0.5,0.5))
res_hist = hist(stfdf$TEMP, breaks=50, plot=F)
cuts = cut(res_hist$breaks, c(-Inf, 0, Inf))
# multiplier <- res_hist$counts / res_hist$density
# hist_density <- density(stfdf$TEMP, na.rm=T)
# hist_density$y <- hist_density$y * multiplier[1]
# lines(hist_density, col="red", lwd=2) 
plot(res_hist, main=NULL, ylim=c(0, 5000), xlab=expression(paste("Temperature [",degree,"C]"))) # , col=c("red","white")[cuts]
# text(x=35, y=4800, labels=paste(res_hist$counts[1], "\n(< 5 mm)", sep=""), col="red")
# text(x=150, y=4600, labels="3rd Quartile: 0.2 mm\nMean: 2.0 mm\nMaximum: 220.9 mm")
dev.off()

source("../stratfolds.R")
test = as.data.frame(stfdf@sp)
test$ID = index(stfdf@sp)
test$weight <-rep(1,dim(test)[1])
strat <- stratfold3d(target.name = c("lon", "lat"), # stavi TEMP
                     other.names = c("lon", "lat"),
                     data = test,
                     num.folds = 10,
                     num.means = 10,
                     seed = 42,
                     cum.prop = 0.9)


# set.seed(42)
# indices <- CreateSpacetimeFolds(test, spacevar = "ST_ID", k=10, seed = 42)
# strat$obs.fold.list
# indices

# save(strat, file="folds.rda")
load(file="folds.rda")

data <- strat$data
by(data$dem, data$fold, summary)
by(data$lon, data$fold, length)

### function for folds plot ###
plotfolds<-function(folds){
  allData<-folds$data
  allData.unique<-ddply(allData,.(ID),here(summarize),longitude=lon[1],latitude=lat[1],fold=fold[1])
  q <- ggplot(allData.unique,aes(x = longitude, y = latitude))
  r <- q + geom_polygon(data = border, aes(x = long, y = lat, group = group), alpha = 0.8, color = "black", fill="white", size = 0.1)
  r <- r + geom_point(pch = 21, alpha=0.6, size=1)# + scale_size_continuous(range=c(1,10))
  r <- r + facet_wrap(~ fold)
  r <- r + aes(fill = fold) + labs(fill = "Fold")
  
  r <- r + theme(plot.title = element_text(hjust = 0.5),
                 axis.text = element_text(size = 4, color="black"),
                 axis.title = element_text(size = 7),
                 text = element_text(size = 7),
                 legend.position = "bottom",
                 legend.key.size= unit(0.2, "cm"),
                 # legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=7, face="bold")) +
    labs(x = "Longitude", y = "Latitude")
  plot(r)
}

# tiff("../plot/folds.tiff", width = 70, height = 70, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/folds.jpeg", width = 70, height = 70, units = 'mm', res = 600)
plotfolds(folds = strat)
dev.off()


###### create IDW ######

load(file = "stfdf.rda")

indices <- CreateSpacetimeFolds(stfdf@sp,spacevar = "ST_ID",
                                k=5, seed = 42)

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

# take all days
# tune p, 1 to 5 step 0.1 for 50 nearest points
# for each p
# for each fold 1-5
# for each day 1-50
n <- 11 # 10 # 50 # found by itterration
p_list <- seq(1,3,0.1)
results <- data.frame(p=p_list, rmse=NA)
for (p in p_list) {
  print(paste("P = ", p, sep = ""))
  
  registerDoParallel(cores=detectCores()-1)
  res <- foreach(val_fold = 1:5, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
    idw_obs <- c()
    idw_pred <- c()
    # for (val_fold in 1:5) {
    # print(paste("Fold ", val_fold, sep = ""))
    val_ids = indices$indexOut[[val_fold]]
    dev_ids = indices$index[[val_fold]]
    dev = stfdf[dev_ids, ]
    val = stfdf[val_ids, ]
    for (day in days){
      # for (day in 1:daysNum){
      train <- dev[, day]
      train <- train[!is.na(train$TEMP), ]
      test <- val[, day]
      test <- test[!is.na(test$TEMP), ]
      gs <- gstat(formula=TEMP~1, locations=train, nmax = n, set=list(idp = p))
      pred <- predict(gs, test)
      idw_obs <- c(idw_obs, test$TEMP)
      idw_pred <- c(idw_pred, pred$var1.pred)
    }
    data.frame(obs=idw_obs, pred=idw_pred)
  }
  res <- do.call("rbind", res)
  
  rmse <- RMSE(res$obs, res$pred)
  results[results$p == p, 2] <- rmse
  print(paste("RMSE = ", rmse, sep = ""))
}
results[which.min(results$rmse), ]
# p     rmse
# 9 1.8 1.784149

# tune number of points, 10 to 70, step 5 for tuned p
p <- results[which.min(results$rmse), 1] # 2.2 # 1.8
n_list <- seq(5,30,1)
results <- data.frame(n=n_list, rmse=NA)
for (n in n_list) {
  print(paste("N = ", n, sep = ""))
  registerDoParallel(cores=detectCores()-1)
  res <- foreach(val_fold = 1:5, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
    print(paste("Fold ", val_fold, sep = ""))
    idw_obs <- c()
    idw_pred <- c()
    val_ids = indices$indexOut[[val_fold]]
    dev_ids = indices$index[[val_fold]]
    dev = stfdf[dev_ids, ]
    val = stfdf[val_ids, ]
    for (day in days){
      # for (day in 1:daysNum){
      train <- dev[, day]
      train <- train[!is.na(train$TEMP), ]
      test <- val[, day]
      test <- test[!is.na(test$TEMP), ]
      gs <- gstat(formula=TEMP~1, locations=train, nmax = n, set=list(idp = p))
      pred <- predict(gs, test)
      idw_obs <- c(idw_obs, test$TEMP)
      idw_pred <- c(idw_pred, pred$var1.pred)
    }
    data.frame(obs=idw_obs, pred=idw_pred)
  }
  res <- do.call("rbind", res)
  
  rmse <- RMSE(res$obs, res$pred)
  results[results$n == n, 2] <- rmse
  print(paste("RMSE = ", rmse, sep = ""))
}
results[which.min(results$rmse), ]
# n     rmse
# 7 11 1.784149

# tune p again for tuned number of points
n <- results[which.min(results$rmse), 1]
p_list <- seq(1,3,0.1)
results <- data.frame(p=p_list, rmse=NA)
for (p in p_list) {
  print(paste("P = ", p, sep = ""))
  idw_obs <- c()
  idw_pred <- c()
  registerDoParallel(cores=detectCores()-1)
  res <- foreach(val_fold = 1:5, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
    idw_obs <- c()
    idw_pred <- c()
    # for (val_fold in 1:5) {
    # print(paste("Fold ", val_fold, sep = ""))
    val_ids = indices$indexOut[[val_fold]]
    dev_ids = indices$index[[val_fold]]
    dev = stfdf[dev_ids, ]
    val = stfdf[val_ids, ]
    for (day in days){
      # for (day in 1:daysNum){
      train <- dev[, day]
      train <- train[!is.na(train$TEMP), ]
      test <- val[, day]
      test <- test[!is.na(test$TEMP), ]
      gs <- gstat(formula=TEMP~1, locations=train, nmax = n, set=list(idp = p))
      pred <- predict(gs, test)
      idw_obs <- c(idw_obs, test$TEMP)
      idw_pred <- c(idw_pred, pred$var1.pred)
    }
    data.frame(obs=idw_obs, pred=idw_pred)
  }
  res <- do.call("rbind", res)
  
  rmse <- RMSE(res$obs, res$pred)
  results[results$p == p, 2] <- rmse
  print(paste("RMSE = ", rmse, sep = ""))
}
results[which.min(results$rmse), ]
# p     rmse
# 9 1.8 1.784149

# p = 1.8
# n = 11

###### IDW - 10-fold CV ######

p = 1.8
n = 11

load(file = "stfdf.rda")
load(file="folds.rda")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

IDW_10f_obs <- c()
IDW_10f_pred <- c()

for(val_fold in 1:10){
  print(paste("Fold ", val_fold, sep = ""))
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  dev = stfdf[stfdf@sp$ind==1, ]
  val = stfdf[stfdf@sp$ind==2, ]
  
  indices <- CreateSpacetimeFolds(dev@sp,spacevar = "ST_ID",
                                  k=5, seed = 42)
  
  # tune 9 folds
  n <- 10 # 10 # 50 # found by itterration
  p_list <- seq(1,3,0.1)
  results <- data.frame(p=p_list, rmse=NA)
  for (p in p_list) {
    print(paste("P = ", p, sep = ""))
    
    registerDoParallel(cores=detectCores()-1)
    res <- foreach(f = 1:5, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
      idw_obs <- c()
      idw_pred <- c()
      # for (f in 1:5) {
      # print(paste("Fold ", f, sep = ""))
      val_ids1 = indices$indexOut[[f]]
      dev_ids1 = indices$index[[f]]
      dev1 = stfdf[dev_ids1, ]
      val1 = stfdf[val_ids1, ]
      for (day in days){
        # for (day in 1:daysNum){
        train <- dev1[, day]
        train <- train[!is.na(train$TEMP), ]
        test <- val1[, day]
        test <- test[!is.na(test$TEMP), ]
        gs <- gstat(formula=TEMP~1, locations=train, nmax = n, set=list(idp = p))
        pred <- predict(gs, test)
        idw_obs <- c(idw_obs, test$TEMP)
        idw_pred <- c(idw_pred, pred$var1.pred)
      }
      data.frame(obs=idw_obs, pred=idw_pred)
    }
    res <- do.call("rbind", res)
    
    rmse <- RMSE(res$obs, res$pred)
    results[results$p == p, 2] <- rmse
    print(paste("RMSE = ", rmse, sep = ""))
  }
  print(results[which.min(results$rmse), ])
  # p     rmse
  # 8 1.7 1.785192
  
  # tune number of points, 10 to 70, step 5 for tuned p
  p <- results[which.min(results$rmse), 1] # 2.2 # 1.8
  n_list <- seq(5,15,1)
  results <- data.frame(n=n_list, rmse=NA)
  for (n in n_list) {
    print(paste("N = ", n, sep = ""))
    registerDoParallel(cores=detectCores()-1)
    res <- foreach(f = 1:5, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
      print(paste("Fold ", f, sep = ""))
      idw_obs <- c()
      idw_pred <- c()
      val_ids1 = indices$indexOut[[f]]
      dev_ids1 = indices$index[[f]]
      dev1 = stfdf[dev_ids1, ]
      val1 = stfdf[val_ids1, ]
      for (day in days){
        # for (day in 1:daysNum){
        train <- dev1[, day]
        train <- train[!is.na(train$TEMP), ]
        test <- val1[, day]
        test <- test[!is.na(test$TEMP), ]
        gs <- gstat(formula=TEMP~1, locations=train, nmax = n, set=list(idp = p))
        pred <- predict(gs, test)
        idw_obs <- c(idw_obs, test$TEMP)
        idw_pred <- c(idw_pred, pred$var1.pred)
      }
      data.frame(obs=idw_obs, pred=idw_pred)
    }
    res <- do.call("rbind", res)
    
    rmse <- RMSE(res$obs, res$pred)
    results[results$n == n, 2] <- rmse
    print(paste("RMSE = ", rmse, sep = ""))
  }
  print(results[which.min(results$rmse), ])
  # n     rmse
  # 2 6 1.973889
  n <- results[which.min(results$rmse), 1]
  
  registerDoParallel(cores=detectCores()-1)
  results <- foreach(d = 1:length(days), .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
    day <- days[d]
    # for (day in days){
    # for (day in 1:daysNum){
    train <- dev[, day]
    train <- train[!is.na(train$TEMP), ]
    test <- val[, day]
    test <- test[!is.na(test$TEMP), ]
    gs <- gstat(formula=TEMP~1, locations=train, nmax = n, set=list(idp = p))
    pred <- predict(gs, test)
    # idw_obs <- c(idw_obs, test$TEMP)
    # idw_pred <- c(idw_pred, pred$var1.pred)
    data.frame(obs=test$TEMP, pred=pred$var1.pred)
  }
  stopImplicitCluster()
  results <- do.call("rbind", results)
  
  IDW_10f_obs <- c(IDW_10f_obs, results$obs)
  IDW_10f_pred <- c(IDW_10f_pred, results$pred)
  
}

# save(IDW_10f_pred, file = "IDW_10f_pred.rda")
# save(IDW_10f_obs, file = "IDW_10f_obs.rda")
load(file = "IDW_10f_pred.rda")
load(file = "IDW_10f_obs.rda")

## RMSE
sqrt(mean((IDW_10f_obs - IDW_10f_pred)^2, na.rm = T))
# 1.753416
## CCC
DescTools::CCC(IDW_10f_obs, IDW_10f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.9742897
## MAE
mean(abs(IDW_10f_obs - IDW_10f_pred), na.rm=TRUE)
# 1.158677
## R2
1 - (t(IDW_10f_obs - IDW_10f_pred) %*% (IDW_10f_obs - IDW_10f_pred)) / (t(IDW_10f_obs - mean(IDW_10f_obs)) %*% (IDW_10f_obs - mean(IDW_10f_obs)))
# 0.9497547
## ME
mean((IDW_10f_obs - IDW_10f_pred), na.rm=TRUE)
# -0.03289852

###### create RF model ###################################################################

load(file = "stfdf.rda")

temp_df = as.data.frame(stfdf)
temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", covariates)]

temp_df = temp_df[complete.cases(temp_df), ]

time=sort(unique(temp_df$time))
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

min.node.size <- 2:20
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
# splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:7

indices <- CreateSpacetimeFolds(temp_df,spacevar = "sp.ID",
                                k=5, seed = 42)

hp <- expand.grid(min.node.size=min.node.size,
                  mtry=mtry, sf=sample.fraction)
hp <- hp[sample(nrow(hp), 50),]
hp <- hp[order(hp$mtry),]
rmse_hp <- rep(NA, nrow(hp))

for (h in 1:nrow(hp)) {
  comb <- hp[h, ]
  print(paste("combination: ", h, sep=""))
  print(comb)
  
  fold_obs <- c()
  fold_pred <- c()
  
  for (f in 1:length(indices$index)) {
    print(paste("fold: ", f, sep=""))
    dev_df1 <- temp_df[indices$index[[f]], covariates]
    val_df1 <- temp_df[indices$indexOut[[f]], covariates]
    
    model <- ranger(TEMP ~ ., data = dev_df1, importance = "none", seed = 42,
                    num.trees = ntree, mtry = comb$mtry,
                    splitrule = "variance",
                    min.node.size = comb$min.node.size,
                    sample.fraction = comb$sf,
                    oob.error = FALSE)
    fold_obs <- c(fold_obs, val_df1$TEMP)
    fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
    
  }
  rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
  
  print(paste("rmse: ", rmse_hp[h], sep=""))
  
}

dev_parameters <- hp[which.min(rmse_hp), ]
print(dev_parameters)
# min.node.size mtry   sf
# 401             3    6 0.85

temp_df <- temp_df[, covariates]

rf_model <- ranger(TEMP ~ ., data = temp_df, importance = "impurity", seed = 42,
                   num.trees = ntree, mtry = dev_parameters$mtry,
                   splitrule = "variance",
                   min.node.size = dev_parameters$min.node.size,
                   sample.fraction = dev_parameters$sf,
                   quantreg = TRUE) ### quantreg???

# save(rf_model, file = "../models/RF.rda")
load(file = "../models/RF.rda")

rf_model$r.squared
# 0.9770036
sqrt(rf_model$prediction.error)
# 1.156616

###### RF - 10-fold CV ###################################################################

load(file = "stfdf.rda")
load(file="folds.rda")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

rf_10f_obs <- c()
rf_10f_pred <- c()

# min.node.size mtry   sf
#            3    6 0.85
min.node.size <- 1:10
sample.fraction <- seq(1, 0.70, -0.05) # 0.632 without / 1 with replacement
# splitrule <- "variance"
ntree <- 250 # 500
mtry <- 3:7

for(val_fold in 1:10){
  cat("fold ")
  cat(val_fold)
  cat("\n")
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  temp_df$HRdem <- as.numeric(as.character(temp_df$HRdem))
  temp_df$HRdsea <- as.numeric(as.character(temp_df$HRdsea))
  temp_df$HRtwi <- as.numeric(as.character(temp_df$HRtwi))
  temp_df$Lat <- as.numeric(as.character(temp_df$Lat))
  temp_df$Lon <- as.numeric(as.character(temp_df$Lon))
  
  temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", "ind", covariates)]
  temp_df = temp_df[complete.cases(temp_df), ]
  
  dev_df <- temp_df[temp_df$ind==1, ]
  val_df <- temp_df[temp_df$ind==2, ]
  
  ### Tune RF ###
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                  k=5, seed = 42)
  
  hp <- expand.grid(min.node.size=min.node.size,
                    mtry=mtry, sf=sample.fraction)
  hp <- hp[sample(nrow(hp), 50),]
  hp <- hp[order(hp$mtry),]
  rmse_hp <- rep(NA, nrow(hp))
  
  for (h in 1:nrow(hp)) {
    comb <- hp[h, ]
    print(paste("combination: ", h, sep=""))
    print(comb)
    
    fold_obs <- c()
    fold_pred <- c()
    
    for (f in 1:length(indices$index)) {
      print(paste("fold: ", f, sep=""))
      dev_df1 <- dev_df[indices$index[[f]], covariates]
      val_df1 <- dev_df[indices$indexOut[[f]], covariates]
      
      model <- ranger(TEMP ~ ., data = dev_df1, importance = "none", seed = 42,
                      num.trees = ntree, mtry = comb$mtry,
                      splitrule = "variance",
                      min.node.size = comb$min.node.size,
                      sample.fraction = comb$sf,
                      oob.error = FALSE)
      fold_obs <- c(fold_obs, val_df1$TEMP)
      fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
      
    }
    rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
    
    print(paste("rmse: ", rmse_hp[h], sep=""))
    
  }
  
  dev_parameters <- hp[which.min(rmse_hp), ]
  print(dev_parameters)
  cat("\n")
  
  dev_df <- dev_df[, covariates]
  
  dev_model <- ranger(TEMP ~ ., data = dev_df, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = dev_parameters$mtry,
                      splitrule = "variance",
                      min.node.size = dev_parameters$min.node.size,
                      sample.fraction = dev_parameters$sf)
  
  ### validation ###
  val_df <- val_df[, covariates]
  rf_pred <- predict(dev_model, val_df)
  
  rf_10f_obs <- c(rf_10f_obs, val_df$TEMP)
  rf_10f_pred <- c(rf_10f_pred, rf_pred$predictions)
}

# save(rf_10f_pred, file="rf_10f_pred.rda")
# save(rf_10f_obs, file="rf_10f_obs.rda")

load(file="rf_10f_pred.rda")
load(file="rf_10f_obs.rda")

rf_10f_pred <- rf_10f_pred[!is.na(rf_10f_obs)]
rf_10f_obs <- rf_10f_obs[!is.na(rf_10f_obs)]

## RMSE
sqrt(mean((rf_10f_obs - rf_10f_pred)^2, na.rm = T))
# 1.574023
## CCC
DescTools::CCC(rf_10f_obs, rf_10f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.9779913
## MAE
mean(abs(rf_10f_obs - rf_10f_pred), na.rm=TRUE)
# 1.100526
## R2
1 - (t(rf_10f_obs - rf_10f_pred) %*% (rf_10f_obs - rf_10f_pred)) / (t(rf_10f_obs - mean(rf_10f_obs)) %*% (rf_10f_obs - mean(rf_10f_obs)))
# 0.9574095
## ME
mean((rf_10f_obs - rf_10f_pred), na.rm=TRUE)
# -0.006456036

###### create RFsp model ###################################################################

load(file = "stfdf.rda")

### Create BBOX ###
xmin=min(stfdf@sp@coords[, "lon"])-1000; xmax=max(stfdf@sp@coords[, "lon"])+1000;
ymin=min(stfdf@sp@coords[, "lat"])-1000 ; ymax=max(stfdf@sp@coords[, "lat"])+1000 # Serbia (proj bbox is 166021.4431, 833978.5569, 0.0000, 9329005.1825)
bbox = extent(xmin, xmax, ymin, ymax)
bbox <- as(bbox, "SpatialPolygons")
bbox@proj4string = CRS(utm33)
fishnet <- raster(bbox)
res(fishnet) <- c(1000, 1000) # 1 x 1 km
fishnet[] <- runif(length(fishnet), -10, 10)
fishnet = as(fishnet, "SpatialPixelsDataFrame")
# plot(fishnet)
st <- as(fishnet, "SpatialPoints")

## Geographic distances only:
grid.distP <- GSIF::buffer.dist(stfdf@sp["staid"], fishnet[1], as.factor(1:nrow(stfdf@sp)))

temp_df = as.data.frame(stfdf)
temp_df$staid <- as.integer(as.character(temp_df$staid))
ov_toma <- do.call(cbind, list(stfdf@sp@data, over(stfdf@sp, grid.distP)))
temp_df <- plyr::join(temp_df, ov_toma, by="staid")

temp_df <- temp_df[, c(covariates, paste("layer.", sort(unique(as.integer(as.character(temp_df$staid)))), sep=""))]

temp_df = temp_df[complete.cases(temp_df), ]

temp_df$HRdem <- as.numeric(as.character(temp_df$HRdem))
temp_df$HRdsea <- as.numeric(as.character(temp_df$HRdsea))
temp_df$HRtwi <- as.numeric(as.character(temp_df$HRtwi))
temp_df$Lat <- as.numeric(as.character(temp_df$Lat))
temp_df$Lon <- as.numeric(as.character(temp_df$Lon))

time=sort(unique(temp_df$time))
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

ntree <- 250 # 500

# tune RFsp as Hengl et al. 2018 did
rt <- makeRegrTask(data = temp_df[sample(nrow(temp_df), 10000),], target = "TEMP")
estimateTimeTuneRanger(rt, num.threads = 88)
# Time consuming >> do on number cruncher
set.seed(42)
t.rfsp <- tuneRanger(rt, num.trees = ntree, build.final.model = FALSE)

dev_parameters <- t.rfsp$recommended.pars
print(dev_parameters)
# mtry min.node.size sample.fraction     mse exec.time
# 1  154             2       0.7670135 3.31092   13.0664

rfsp_model <- ranger(TEMP ~ ., data = temp_df, importance = "impurity", seed = 42,
                     num.trees = ntree, mtry = dev_parameters$mtry,
                     splitrule = "variance",
                     min.node.size = dev_parameters$min.node.size,
                     sample.fraction = dev_parameters$sample.fraction,
                     quantreg = TRUE) ### quantreg???

# save(rfsp_model, file = "../models/RFsp.rda")
load(file = "../models/RFsp.rda")

rfsp_model$r.squared
# 0.9767293
sqrt(rfsp_model$prediction.error)
# 1.163492

###### RFsp- 10-fold CV ###################################################################

load(file = "stfdf.rda")
load(file="folds.rda")

### Create BBOX ###
xmin=min(stfdf@sp@coords[, "lon"])-1000; xmax=max(stfdf@sp@coords[, "lon"])+1000;
ymin=min(stfdf@sp@coords[, "lat"])-1000 ; ymax=max(stfdf@sp@coords[, "lat"])+1000 # Serbia (proj bbox is 166021.4431, 833978.5569, 0.0000, 9329005.1825)
bbox = extent(xmin, xmax, ymin, ymax)
bbox <- as(bbox, "SpatialPolygons")
bbox@proj4string = CRS(utm33)
fishnet <- raster(bbox)
res(fishnet) <- c(1000, 1000) # 1 x 1 km
fishnet[] <- runif(length(fishnet), -10, 10)
fishnet = as(fishnet, "SpatialPixelsDataFrame")
plot(fishnet)

st <- as(fishnet, "SpatialPoints")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

ntree <- 250


rfsp_10f_obs <- c()
rfsp_10f_pred <- c()

## Buffer distances only:
grid.distP <- GSIF::buffer.dist(stfdf@sp["staid"], fishnet[1], as.factor(1:nrow(stfdf@sp)))

for(val_fold in 1:10){
  cat(val_fold)
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  temp_df$HRdem <- as.numeric(as.character(temp_df$HRdem))
  temp_df$HRdsea <- as.numeric(as.character(temp_df$HRdsea))
  temp_df$HRtwi <- as.numeric(as.character(temp_df$HRtwi))
  temp_df$Lat <- as.numeric(as.character(temp_df$Lat))
  temp_df$Lon <- as.numeric(as.character(temp_df$Lon))
  
  temp_df$staid <- as.integer(as.character(temp_df$staid))
  
  ov_toma <- do.call(cbind, list(stfdf@sp@data, over(stfdf@sp, grid.distP)))
  temp_df <- plyr::join(temp_df, ov_toma, by="staid")
  
  # temp_df <- temp_df[, c(covariates, paste("layer.", sort(unique(as.integer(as.character(temp_df$staid)))), sep=""))]
  
  temp_df = temp_df[complete.cases(temp_df), ]
  
  dev_df <- temp_df[temp_df$ind==1, c(covariates, paste("layer.", sort(unique(as.integer(as.character(temp_df$staid[temp_df$ind==1])))), sep=""))]
  val_df <- temp_df[temp_df$ind==2, c(covariates, paste("layer.", sort(unique(as.integer(as.character(temp_df$staid[temp_df$ind==1])))), sep=""))]
  
  # # tune RFsp as Hengl et al. 2018 did
  rt <- makeRegrTask(data = dev_df[sample(nrow(dev_df), 10000),], target = "TEMP")
  # estimateTimeTuneRanger(rt, num.threads = 88)
  # # Time consuming >> do on number cruncher
  set.seed(42)
  t.rfsp <- tuneRanger(rt, num.trees = ntree, build.final.model = FALSE)
  dev_parameters <- t.rfsp$recommended.pars
  print(dev_parameters)
  # dev_parameters <- data.frame(mtry=154, min.node.size=2, sample.fraction=0.7670135)
  # mtry min.node.size sample.fraction     mse exec.time
  # 1  154             2       0.7670135 3.31092   13.0664
  
  dev_model <- ranger(TEMP ~ ., data = dev_df, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = dev_parameters$mtry,
                      splitrule = "variance",
                      min.node.size = dev_parameters$min.node.size, # always 2
                      sample.fraction = dev_parameters$sample.fraction)
  
  ### validation ###
  rfsp_pred <- predict(dev_model, val_df)
  
  rfsp_10f_obs <- c(rfsp_10f_obs, val_df$TEMP)
  rfsp_10f_pred <- c(rfsp_10f_pred, rfsp_pred$predictions)
  
}

# save(rfsp_10f_pred, file="rfsp_10f_pred.rda")
# save(rfsp_10f_obs, file="rfsp_10f_obs.rda")
load(file="rfsp_10f_pred.rda")
load(file="rfsp_10f_obs.rda")

rfsp_10f_pred <- rfsp_10f_pred[!is.na(rfsp_10f_obs)]
rfsp_10f_obs <- rfsp_10f_obs[!is.na(rfsp_10f_obs)]

## RMSE
sqrt(mean((rfsp_10f_obs - rfsp_10f_pred)^2, na.rm = T))
# 1.626198
## CCC
DescTools::CCC(rfsp_10f_obs, rfsp_10f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.9763768
## MAE
mean(abs(rfsp_10f_obs - rfsp_10f_pred), na.rm=TRUE)
# 1.132449
## R2
1 - (t(rfsp_10f_obs - rfsp_10f_pred) %*% (rfsp_10f_obs - rfsp_10f_pred)) / (t(rfsp_10f_obs - mean(rfsp_10f_obs)) %*% (rfsp_10f_obs - mean(rfsp_10f_obs)))
# 0.9545391
## ME
mean((rfsp_10f_obs - rfsp_10f_pred), na.rm=TRUE)
# -0.003914756

###### create RFSI model ###################################################################

load("stfdf.rda")

temp_df = as.data.frame(stfdf)
temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", covariates)]

temp_df = temp_df[complete.cases(temp_df), ]

time=sort(unique(temp_df$time))
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

nos <- 6:14#c(6,8,10,12,14)
min.node.size <- 2:20
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
# splitrule <- "variance"
ntree <- 250 # 500

n_obs = max(nos)
mtry <- 3:(9+2*n_obs)

indices <- CreateSpacetimeFolds(temp_df,spacevar = "sp.ID",
                                k=5, seed = 42)

hp <- expand.grid(min.node.size=min.node.size,
                  mtry=mtry, no=nos, sf=sample.fraction)
hp <- hp[hp$mtry < (9+2*hp$no-1), ]
hp <- hp[sample(nrow(hp), 100),]
hp <- hp[order(hp$no),]
rmse_hp <- rep(NA, nrow(hp))

for (h in 1:nrow(hp)) {
  comb <- hp[h, ]
  print(paste("combination: ", h, sep=""))
  print(comb)
  
  dev_df <- temp_df
  
  fold_obs <- c()
  fold_pred <- c()
  
  for (f in 1:length(indices$index)) {
    print(paste("fold: ", f, sep=""))
    cols <- c(covariates, paste("dist", 1:(comb$no), sep=""), paste("obs", 1:(comb$no), sep=""))
    dev_df1 <- dev_df[indices$index[[f]], ]
    
    cpus <- detectCores()-1
    registerDoParallel(cores=cpus)
    nearest_obs <- foreach (t = time) %dopar% {
      
      dev_day_df <- dev_df1[dev_df1$time==t, c("lon", "lat", "TEMP")]
      day_df <- dev_df[dev_df$time==t, c("lon", "lat", "TEMP")]
      
      if (nrow(day_df)==0) {
        return(NULL)
      }
      return(near.obs(
        locations = day_df,
        observations = dev_day_df,
        zcol = "TEMP",
        n.obs = comb$no
      ))
      
    }
    stopImplicitCluster()
    
    nearest_obs <- do.call("rbind", nearest_obs)
    dev_df <- cbind(dev_df, nearest_obs)
    dev_df = dev_df[complete.cases(dev_df), ]
    
    
    dev_df1 <- dev_df[indices$index[[f]], cols]
    dev_df1 <- dev_df1[complete.cases(dev_df1), ]
    val_df1 <- dev_df[indices$indexOut[[f]], cols]
    val_df1 <- val_df1[complete.cases(val_df1), ]
    
    model <- ranger(TEMP ~ ., data = dev_df1, importance = "none", seed = 42,
                    num.trees = ntree, mtry = comb$mtry,
                    splitrule = "variance",
                    min.node.size = comb$min.node.size,
                    sample.fraction = comb$sf,
                    oob.error = FALSE)
    fold_obs <- c(fold_obs, val_df1$TEMP)
    fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
    
  }
  rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
  
  print(paste("rmse: ", rmse_hp[h], sep=""))
  
}

dev_parameters <- hp[which.min(rmse_hp), ]
print(dev_parameters)
# min.node.size mtry no   sf
# 5344             6    4  7 0.95

cpus <- detectCores()-1
registerDoParallel(cores=cpus)
nearest_obs <- foreach (t = time) %dopar% {
  
  day_df <- temp_df[temp_df$time==t, c("lon", "lat", "TEMP")]
  if (nrow(day_df)==0) {
    return(NULL)
  }
  return(near.obs(
    locations = day_df,
    observations = day_df,
    zcol = "TEMP",
    n.obs = dev_parameters$no
  ))
  
}
stopImplicitCluster()

nearest_obs <- do.call("rbind", nearest_obs)
temp_df <- cbind(temp_df, nearest_obs)

cols <- c(covariates, paste("dist", 1:(dev_parameters$no), sep=""), paste("obs", 1:(dev_parameters$no), sep=""))
temp_df <- temp_df[, cols]

rfsi_model <- ranger(TEMP ~ ., data = temp_df, importance = "impurity", seed = 42,
                     num.trees = ntree, mtry = dev_parameters$mtry,
                     splitrule = "variance",
                     min.node.size = dev_parameters$min.node.size,
                     sample.fraction = dev_parameters$sf,
                     quantreg = TRUE) ### quantreg???

# save(rfsi_model, file = "../models/RFSI.rda")
load(file = "../models/RFSI.rda")

rfsi_model$r.squared
# 0.736633
sqrt(rfsi_model$prediction.error)
# 3.530753

###### RFSI - 10-fold CV ###################################################################

load(file = "stfdf.rda")
load(file="folds.rda")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

# min.node.size mtry no  sf
#            15    5 10 0.9
nos <- 6:14# 14 # c(6,8,10,12,14)
min.node.size <- 10:20
sample.fraction <- seq(1, 0.75, -0.05) # 0.632 without / 1 with replacement
# splitrule <- "variance"
ntree <- 250 # 500

n_obs = max(nos)
mtry <- 3:(9+2*n_obs)

rfsi_10f_obs <- c()
rfsi_10f_pred <- c()

for(val_fold in 1:10){
  cat("fold ")
  cat(val_fold)
  cat("\n")
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  
  temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", "ind", covariates)]
  temp_df = temp_df[complete.cases(temp_df), ]
  
  time=sort(unique(temp_df$time))
  daysNum = length(time)
  days=gsub("-","",time,fixed=TRUE)
  
  dev_df <- temp_df[temp_df$ind==1, ]
  val_df <- temp_df[temp_df$ind==2, ]
  
  ### Tune RFSI ###
  indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                  k=5, seed = 42)
  
  hp <- expand.grid(min.node.size=min.node.size,
                    mtry=mtry, no=nos, sf=sample.fraction)
  hp <- hp[hp$mtry < (9+2*hp$no-1), ]
  hp <- hp[sample(nrow(hp), 50),]
  hp <- hp[order(hp$no),]
  rmse_hp <- rep(NA, nrow(hp))
  
  for (h in 1:nrow(hp)) {
    comb <- hp[h, ]
    print(paste("combination: ", h, sep=""))
    print(comb)
    
    fold_obs <- c()
    fold_pred <- c()
    
    for (f in 1:length(indices$index)) {
      print(paste("fold: ", f, sep=""))
      cols <- c(covariates, paste("dist", 1:(comb$no), sep=""), paste("obs", 1:(comb$no), sep=""))
      
      dev_df1 <- dev_df[indices$index[[f]], ]
      dev_df1 <- dev_df1[complete.cases(dev_df1), ]
      
      cpus <- detectCores()-1
      registerDoParallel(cores=cpus)
      nearest_obs <- foreach (t = time) %dopar% {
        
        dev_day_df <- dev_df1[dev_df1$time==t, c("lon", "lat", "TEMP", "ind")]
        day_df <- dev_df[dev_df$time==t, c("lon", "lat", "TEMP", "ind")]
        
        if (nrow(day_df)==0) {
          return(NULL)
        }
        return(near.obs(
          locations = day_df,
          observations = dev_day_df,#[day_df$ind==1, ],
          zcol = "TEMP",
          n.obs = comb$no
        ))
        
      }
      stopImplicitCluster()
      
      nearest_obs <- do.call("rbind", nearest_obs)
      dev_df <- cbind(dev_df, nearest_obs)
      dev_df = dev_df[complete.cases(dev_df), ]
      
      
      dev_df1 <- dev_df[indices$index[[f]], cols]
      dev_df1 <- dev_df1[complete.cases(dev_df1), ]
      val_df1 <- dev_df[indices$indexOut[[f]], cols]
      val_df1 <- val_df1[complete.cases(val_df1), ]
      
      model <- ranger(TEMP ~ ., data = dev_df1, importance = "none", seed = 42,
                      num.trees = ntree, mtry = comb$mtry,
                      splitrule = "variance",
                      min.node.size = comb$min.node.size,
                      sample.fraction = comb$sf,
                      oob.error = FALSE)
      fold_obs <- c(fold_obs, val_df1$TEMP)
      fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
      
    }
    rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
    
    print(paste("rmse: ", rmse_hp[h], sep=""))
    
  }
  
  dev_parameters <- hp[which.min(rmse_hp), ]
  print(dev_parameters)
  cat("\n")
  
  dev_df <- temp_df[temp_df$ind==1, ]
  
  cpus <- detectCores()-1
  registerDoParallel(cores=cpus)
  nearest_obs <- foreach (t = time) %dopar% {
    
    dev_day_df <- dev_df[dev_df$time==t, c("lon", "lat", "TEMP", "ind")]
    day_df <- temp_df[temp_df$time==t, c("lon", "lat", "TEMP", "ind")]
    
    if (nrow(day_df)==0) {
      return(NULL)
    }
    return(near.obs(
      locations = day_df,
      observations = dev_day_df,#[day_df$ind==1, ],
      zcol = "TEMP",
      n.obs = dev_parameters$no
    ))
    
  }
  stopImplicitCluster()
  
  nearest_obs <- do.call("rbind", nearest_obs)
  temp_df <- cbind(temp_df, nearest_obs)
  temp_df = temp_df[complete.cases(temp_df), ]
  
  cols <- c(covariates, paste("dist", 1:(dev_parameters$no), sep=""), paste("obs", 1:(dev_parameters$no), sep=""))
  
  dev_df <- temp_df[temp_df$ind==1, cols]
  val_df <- temp_df[temp_df$ind==2, cols]
  
  dev_model <- ranger(TEMP ~ ., data = dev_df, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = dev_parameters$mtry,
                      splitrule = "variance",
                      min.node.size = dev_parameters$min.node.size,
                      sample.fraction = dev_parameters$sf)
  
  ### validation ###
  val_df <- val_df[, cols]
  rfsi_pred <- predict(dev_model, val_df)
  
  rfsi_10f_obs <- c(rfsi_10f_obs, val_df$TEMP)
  rfsi_10f_pred <- c(rfsi_10f_pred, rfsi_pred$predictions)
}

# save(rfsi_10f_pred, file="rfsi_10f_pred.rda")
# save(rfsi_10f_obs, file="rfsi_10f_obs.rda")

load(file="rfsi_10f_pred.rda")
load(file="rfsi_10f_obs.rda")

rfsi_10f_pred <- rfsi_10f_pred[!is.na(rfsi_10f_obs)]
rfsi_10f_obs <- rfsi_10f_obs[!is.na(rfsi_10f_obs)]

## RMSE
sqrt(mean((rfsi_10f_obs - rfsi_10f_pred)^2, na.rm = T))
# 1.398816
## CCC
DescTools::CCC(rfsi_10f_obs, rfsi_10f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.9827464
## MAE
mean(abs(rfsi_10f_obs - rfsi_10f_pred), na.rm=TRUE)
# 0.9550509
## R2
1 - (t(rfsi_10f_obs - rfsi_10f_pred) %*% (rfsi_10f_obs - rfsi_10f_pred)) / (t(rfsi_10f_obs - mean(rfsi_10f_obs)) %*% (rfsi_10f_obs - mean(rfsi_10f_obs)))
# 0.9663635
## ME
mean((rfsi_10f_obs - rfsi_10f_pred), na.rm=TRUE)
# -0.02586521
###### create RFSI model without covariates ###################################################################

load("stfdf.rda")

temp_df = as.data.frame(stfdf)
temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", covariates)]

temp_df = temp_df[complete.cases(temp_df), ]

time=sort(unique(temp_df$time))
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

nos <- 6:14#c(6,8,10,12,14)
min.node.size <- 2:20
sample.fraction <- seq(1, 0.632, -0.05) # 0.632 without / 1 with replacement
# splitrule <- "variance"
ntree <- 250 # 500

n_obs = max(nos)
mtry <- 3:(2*n_obs)

indices <- CreateSpacetimeFolds(temp_df,spacevar = "sp.ID",
                                k=5, seed = 42)

hp <- expand.grid(min.node.size=min.node.size,
                  mtry=mtry, no=nos, sf=sample.fraction)
hp <- hp[hp$mtry < (2*hp$no-1), ]
hp <- hp[sample(nrow(hp), 100),]
hp <- hp[order(hp$no),]
rmse_hp <- rep(NA, nrow(hp))

for (h in 1:nrow(hp)) {
  comb <- hp[h, ]
  print(paste("combination: ", h, sep=""))
  print(comb)
  
  dev_df <- temp_df
  
  fold_obs <- c()
  fold_pred <- c()
  
  for (f in 1:length(indices$index)) {
    print(paste("fold: ", f, sep=""))
    cols <- c("TEMP", paste("dist", 1:(comb$no), sep=""), paste("obs", 1:(comb$no), sep=""))
    dev_df1 <- dev_df[indices$index[[f]], ]
    
    cpus <- detectCores()-1
    registerDoParallel(cores=cpus)
    nearest_obs <- foreach (t = time, .packages = c("meteo")) %dopar% {
      
      dev_day_df <- dev_df1[dev_df1$time==t, c("lon", "lat", "TEMP")]
      day_df <- dev_df[dev_df$time==t, c("lon", "lat", "TEMP")]
      
      if (nrow(day_df)==0) {
        return(NULL)
      }
      return(near.obs(
        locations = day_df,
        observations = dev_day_df,
        zcol = "TEMP",
        n.obs = comb$no
      ))
      
    }
    stopImplicitCluster()
    
    nearest_obs <- do.call("rbind", nearest_obs)
    dev_df <- cbind(dev_df, nearest_obs)
    dev_df = dev_df[complete.cases(dev_df), ]
    
    
    dev_df1 <- dev_df[indices$index[[f]], cols]
    dev_df1 <- dev_df1[complete.cases(dev_df1), ]
    val_df1 <- dev_df[indices$indexOut[[f]], cols]
    val_df1 <- val_df1[complete.cases(val_df1), ]
    
    model <- ranger(TEMP ~ ., data = dev_df1, importance = "none", seed = 42,
                    num.trees = ntree, mtry = comb$mtry,
                    splitrule = "variance",
                    min.node.size = comb$min.node.size,
                    sample.fraction = comb$sf,
                    oob.error = FALSE)
    fold_obs <- c(fold_obs, val_df1$TEMP)
    fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
    
  }
  rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
  
  print(paste("rmse: ", rmse_hp[h], sep=""))
  
}

dev_parameters <- hp[which.min(rmse_hp), ]
print(dev_parameters)
# min.node.size mtry no  sf
# 19820             4    6 10 0.8

cpus <- detectCores()-1
registerDoParallel(cores=cpus)
nearest_obs <- foreach (t = time, .packages = c("meteo")) %dopar% {
  
  day_df <- temp_df[temp_df$time==t, c("lon", "lat", "TEMP")]
  if (nrow(day_df)==0) {
    return(NULL)
  }
  return(near.obs(
    locations = day_df,
    observations = day_df,
    zcol = "TEMP",
    n.obs = dev_parameters$no
  ))
  
}
stopImplicitCluster()

nearest_obs <- do.call("rbind", nearest_obs)
temp_df <- cbind(temp_df, nearest_obs)

cols <- c("TEMP", paste("dist", 1:(dev_parameters$no), sep=""), paste("obs", 1:(dev_parameters$no), sep=""))
temp_df <- temp_df[, cols]

rfsi2_model <- ranger(TEMP ~ ., data = temp_df, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = dev_parameters$mtry,
                      splitrule = "variance",
                      min.node.size = dev_parameters$min.node.size,
                      sample.fraction = dev_parameters$sf,
                      quantreg = TRUE) ### quantreg???

# save(rfsi2_model, file = "../models/RFSI2.rda")
load(file = "../models/RFSI2.rda")

rfsi2_model$r.squared
# 0.9808053
sqrt(rfsi2_model$prediction.error)
# 1.056696

###### RFSI without covariates - 10-fold CV ###################################################################

load(file = "stfdf.rda")
load(file="folds.rda")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

# min.node.size mtry no  sf
#            15    5 10 0.9
nos <- 6:14# 14 # c(6,8,10,12,14)
min.node.size <- 10:20
sample.fraction <- seq(1, 0.75, -0.05) # 0.632 without / 1 with replacement
# splitrule <- "variance"
ntree <- 250 # 500

n_obs = max(nos)
mtry <- 3:(2*n_obs-1)

rfsi2_10f_obs <- c()
rfsi2_10f_pred <- c()

for(val_fold in 1:10){
  cat("fold ")
  cat(val_fold)
  cat("\n")
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  
  temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", "ind", "TEMP")]
  temp_df = temp_df[complete.cases(temp_df), ]
  
  time=sort(unique(temp_df$time))
  daysNum = length(time)
  days=gsub("-","",time,fixed=TRUE)
  
  dev_df <- temp_df[temp_df$ind==1, ]
  val_df <- temp_df[temp_df$ind==2, ]
  
  ### Tune RFSI ###
  indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                  k=5, seed = 42)
  
  hp <- expand.grid(min.node.size=min.node.size,
                    mtry=mtry, no=nos, sf=sample.fraction)
  hp <- hp[hp$mtry < (2*hp$no-1), ]
  hp <- hp[sample(nrow(hp), 50),]
  hp <- hp[order(hp$no),]
  rmse_hp <- rep(NA, nrow(hp))
  
  for (h in 1:nrow(hp)) {
    comb <- hp[h, ]
    print(paste("combination: ", h, sep=""))
    print(comb)
    
    fold_obs <- c()
    fold_pred <- c()
    
    for (f in 1:length(indices$index)) {
      print(paste("fold: ", f, sep=""))
      cols <- c("TEMP", paste("dist", 1:(comb$no), sep=""), paste("obs", 1:(comb$no), sep=""))
      
      dev_df1 <- dev_df[indices$index[[f]], ]
      dev_df1 <- dev_df1[complete.cases(dev_df1), ]
      
      cpus <- detectCores()-1
      registerDoParallel(cores=cpus)
      nearest_obs <- foreach (t = time) %dopar% {
        
        dev_day_df <- dev_df1[dev_df1$time==t, c("lon", "lat", "TEMP", "ind")]
        day_df <- dev_df[dev_df$time==t, c("lon", "lat", "TEMP", "ind")]
        
        if (nrow(day_df)==0) {
          return(NULL)
        }
        return(near.obs(
          locations = day_df,
          observations = dev_day_df,#[day_df$ind==1, ],
          zcol = "TEMP",
          n.obs = comb$no
        ))
        
      }
      stopImplicitCluster()
      
      nearest_obs <- do.call("rbind", nearest_obs)
      dev_df <- cbind(dev_df, nearest_obs)
      dev_df = dev_df[complete.cases(dev_df), ]
      
      
      dev_df1 <- dev_df[indices$index[[f]], cols]
      dev_df1 <- dev_df1[complete.cases(dev_df1), ]
      val_df1 <- dev_df[indices$indexOut[[f]], cols]
      val_df1 <- val_df1[complete.cases(val_df1), ]
      
      model <- ranger(TEMP ~ ., data = dev_df1, importance = "none", seed = 42,
                      num.trees = ntree, mtry = comb$mtry,
                      splitrule = "variance",
                      min.node.size = comb$min.node.size,
                      sample.fraction = comb$sf,
                      oob.error = FALSE)
      fold_obs <- c(fold_obs, val_df1$TEMP)
      fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
      
    }
    rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm=TRUE))
    
    print(paste("rmse: ", rmse_hp[h], sep=""))
    
  }
  
  dev_parameters <- hp[which.min(rmse_hp), ]
  print(dev_parameters)
  cat("\n")
  
  dev_df <- temp_df[temp_df$ind==1, ]
  
  cpus <- detectCores()-1
  registerDoParallel(cores=cpus)
  nearest_obs <- foreach (t = time) %dopar% {
    
    dev_day_df <- dev_df[dev_df$time==t, c("lon", "lat", "TEMP", "ind")]
    day_df <- temp_df[temp_df$time==t, c("lon", "lat", "TEMP", "ind")]
    
    if (nrow(day_df)==0) {
      return(NULL)
    }
    return(near.obs(
      locations = day_df,
      observations = dev_day_df,#[day_df$ind==1, ],
      zcol = "TEMP",
      n.obs = dev_parameters$no
    ))
    
  }
  stopImplicitCluster()
  
  nearest_obs <- do.call("rbind", nearest_obs)
  temp_df <- cbind(temp_df, nearest_obs)
  temp_df = temp_df[complete.cases(temp_df), ]
  
  cols <- c("TEMP", paste("dist", 1:(dev_parameters$no), sep=""), paste("obs", 1:(dev_parameters$no), sep=""))
  
  dev_df <- temp_df[temp_df$ind==1, cols]
  val_df <- temp_df[temp_df$ind==2, cols]
  
  dev_model <- ranger(TEMP ~ ., data = dev_df, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = dev_parameters$mtry,
                      splitrule = "variance",
                      min.node.size = dev_parameters$min.node.size,
                      sample.fraction = dev_parameters$sf)
  
  ### validation ###
  val_df <- val_df[, cols]
  rfsi_pred <- predict(dev_model, val_df)
  
  rfsi2_10f_obs <- c(rfsi2_10f_obs, val_df$TEMP)
  rfsi2_10f_pred <- c(rfsi2_10f_pred, rfsi_pred$predictions)
}

# save(rfsi2_10f_pred, file="rfsi2_10f_pred.rda")
# save(rfsi2_10f_obs, file="rfsi2_10f_obs.rda")

load(file="rfsi2_10f_pred.rda")
load(file="rfsi2_10f_obs.rda")

rfsi2_10f_pred <- rfsi2_10f_pred[!is.na(rfsi2_10f_obs)]
rfsi2_10f_obs <- rfsi2_10f_obs[!is.na(rfsi2_10f_obs)]

## RMSE
sqrt(mean((rfsi2_10f_obs - rfsi2_10f_pred)^2, na.rm = T))
# 1.764938
## CCC
DescTools::CCC(rfsi2_10f_obs, rfsi2_10f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.9738406
## MAE
mean(abs(rfsi2_10f_obs - rfsi2_10f_pred), na.rm=TRUE)
# 1.174279
## R2
1 - (t(rfsi2_10f_obs - rfsi2_10f_pred) %*% (rfsi2_10f_obs - rfsi2_10f_pred)) / (t(rfsi2_10f_obs - mean(rfsi2_10f_obs)) %*% (rfsi2_10f_obs - mean(rfsi2_10f_obs)))
# 0.9490922
## ME
mean((rfsi2_10f_obs - rfsi2_10f_pred), na.rm=TRUE)
# 0.01721196
###### importance plot ###################################################################

load(file = "../models/RFSI.rda")

xlP.g <- as.list(round(rfsi_model$variable.importance))
df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:15,]/100000
print(df)
# HRdem+HRdsea+Lat+Lon+HRtwi+INSOL+cos((cday-theta)*pi/180)+MODIS.LST
names(df)[names(df)=="HRdem"] <- "DEM"
names(df)[names(df)=="HRdsea"] <- "Distance-to-\n-coastline"
names(df)[names(df)=="HRtwi"] <- "TWI"
names(df)[names(df)=="INSOL"] <- "Insolation"
names(df)[names(df)=="Lat"] <- "Latitude"
names(df)[names(df)=="Lon"] <- "Longitude"
names(df)[names(df)=="ctd"] <- "Seasonal\nfluctuation" # cos((cday-theta)*pi/180 - model seasonal fluctuation of daily temperature
pr <- 15:1
df <- as.data.frame(cbind(cbind(df, names(df)), pr))
df$V2 <- as.character(df$V2)
names(df) <- c("importance", "covariate", "order")

df$importance <- as.numeric(as.character(df$importance))
# df$importance <- (df$importance-min(df$importance))/(max(df$importance)-min(df$importance))
df$importance <- df$importance / max(df$importance)
summary(df$importance)

# Reorder the data
df <- df %>%
  arrange(importance) %>%
  mutate(covariate=factor(covariate,covariate))

theme = theme_set(theme_minimal())
# tiff("../plot/importance_rfsi.tiff", width = 35, height = 70, units = 'mm', res = 600, compression = "lzw")
# jpeg("../plot/importance_rfsi.jpeg", width = 35, height = 70, units = 'mm', res = 600)
# Plot
cov_names <- c("DEM","Distance-to-\n-coastline","TWI", "Insolation", "Seasonal\nfluctuation", "MODIS.LST", "Longitude", "Latitude")
rfsi_importance <- 
  ggplot(df, aes(x=covariate, y=importance)) +
  geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% cov_names, "red", "black"), size=ifelse(df$covariate %in% cov_names, .5, .25) ) +
  geom_point( color=ifelse(df$covariate %in% cov_names, "red", "black"), size=ifelse(df$covariate %in% cov_names, .5, .25) ) +
  # theme_ipsum() +
  coord_flip() +
  theme(
    legend.position="none"
  ) +
  xlab("Covariate") +
  ylab("Importance") +
  ggtitle("RFSI") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text = element_text(size = 7))
# dev.off()

###
load(file = "../models/RF.rda")

xlP.g <- as.list(round(rf_model$variable.importance))
df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:8,]/100000
print(df)
# HRdem+HRdsea+Lat+Lon+HRtwi+INSOL+cos((cday-theta)*pi/180)+MODIS.LST
names(df)[names(df)=="HRdem"] <- "DEM"
names(df)[names(df)=="HRdsea"] <- "Distance-to-\n-coastline"
names(df)[names(df)=="HRtwi"] <- "TWI"
names(df)[names(df)=="INSOL"] <- "Insolation"
names(df)[names(df)=="Lat"] <- "Latitude"
names(df)[names(df)=="Lon"] <- "Longitude"
names(df)[names(df)=="ctd"] <- "Seasonal\nfluctuation" # cos((cday-theta)*pi/180 - model seasonal fluctuation of daily temperature
pr <- 8:1
df <- as.data.frame(cbind(cbind(df, names(df)), pr))
df$V2 <- as.character(df$V2)
names(df) <- c("importance", "covariate", "order")

df$importance <- as.numeric(as.character(df$importance))

df$importance <- as.numeric(as.character(df$importance))
# df$importance <- (df$importance-min(df$importance))/(max(df$importance)-min(df$importance))
df$importance <- df$importance / max(df$importance)
summary(df$importance)

# Reorder the data
df <- df %>%
  arrange(importance) %>%
  mutate(covariate=factor(covariate,covariate))

# tiff("../plot/importance_rfsi.tiff", width = 35, height = 70, units = 'mm', res = 600, compression = "lzw")
# jpeg("../plot/importance_rfsi.jpeg", width = 35, height = 70, units = 'mm', res = 600)
# Plot
rf_importance <- 
  ggplot(df, aes(x=covariate, y=importance)) +
  geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% cov_names, "red", "black"), size=ifelse(df$covariate %in% cov_names, .5, .25) ) +
  geom_point( color=ifelse(df$covariate %in% cov_names, "red", "black"), size=ifelse(df$covariate %in% cov_names, .5, .25) ) +
  # theme_ipsum() +
  coord_flip() +
  theme(
    legend.position="none"
  ) +
  xlab("Covariate") +
  ylab("Importance") +
  ggtitle("RF") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text = element_text(size = 7))
# dev.off()

###
load(file = "../models/RFsp.rda")

xlP.g <- as.list(round(rfsp_model$variable.importance))
df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:15,]/100000
print(df)
# HRdem+HRdsea+Lat+Lon+HRtwi+INSOL+cos((cday-theta)*pi/180)+MODIS.LST
names(df)[names(df)=="HRdem"] <- "DEM"
names(df)[names(df)=="HRdsea"] <- "Distance-to-\n-coastline"
names(df)[names(df)=="HRtwi"] <- "TWI"
names(df)[names(df)=="INSOL"] <- "Insolation"
names(df)[names(df)=="ctd"] <- "Seasonal\nfluctuation" # cos((cday-theta)*pi/180 - model seasonal fluctuation of daily temperature
names(df)[names(df)=="Lat"] <- "Latitude"
names(df)[names(df)=="Lon"] <- "Longitude"
pr <- 15:1
df <- as.data.frame(cbind(cbind(df, names(df)), pr))
df$V2 <- as.character(df$V2)
names(df) <- c("importance", "covariate", "order")

df$importance <- as.numeric(as.character(df$importance))
# df$importance <- (df$importance-min(df$importance))/(max(df$importance)-min(df$importance))
df$importance <- df$importance / max(df$importance)
summary(df$importance)

# Reorder the data
df <- df %>%
  arrange(importance) %>%
  mutate(covariate=factor(covariate,covariate))

# tiff("../plot/importance_rfsi.tiff", width = 35, height = 70, units = 'mm', res = 600, compression = "lzw")
# jpeg("../plot/importance_rfsi.jpeg", width = 35, height = 70, units = 'mm', res = 600)
# Plot
rfsp_importance <- 
  ggplot(df, aes(x=covariate, y=importance)) +
  geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% cov_names, "red", "black"), size=ifelse(df$covariate %in% cov_names, .5, .25) ) +
  geom_point( color=ifelse(df$covariate %in% cov_names, "red", "black"), size=ifelse(df$covariate %in% cov_names, .5, .25) ) +
  # theme_ipsum() +
  coord_flip() +
  theme(
    legend.position="none"
  ) +
  xlab("Covariate") +
  ylab("Importance") +
  ggtitle("RFsp") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text = element_text(size = 7))
# dev.off()

### Fig12 ###
# tiff("../plot/importance.tiff", width = 110, height = 70, units = 'mm', res = 1200, compression = "lzw")
jpeg("../plot/importance.jpeg", width = 110, height = 70, units = 'mm', res = 1200)
ggarrange(rf_importance, rfsp_importance, rfsi_importance, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

###### CV error histogram, scatter plot ###################################################################

# load(file = "STRK_10f_pred.rda")
# load(file = "STRK_10f_obs.rda")
# summary(STRK_10f_pred)
# STRK_10f_pred <- ifelse(STRK_10f_pred<0, 0, STRK_10f_pred)

load(file = "IDW_10f_obs.rda")
load(file = "IDW_10f_pred.rda")

load(file = "rfsi_10f_obs.rda")
load(file = "rfsi_10f_pred.rda")

load(file = "rfsp_10f_obs.rda")
load(file = "rfsp_10f_pred.rda")

load(file = "rf_10f_obs.rda")
load(file = "rf_10f_pred.rda")

# # tiff("../plot/scaterplot.tiff", width = 174, height = 65, units = 'mm', res = 600, compression = "lzw")
# jpeg("../plot/scaterplot.jpeg", width = 174, height = 65, units = 'mm', res = 600)
# par(mfrow=c(1,4), cex = 0.7, mar=c(4.5,4.5,1,0.5))
# plot(STRK_10f_obs ~ STRK_10f_pred, main="STRK", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "Observations [mm]", cex = 0.5)
# abline(a=0, b=1, col="red")
# plot(rf_10f_obs ~ rf_10f_pred, main="RF", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "", cex = 0.5)
# abline(a=0, b=1, col="red")
# plot(rfsi_10f_obs ~ rfsi_10f_pred, main="RFSI", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "", cex = 0.5)
# abline(a=0, b=1, col="red")
# plot(rfsp_10f_obs ~ rfsp_10f_pred, main="RFsp", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "", cex = 0.5)
# abline(a=0, b=1, col="red")
# dev.off()

# STRK_10f_pred <- STRK_10f_pred[!is.na(STRK_10f_obs)]
# STRK_10f_obs <- STRK_10f_obs[!is.na(STRK_10f_obs)]
# STRK_res <- STRK_10f_obs-STRK_10f_pred
# # STRK_res <- STRK_res[abs(STRK_res)>50]
IDW_10f_res <- IDW_10f_obs-IDW_10f_pred
rf_10f_res <- rf_10f_obs-rf_10f_pred
rfsi_10f_res <- rfsi_10f_obs-rfsi_10f_pred
rfsp_10f_res <- rfsp_10f_obs-rfsp_10f_pred

summary(STRK_res)
summary(IDW_10f_res)
summary(rf_10f_res)
summary(rfsi_10f_res)
summary(rfsp_10f_res)

STRK_res[which.max(STRK_10f_obs)]
rf_10f_res[which.max(rf_10f_obs)]
rfsi_10f_res[which.max(rfsi_10f_obs)]
rfsp_10f_res[which.max(rfsp_10f_obs)]

STRK_10f_obs[which.min(STRK_res)]
STRK_10f_pred[which.min(STRK_res)]
rf_10f_pred[which.min(STRK_res)]
rfsi_10f_pred[which.min(STRK_res)]
rfsp_10f_pred[which.min(STRK_res)]

summary(STRK_res[which(STRK_10f_obs>100)])
summary(rf_10f_res[which(rf_10f_obs>100)])
summary(rfsi_10f_res[which(rfsi_10f_obs>100)])
summary(rfsp_10f_res[which(rfsp_10f_obs>100)])

summary(STRK_res[which(STRK_10f_obs<10)])
summary(rf_10f_res[which(rf_10f_obs<10)])
summary(rfsi_10f_res[which(rfsi_10f_obs<10)])
summary(rfsp_10f_res[which(rfsp_10f_obs<10)])

summary(STRK_res[which(STRK_10f_obs>10 & STRK_10f_obs<100)])
summary(rf_10f_res[which(rf_10f_obs>10 & rf_10f_obs<100)])
summary(rfsi_10f_res[which(rfsi_10f_obs>10 & rfsi_10f_obs<100)])
summary(rfsp_10f_res[which(rfsp_10f_obs>10 & rfsp_10f_obs<100)])

# scaterplots
my_colors=colorRampPalette(rev(bpy.colors()))(50)
# data <- as.data.frame(cbind(STRK_10f_obs, STRK_10f_pred))
# p1 <- ggplot(data, aes(x=STRK_10f_pred, y=STRK_10f_obs) ) +
#   geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 10, 50, 100, 200, 400, Inf)), labels = F, right = T, include.lowest = T)))) +
#   scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('10', '50', '100', '200', '400', '600+'))+
#   theme(plot.title = element_text(hjust = 0.5, face="bold"),
#         axis.text = element_text(size = 8),
#         axis.title = element_text(size = 8),
#         text = element_text(size = 8),
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.key.height = unit(0.3, "cm"),
#         legend.key.width = unit(2, "cm"),
#         legend.margin = unit(0, "cm"),
#         legend.title = element_text(size=8, face="bold"),
#         legend.text=element_text(size=8)) +
#   labs(x = "Residuals [mm]", y = "Observations [mm]", title = "STRK") + coord_fixed()+
#   xlim(0, 35) +
#   geom_abline(slope=1, intercept=0, size = 0.1)
# p1

data <- as.data.frame(cbind(IDW_10f_obs, IDW_10f_pred))
p2 <- ggplot(data, aes(x=IDW_10f_pred, y=IDW_10f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 10, 50, 100, 200, 400, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('10', '50', '100', '200', '400', '600+'))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=8)) +
  labs(x = "Residuals [mm]", y = "", title = "IDW") + coord_fixed()+
  xlim(0, 35) +
  geom_abline(slope=1, intercept=0, size = 0.1)

data <- as.data.frame(cbind(rf_10f_obs, rf_10f_pred))
p3 <- ggplot(data, aes(x=rf_10f_pred, y=rf_10f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 10, 50, 100, 200, 400, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('10', '50', '100', '200', '400', '600+'))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=8)) +
  labs(x = "Residuals [mm]", y = "", title = "RF") + coord_fixed()+
  xlim(0, 35) +
  geom_abline(slope=1, intercept=0, size = 0.1)

data <- as.data.frame(cbind(rfsi_10f_obs, rfsi_10f_pred))
p4 <- ggplot(data, aes(x=rfsi_10f_pred, y=rfsi_10f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 10, 50, 100, 200, 400, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('10', '50', '100', '200', '400', '600+'))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=8)) +
  labs(x = "Residuals [mm]", y = "", title = "RFSI") + coord_fixed()+
  xlim(0, 35) +
  geom_abline(slope=1, intercept=0, size = 0.1)

data <- as.data.frame(cbind(rfsp_10f_obs, rfsp_10f_pred))
p5 <- ggplot(data, aes(x=rfsp_10f_pred, y=rfsp_10f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 10, 50, 100, 200, 400, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('10', '50', '100', '200', '400', '600+'))+
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        text = element_text(size = 8),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=8)) +
  labs(x = "Residuals [mm]", y = "", title = "RFsp") + coord_fixed()+
  xlim(0, 35) +
  geom_abline(slope=1, intercept=0, size = 0.1)

# # histograms
# res_hist1 = ggplot(as.data.frame(STRK_res), aes(x=STRK_res)) +
#   geom_histogram(binwidth=1, fill="white", color="black", alpha=1,size=0.2) +
#   scale_x_continuous(limits = c(-20,20)) +
#   scale_y_continuous(limits = c(0,70000)) +
#   labs(x = "Residuals [mm]", y = "Frequency") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 9),
#         text = element_text(size = 7),
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.key.height = unit(0.3, "cm"),
#         legend.key.width = unit(2, "cm"),
#         legend.margin = unit(0, "cm"),
#         legend.title = element_text(size=7))
# res_hist1
# # rug(STRK_res)
# # res_hist2 = histogram(rf_10f_res, main="", xlim=c(-20,20), ylim=c(0,80000), xlab="Residuals [mm]", ylab="", breaks=300, type = "count", col="white",
# #                       scales=(list(cex=0.5)))
# 
# res_hist2 = ggplot(as.data.frame(rf_10f_res), aes(x=rf_10f_res)) +
#   geom_histogram(binwidth=1, fill="white", color="black", alpha=1,size=0.2) +
#   scale_x_continuous(limits = c(-20,20)) +
#   scale_y_continuous(limits = c(0,70000)) +
#   labs(x = "Residuals [mm]", y = "") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 9),
#         text = element_text(size = 7),
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.key.height = unit(0.3, "cm"),
#         legend.key.width = unit(2, "cm"),
#         legend.margin = unit(0, "cm"),
#         legend.title = element_text(size=7))
# 
# res_hist3 = ggplot(as.data.frame(rfsi_10f_res), aes(x=rfsi_10f_res)) +
#   geom_histogram(binwidth=1, fill="white", color="black", alpha=1,size=0.2) +
#   scale_x_continuous(limits = c(-20,20)) +
#   scale_y_continuous(limits = c(0,70000)) +
#   labs(x = "Residuals [mm]", y = "") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 9),
#         text = element_text(size = 7),
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.key.height = unit(0.3, "cm"),
#         legend.key.width = unit(2, "cm"),
#         legend.margin = unit(0, "cm"),
#         legend.title = element_text(size=7))
# 
# res_hist4 = ggplot(as.data.frame(rfsp_10f_res), aes(x=rfsp_10f_res)) +
#   geom_histogram(binwidth=1, fill="white", color="black", alpha=1,size=0.2) +
#   scale_x_continuous(limits = c(-20,20)) +
#   scale_y_continuous(limits = c(0,70000)) +
#   labs(x = "Residuals [mm]", y = "") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text = element_text(size = 7),
#         axis.title = element_text(size = 9),
#         text = element_text(size = 7),
#         legend.position = "bottom",
#         legend.direction = "horizontal",
#         legend.key.height = unit(0.3, "cm"),
#         legend.key.width = unit(2, "cm"),
#         legend.margin = unit(0, "cm"),
#         legend.title = element_text(size=7))

# tiff("../plot/scatter_cv_hist_res.tiff", width = 174, height = 130, units = 'mm', res = 600)
jpeg("../plot/scatter_cv_hist_res.jpeg", width = 174, height = 70, units = 'mm', res = 600)
# ggarrange(
# ggarrange(p1, p2, p3, p4, p5, ncol=5, nrow=1, common.legend = TRUE, legend="bottom") #,
ggarrange(p2, p3, p4, p5, ncol=4, nrow=1, common.legend = TRUE, legend="bottom") #,
# ggarrange(res_hist1, res_hist2, res_hist3, res_hist4, ncol=4, nrow=1),
# ncol=1, nrow=2)
dev.off()

par(mfrow=c(1,1))












###### Bubble plots ######

load(file = "stfdf.rda")
load(file="folds.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=longlat +datum=WGS84"))
st_proj <- stfdf@sp

load(file = "IDW_10f_obs.rda")
load(file = "IDW_10f_pred.rda")

load(file = "rfsi_10f_obs.rda")
load(file = "rfsi_10f_pred.rda")

load(file = "rfsp_10f_obs.rda")
load(file = "rfsp_10f_pred.rda")

load(file = "rf_10f_obs.rda")
load(file = "rf_10f_pred.rda")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

IDW_10f_staid <- c()
IDW_10f_days <- c()


for(val_fold in 1:10){
  print(paste("Fold ", val_fold, sep = ""))
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  dev = stfdf[stfdf@sp$ind==1, ]
  val = stfdf[stfdf@sp$ind==2, ] 
  val_df <- as.data.frame(val)
  val_df <- val_df[!is.na(val_df$TEMP), ]
  val_df$staid <- as.integer(as.character(val_df$staid))
  
  IDW_10f_staid <- c(IDW_10f_staid, val_df$staid)
  IDW_10f_days <- c(IDW_10f_days, as.character(as.Date(val_df$time)))
  
  
}

rfsi_10f_staid <- c()
rfsi_10f_days <- c()


for(val_fold in 1:10){
  print(paste("Fold ", val_fold, sep = ""))
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  dev = stfdf[stfdf@sp$ind==1, ]
  val = stfdf[stfdf@sp$ind==2, ] 
  val_df <- as.data.frame(val)
  val_df <- val_df[complete.cases(val_df), ]
  val_df$staid <- as.integer(as.character(val_df$staid))
  
  rfsi_10f_staid <- c(rfsi_10f_staid, val_df$staid)
  rfsi_10f_days <- c(rfsi_10f_days, as.character(as.Date(val_df$time)))
  
  
}

# IDW #
# calculate avg rmse per station #
staid <- unique(stfdf@sp$staid)
st = matrix(nrow = length(staid), ncol = 4)
for (i in staid){
  coord <- coordinates(stfdf@sp)[stfdf@sp$staid == i, ]
  rmse = RMSE(IDW_10f_obs[IDW_10f_staid == i], IDW_10f_pred[IDW_10f_staid == i])
  st[i,] = c(i, coord[1], coord[2], rmse)
}
colnames(st) = c("staid", "x", "y", "rmse")
st <- as.data.frame(st)
# coordinates(st) <- c("x", "y")
# st@proj4string <- CRS(utm33)
summary(st$rmse)

theme = theme_set(theme_minimal())
idw_plot <- ggplot(st, aes(x = x, y = y, size = rmse)) +
  geom_polygon(data = border_wgs84, aes(x = long, y = lat, group = group), alpha = 0.8, color = "black", fill="white", size = 0.1) +
  geom_point(shape = 21, colour = "mediumvioletred", fill = "springgreen", alpha = 0.8) +
  # geom_point(data = st[st$rmse<=0, ], aes(x = lon, y = lat, color = "red", size = rmse), stroke=0.1, alpha=0.5) +
  # geom_point(data = st[st$rmse>0, ], aes(x = lon, y = lat, color = "blue", size = rmse), stroke=0.1, alpha=0.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),#element_text(size = 6),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=6, face="bold")) +
  # labs(x = NULL, y = NULL) +
  labs(size = "RMSE") +
  # xlab("Number of nearest locations") +
  # ylab("RMSE") +
  ggtitle("IDW") +
  coord_fixed() +
  # scale_colour_manual(name = "",
  #                     labels = c("positive", "negative"),
  #                     values = c("blue"="blue", "red"="red")) +
  # scale_size_identity(guide="legend", breaks = c(1, 2, 4, 6, 8, 10)) +
  scale_size(range = c(.1, 2.5), breaks = c(1, 2, 4, 6, 8, 10), name="RMSE") +
  scale_x_longitude(xmin=13, xmax=20, step=1) +
  scale_y_latitude(ymin=42, ymax=47, step=1)
  # scale_x_continuous(breaks = seq(300000, 900000, 100000), labels = seq(300000, 900000, 100000)/100000, expand = c(0, 0)) + 
  # scale_y_continuous(breaks = seq(4600000, 5200000, 100000), labels = seq(4600000, 5200000, 100000)/100000, expand = c(0, 0))
  
# RF #
# calculate avg rmse per station #
staid <- unique(stfdf@sp$staid)
st = matrix(nrow = length(staid), ncol = 4)
for (i in staid){
  coord <- coordinates(stfdf@sp)[stfdf@sp$staid == i, ]
  rmse = RMSE(rf_10f_obs[rfsi_10f_staid == i], rf_10f_pred[rfsi_10f_staid == i])
  st[i,] = c(i, coord[1], coord[2], rmse)
}
colnames(st) = c("staid", "x", "y", "rmse")
st <- as.data.frame(st)
# coordinates(st) <- c("x", "y")
# st@proj4string <- CRS(utm33)
summary(st$rmse)

rf_plot <- ggplot(st, aes(x = x, y = y, size = rmse)) +
  geom_polygon(data = border_wgs84, aes(x = long, y = lat, group = group), alpha = 0.8, color = "black", fill="white", size = 0.1) +
  geom_point(shape = 21, colour = "mediumvioletred", fill = "springgreen", alpha = 0.8) +
  # geom_point(data = st[st$rmse<=0, ], aes(x = lon, y = lat, color = "red", size = rmse), stroke=0.1, alpha=0.5) +
  # geom_point(data = st[st$rmse>0, ], aes(x = lon, y = lat, color = "blue", size = rmse), stroke=0.1, alpha=0.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),#element_text(size = 6),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=6, face="bold")) +
  labs(x = NULL, y = NULL) + labs(size = "RMSE") +
  ggtitle("RF") +
  north(x.min = 13, x.max = 20,
        y.min = 42, y.max = 47,
        scale = 0.2, symbol = 3, location="bottomright",
        anchor=c(x=19, y=43.5)) +
  coord_fixed() +
  # north(x.min = 400000, x.max = 800000,
  #       y.min = 4700000, y.max = 5150000,
  #       scale = 0.2, symbol = 3, location="bottomright",
  #       anchor=c(x=800000, y=4850000)) +
  # scale_colour_manual(name = "",
  #                     labels = c("positive", "negative"),
  #                     values = c("blue"="blue", "red"="red")) +
  # scale_size_identity(guide="legend", breaks = c(1, 2, 4, 6, 8, 10)) +
  scale_size(range = c(.1, 2.5), breaks = c(1, 2, 4, 6, 8, 10), name="RMSE") +
  scale_x_longitude(xmin=13, xmax=20, step=1) +
  scale_y_latitude(ymin=42, ymax=47, step=1)
  # scale_x_continuous(breaks = seq(300000, 900000, 100000), labels = seq(300000, 900000, 100000)/100000, expand = c(0, 0)) + 
  # scale_y_continuous(breaks = seq(4600000, 5200000, 100000), labels = seq(4600000, 5200000, 100000)/100000, expand = c(0, 0))

# RFSI #
# calculate avg rmse per station #
staid <- unique(stfdf@sp$staid)
st = matrix(nrow = length(staid), ncol = 4)
for (i in staid){
  coord <- coordinates(stfdf@sp)[stfdf@sp$staid == i, ]
  rmse = RMSE(rfsi_10f_obs[rfsi_10f_staid == i], rfsi_10f_pred[rfsi_10f_staid == i])
  st[i,] = c(i, coord[1], coord[2], rmse)
}
colnames(st) = c("staid", "x", "y", "rmse")
st <- as.data.frame(st)
# coordinates(st) <- c("x", "y")
# st@proj4string <- CRS(utm33)
summary(st$rmse)

rfsi_plot <- ggplot(st, aes(x = x, y = y, size = rmse)) +
  geom_polygon(data = border_wgs84, aes(x = long, y = lat, group = group), alpha = 0.8, color = "black", fill="white", size = 0.1) +
  geom_point(shape = 21, colour = "mediumvioletred", fill = "springgreen", alpha = 0.8) +
  # geom_point(data = st[st$rmse<=0, ], aes(x = lon, y = lat, color = "red", size = rmse), stroke=0.1, alpha=0.5) +
  # geom_point(data = st[st$rmse>0, ], aes(x = lon, y = lat, color = "blue", size = rmse), stroke=0.1, alpha=0.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),#element_text(size = 6),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=6, face="bold")) +
  labs(x = NULL, y = NULL) + labs(size = "RMSE") +
  ggtitle("RFSI") +
  coord_fixed() +
  scalebar(x.min = 13, x.max = 20,
           y.min = 42, y.max = 47,
           st.size = 2.5, location="bottomright", st.dist=0.04, border.size=0.5,
           dist = 75, dist_unit = "km",
           transform = T, model = "WGS84",
           anchor=c(x=15.5, y=43)) +
  # scalebar(x.min = 400000, x.max = 800000,
  #          y.min = 4700000, y.max = 5150000,
  #          st.size = 2.5, location="bottomright", st.dist=0.05, border.size=0.5,
  #          dist = 75, dist_unit = "km",
  #          transform = F, model = "WGS84",
  #          anchor=c(x=550000, y=4750000)) +
  # scale_colour_manual(name = "",
  #                     labels = c("positive", "negative"),
  #                     values = c("blue"="blue", "red"="red")) +
  # scale_size_identity(guide="legend", breaks = c(1, 2, 4, 6, 8, 10)) +
  scale_size(range = c(.1, 2.5), breaks = c(1, 2, 4, 6, 8, 10), name="RMSE") +
  scale_x_longitude(xmin=13, xmax=20, step=1) +
  scale_y_latitude(ymin=42, ymax=47, step=1)
  # scale_x_continuous(breaks = seq(300000, 900000, 100000), labels = seq(300000, 900000, 100000)/100000, expand = c(0, 0)) + 
  # scale_y_continuous(breaks = seq(4600000, 5200000, 100000), labels = seq(4600000, 5200000, 100000)/100000, expand = c(0, 0))

# RFsp #
# calculate avg rmse per station #
staid <- unique(stfdf@sp$staid)
st = matrix(nrow = length(staid), ncol = 4)
for (i in staid){
  coord <- coordinates(stfdf@sp)[stfdf@sp$staid == i, ]
  rmse = RMSE(rfsp_10f_obs[rfsi_10f_staid == i], rfsp_10f_pred[rfsi_10f_staid == i])
  st[i,] = c(i, coord[1], coord[2], rmse)
}
colnames(st) = c("staid", "x", "y", "rmse")
st <- as.data.frame(st)
# coordinates(st) <- c("x", "y")
# st@proj4string <- CRS(utm33)
summary(st$rmse)

rfsp_plot <- ggplot(st, aes(x = x, y = y, size = rmse)) +
  geom_polygon(data = border_wgs84, aes(x = long, y = lat, group = group), alpha = 0.8, color = "black", fill="white", size = 0.1) +
  geom_point(shape = 21, colour = "mediumvioletred", fill = "springgreen", alpha = 0.8) +
  # geom_point(data = st[st$rmse<=0, ], aes(x = lon, y = lat, color = "red", size = rmse), stroke=0.1, alpha=0.5) +
  # geom_point(data = st[st$rmse>0, ], aes(x = lon, y = lat, color = "blue", size = rmse), stroke=0.1, alpha=0.5) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),#element_text(size = 6),
        text = element_text(size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        legend.margin = unit(0, "cm"),
        legend.title = element_text(size=6, face="bold")) +
  labs(x = NULL, y = NULL) + labs(size = "RMSE") +
  ggtitle("RFsp") +
  coord_fixed() +
  # scale_colour_manual(name = "",
  #                     labels = c("positive", "negative"),
  #                     values = c("blue"="blue", "red"="red")) +
  # scale_size_identity(guide="legend", breaks = c(1, 2, 4, 6, 8, 10)) +
  scale_size(range = c(.1, 2.5), breaks = c(1, 2, 4, 6, 8, 10), name="RMSE") +
  scale_x_longitude(xmin=13, xmax=20, step=1) +
  scale_y_latitude(ymin=42, ymax=47, step=1)
  # scale_x_continuous(breaks = seq(300000, 900000, 100000), labels = seq(300000, 900000, 100000)/100000, expand = c(0, 0)) + 
  # scale_y_continuous(breaks = seq(4600000, 5200000, 100000), labels = seq(4600000, 5200000, 100000)/100000, expand = c(0, 0))

# tiff("../plot/bubbles.tiff", width = 130, height = 110, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/bubbles.jpeg", width = 130, height = 100, units = 'mm', res = 600)
ggarrange(idw_plot, rf_plot, rfsp_plot, rfsi_plot, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()





