### install latest version of R meteo package ###
install.packages("meteo", repos="http://R-Forge.R-project.org")
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

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

border <- readOGR("catalonia_border/union_of_selected_boundaries_AL4-AL4.shp")

dir.create("plot")
dir.create("temp_data")
setwd(paste(wd, "/temp_data/", sep = ""))

v = "prcp"
years = 2016:2018
time<-seq(as.Date(paste(as.numeric(min(years)),"-01-01", sep="")), as.Date(paste(max(years),"-12-31", sep="")), by="day")
days<-gsub("-","",time,fixed=TRUE)
daysNum <- length(time)

el = 'PRCP'
names <- c('staid','date','prcp')

###### prepare Catalunya dataset (3 years) ###################################################################

ob <- c()

for (year in years){  
  time<-seq(as.Date(paste(as.numeric(year),"-01-01", sep="")), as.Date(paste(year,"-12-31", sep="")), by="day")
  days<-gsub("-","",time,fixed=TRUE)
  daysNum <- length(time)
  days_rev <- paste(substr(days, 7, 8), substr(days, 5, 6), substr(days, 1, 4), sep = "")
  
  ################# GHCND download ######################
  url <- paste("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/by_year/", year, ".csv.gz", sep = "")
  destfile <- paste("ghcnd_", year, ".csv.gz", sep = "")
  download.file(url, destfile)
  gunzip(destfile)
  
  ghcnd <- read.csv(paste('ghcnd_', year,".csv", sep = ""), header=F)
  names(ghcnd) <- c('staid', 'date', 'elem', 'temp','mflag', 'qflag', 'sflag', 'obstime')
  
  ################# GHCND ######################
  
  ghcnd <- ghcnd[ghcnd$elem == el, ]
  ghcnd <- ghcnd[ghcnd$qflag == "", ]
  ghcnd$temp  <- ghcnd$temp /10
  ghcnd= ghcnd[, c('staid','date','temp')]
  names(ghcnd) = names
  ghcnd$date  <- as.numeric(ghcnd$date)
  ghcnd$date <- as.Date(paste(substr(ghcnd$date, 1, 4), substr(ghcnd$date, 5, 6), substr(ghcnd$date, 7, 8), sep="-"))
  ghcnd$staid = as.character(ghcnd$staid)
  
  ob <- rbind(ob, ghcnd)
  
}

# save(ob, file = "ob.rda")
load(file = "ob.rda")

### GHCND stations ###
url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
destfile <- paste("ghcnd_stations.txt", sep = "")
download.file(url, destfile)

sta_ghcnd <- read.fwf("ghcnd_stations.txt", width = c(11,9,10,7), header=F) # 11 ILI 12 !!!!!!!!!!!!!!!!!!!!!!!!!!
names(sta_ghcnd) <- c('staid', 'lat', 'lon', 'h')
sta_ghcnd$staid <- as.character(sta_ghcnd$staid)
sta_ghcnd$h <- ifelse(sta_ghcnd$h==-999.9, NA, sta_ghcnd$h)
sta_ghcnd = sta_ghcnd[!is.na(sta_ghcnd$lon),]

sta <- sta_ghcnd
coordinates(sta_ghcnd) <- ~ lon + lat
sta_ghcnd@proj4string = border@proj4string
sta <- sta[!is.na(over(sta_ghcnd, border)$id), ]
rm(sta_ghcnd)

# plot(sta)
# plot(border, add=T)

# sta <- sta[sta$staid %in% stfdf@sp$staid, ]
ob <- ob[ob$staid %in% sta$staid, ]

stfdf <- meteo2STFDF ( obs      = ob,
                       stations = sta,
                       crs      = CRS("+proj=longlat +datum=WGS84"),
                       obs.staid.time=c(1,2),
                       stations.staid.lon.lat=c(1,3,2)
)

stfdf = rm.dupl(stfdf, zcol = 1, zero.tol = 0.6) # 1 stations removed

nrowsp <- length(stfdf@sp)

numNA <- apply(matrix(stfdf@data[,"prcp"],
                      nrow=nrowsp,byrow=F), MARGIN=1,
               FUN=function(x) sum(is.na(x)))

# Remove stations out of covariates
rem <- numNA != daysNum
stfdf <-  stfdf[rem,drop=F] # same

# plot(stfdf@sp)
# plot(border, add=T)
summary(stfdf@data$prcp)
summary(as.factor(ifelse(stfdf@data$prcp==0,0,1)))

### SP OVERLAY DEM and TWI ###

r <- raster("../dem_twi/dem_cat.tif")     
e <- raster::extract(r,stfdf@sp)
stfdf@sp$dem=e

r <- raster("../dem_twi/twi_cat.tif")     
e <- raster::extract(r, stfdf@sp) 
stfdf@sp$twi=e

stfdf <- stfdf[!is.na(stfdf@sp$dem),] # 1 station removed

save(stfdf, file=paste(v, "_", min(years), "_", max(years), '.rda', sep=""))

### IMERG ###

# download imerg first (1_download_gpmdata.R) #

time=index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

registerDoParallel(cores=detectCores()-1)
imerg = foreach (i = 1:daysNum) %dopar% { 
  version <- c('06B', '05B', '04B')
  r <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V06A.1day.tif", sep = "")))
  if(inherits(r, "try-error")) {
    for (ver in version) {
      r <- try(raster(paste("../imerg/", "3B-HHR-L.MS.MRG.3IMERG.", days[i], "-S233000-E235959.1410.V", ver, ".1day.tif", sep = "")))
      if(inherits(r, "try-error")) {
        next
      } else {
        break
      }
    }
    if(inherits(r, "try-error")) {
      return(rep(NA, length(stfdf@sp)))
    }
  }
  
  r@file@nodatavalue = 9999
  e <- raster::extract(r,stfdf@sp)
  return(e)
}
stopImplicitCluster()
imerg <- do.call("cbind", imerg)
imerg <- as.vector(imerg)
# imerg = imerg * 2.4 #/ 10 * 24
imerg[imerg==29999] = NA
stfdf@data$imerg = imerg

nrowsp <- length(stfdf@sp)

numNA <- apply(matrix(stfdf@data[,"imerg"],
                      nrow=nrowsp,byrow=F), MARGIN=2,
               FUN=function(x) sum(is.na(x)))

# Remove stations out of covariates
rem <- numNA != nrowsp
stfdf <-  stfdf[, rem,drop=F] # same

# save(stfdf, file=paste(v, "_", min(years), "_", max(years), '.rda', sep=""))
load(file=paste(v, "_", min(years), "_", max(years), '.rda', sep=""))

### stations and imerg plot ###

r <- raster("../dem_twi/dem_cat.tif") 
r <- as(r, "SpatialPixelsDataFrame")

theme = theme_set(theme_minimal())
sta_dem_plot <- ggplot() + # watch out for attribute name color order
  geom_raster(data=as.data.frame(r), aes(x=x, y=y, fill=dem_cat), alpha=0.8) +
  scale_fill_gradientn(colours = terrain.colors(100), name = "DEM [m]") +
  geom_polygon(data = border, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
  geom_point(data = as.data.frame(stfdf), aes(x = lon, y = lat, color = "red", shape = as.factor("staid")), size = 0.5) +
  theme(plot.title = element_text(hjust = 8),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        text = element_text(size = 8),
        # legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=8, face="bold"),
        legend.text=element_text(size=8)) +
  labs(x = "Longitude", y = "Latitude") +
  scale_colour_manual(name = "Stations",
                      labels = c("GHCN-daily"),
                      values = c("red"="red")) +
  scale_shape_manual(name = "Stations",
                     labels = c("GHCN-daily"),
                     values = c(17)) +
  scalebar(x.min = 0.1594133, x.max = 3.322251,
           y.min = 40.5230524, y.max = 42.861523,
           st.size = 3.5, location="bottomright", st.dist=0.07, border.size=0.5,
           dist = 75, dist_unit = "km",
           transform = T, model = "WGS84",
           anchor=c(x=3.1, y=40.6)) +
  north(x.min = 0.1594133, x.max = 3.322251,
        y.min = 40.5230524, y.max = 42.861523,
        scale = 0.2, symbol = 3, location="bottomright",
        anchor=c(x=3.1, y=40.9))

### Fig7 ###
# tiff("../plot/stations.tiff", width = 100, height = 70, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/stations.jpeg", width = 100, height = 70, units = 'mm', res = 600)
sta_dem_plot
dev.off()

###### TMAX TMIN ###################################################################

load(paste(v, "_", min(years), "_", max(years), '.rda', sep=""))
summary(stfdf@data$prcp)
length(stfdf@data$prcp)
summary(stfdf@data$prcp==0)

time=index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

registerDoParallel(cores=detectCores()-1)
tmax = foreach (i = 1:daysNum) %dopar% { 
  # for (i in 1:daysNum) {
  tmax <- raster(paste("../max/", days[i], ".tif", sep = ""))
  tmax <- raster::extract(tmax,stfdf@sp)
  return(tmax)
}
stopImplicitCluster()
tmax <- do.call("cbind", tmax) # check this
tmax <- as.vector(tmax)
stfdf@data$tmax = tmax / 10

registerDoParallel(cores=detectCores()-1)
tmin = foreach (i = 1:daysNum) %dopar% { 
  # for (i in 1:daysNum) {
  tmin <- raster(paste("../min/", days[i], ".tif", sep = ""))
  tmin <- raster::extract(tmin,stfdf@sp)
  return(tmin)
}
stopImplicitCluster()
tmin <- do.call("cbind", tmin) # check this
tmin <- as.vector(tmin)
stfdf@data$tmin = tmin / 10

stfdf@data$tdif <- stfdf@data$tmax - stfdf@data$tmin

save(stfdf, file=paste(v, "_", min(years), "_", max(years), '.rda', sep=""))

### tmin & tmax plot ###

for (day in days[1:4]) {
  tmin <- raster(paste("../min/", day, ".tif", sep = ""))
  values(tmin) <- values(tmin)/10
  tmin <- as(tmin, "SpatialPixelsDataFrame")
  names(tmin) <- "tmin"
  
  tmax <- raster(paste("../max/", day, ".tif", sep = ""))
  values(tmax) <- values(tmax)/10
  tmax <- as(tmax, "SpatialPixelsDataFrame")
  names(tmax) <- "tmax"
  
  im <- raster(paste('../imerg/3B-HHR-L.MS.MRG.3IMERG.', day, '-S233000-E235959.1410.V05B.1day.tif', sep = ''))
  im <- crop(im, raster::buffer(border, width = 0.1))
  im <- mask(im, raster::buffer(border, width = 0.1))
  # plot(im)
  # plot(border, add=T)
  im <- as(im, "SpatialPixelsDataFrame")
  im[[1]] <- im[[1]] / 10
  names(im) <- "im"
  
  assign(paste("tmin", day, sep=""), ggplot() + # watch out for attribute name color order
           geom_raster(data=as.data.frame(tmin), aes(x=x, y=y, fill=tmin), alpha=0.8) +
           scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "red"))(30), name = expression(paste("Tmin [",degree,"C]", sep = ""))) +
           geom_polygon(data = border, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
           theme(plot.title = element_text(hjust = 8),
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 text = element_text(size = 8),
                 # legend.key.size= unit(0.2, "cm"),
                 # legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=8), #, face="bold"),
                 legend.text=element_text(size=8)) +
           # theme(plot.title = element_text(hjust = 0.5)) +
           # theme(axis.text = element_text(size = 8)) + # changes axis labels
           # theme(axis.title = element_text(size = 8)) + # change axis titles
           # theme(text = element_text(size = 8))+
           labs(x = "Longitude", y = "Latitude") +
           scale_x_continuous(limits = c(im@bbox[1,1], im@bbox[1,2])) +
           scale_y_continuous(limits = c(im@bbox[2,1], im@bbox[2,2]))
  )
  
  assign(paste("tmax", day, sep=""), ggplot() + # watch out for attribute name color order
           geom_raster(data=as.data.frame(tmax), aes(x=x, y=y, fill=tmax), alpha=0.8) +
           scale_fill_gradientn(colours = colorRampPalette(c("orange", "red"))(30), name = expression(paste("Tmax [",degree,"C]", sep = ""))) +
           geom_polygon(data = border, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
           theme(plot.title = element_text(hjust = 8),
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 text = element_text(size = 8),
                 # legend.key.size= unit(0.2, "cm"),
                 # legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=8), #, face="bold"),
                 legend.text=element_text(size=8)) +
           # theme(plot.title = element_text(hjust = 0.5)) +
           # theme(axis.text = element_text(size = 8)) + # changes axis labels
           # theme(axis.title = element_text(size = 8)) + # change axis titles
           # theme(text = element_text(size = 8))+
           labs(x = "Longitude", y = "Latitude") +
           scale_x_continuous(limits = c(im@bbox[1,1], im@bbox[1,2])) +
           scale_y_continuous(limits = c(im@bbox[2,1], im@bbox[2,2]))
  )
  
  assign(paste("im", day, sep=""), ggplot() + # watch out for attribute name color order
           geom_raster(data=as.data.frame(im), aes(x=x, y=y, fill=im), alpha=0.8) +
           scale_fill_continuous(name = "IMERG\n[mm/day]", high = "#132B43", low = "#56B1F7") +
           geom_polygon(data = border, aes(x = long, y = lat, group = group), alpha = 0, color = "black", fill=NA, size = 0.1) +
           theme(plot.title = element_text(hjust = 8),
                 axis.text = element_text(size = 8),
                 axis.title = element_text(size = 8),
                 text = element_text(size = 8),
                 # legend.key.size= unit(0.2, "cm"),
                 # legend.margin = unit(0, "cm"),
                 legend.title = element_text(size=8), #, face="bold"),
                 legend.text=element_text(size=8)) +
           labs(x = "Longitude", y = "Latitude") +
           scale_x_continuous(limits = c(im@bbox[1,1], im@bbox[1,2])) +
           scale_y_continuous(limits = c(im@bbox[2,1], im@bbox[2,2])))
  if (day==days[1]) {
    im20160101 <- im20160101 + north(x.min = 0.1000027, x.max = 3.400003,
                                     y.min = 40.3999993, y.max = 42.999999,
                                     scale = 0.4, symbol = 3, location="bottomright")
  }
  if (day==days[4]) {
    im20160104 <- im20160104 + scalebar(x.min = 0.1594133, x.max = 3.322251,
                                        y.min = 40.5230524, y.max = 42.861523,
                                        st.size = 3.5, location="bottomright", st.dist=0.07, border.size=0.5,
                                        dist = 75, dist_unit = "km",
                                        transform = T, model = "WGS84",
                                        anchor=c(x=3.1, y=40.6))
  }
}

### Fig9 ###
# tiff("../plot/tmax_tmin_im.tiff", width = 174, height = 232, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/tmax_tmin_im.jpeg", width = 174, height = 232, units = 'mm', res = 600)
# ggarrange(tmax20160101, tmin20160101, im20160101,
#           tmax20160102, tmin20160102, im20160102,
#           tmax20160103, tmin20160103, im20160103,
#           tmax20160104, tmin20160104, im20160104,
#           ncol=3, nrow=4)#, common.legend = TRUE, legend="bottom")
f1 <- ggarrange(ggarrange(tmax20160101, tmax20160102, tmax20160103, tmax20160104, ncol=1, nrow=4, common.legend = TRUE, legend="bottom"),
                    # labels = c("1 Jan 2016", "2 Jan 2016", "3 Jan 2016", "4 Jan 2016")),
          ggarrange(tmin20160101, tmin20160102, tmin20160103, tmin20160104, ncol=1, nrow=4, common.legend = TRUE, legend="bottom"),
          ggarrange(im20160101, im20160102, im20160103, im20160104, ncol=1, nrow=4, common.legend = TRUE, legend="bottom"),
          ncol=3, nrow=1)#,
          # labels = c("Max. temperature", "Min. temperature", "IMERG"))#, common.legend = TRUE, legend="bottom")
annotate_figure(f1,
                top = text_grob("Max. temperature                            Min. temperature                                    IMERG", face = "bold", size = 10),
                left = text_grob("               1 Jan 2016                                    2 Jan 2016                                     3 Jan 2016                                     4 Jan 2016", face = "bold", rot = 90, size = 10)
)
dev.off()

##### Correlation ################################################
load(paste(v, "_", min(years), "_", max(years), '.rda', sep=""))
cor(stfdf@data$prcp, stfdf@data$imerg, use = "pairwise.complete.obs", method = "pearson")
# cor(stfdf@data$prcp[stfdf@data$prcp>0], stfdf@data$imerg[stfdf@data$prcp>0], use = "pairwise.complete.obs", method = "pearson")
cor(stfdf@data$prcp, stfdf@data$tmax, use = "pairwise.complete.obs", method = "pearson")
cor(stfdf@data$prcp, stfdf@data$tmin, use = "pairwise.complete.obs", method = "pearson")
###### create 5 space-time folds ###################################################################

load(paste(v, "_", min(years), "_", max(years), '.rda', sep=""))
hist(stfdf@data$prcp)

### Fig8 ###
# tiff("../plot/histogram.tiff", width = 80, height = 75, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/histogram.jpeg", width = 80, height = 75, units = 'mm', res = 600)
par(cex = 0.7, mar=c(4,4,0.5,0.5))
res_hist = hist(stfdf$prcp, breaks=50, plot=F)
cuts = cut(res_hist$breaks, c(-Inf, 0, Inf))
# multiplier <- res_hist$counts / res_hist$density
# hist_density <- density(stfdf$prcp, na.rm=T)
# hist_density$y <- hist_density$y * multiplier[1]
# lines(hist_density, col="red", lwd=2) 
plot(res_hist, col=c("red","white")[cuts], main=NULL, ylim=c(0, 5000), xlab="Precipitation [mm]")
text(x=35, y=4800, labels=paste(res_hist$counts[1], "\n(< 5 mm)", sep=""), col="red")
text(x=150, y=4600, labels="3rd Quartile: 0.2 mm\nMean: 2.0 mm\nMaximum: 220.9 mm")
dev.off()

source("../stratfolds.R")
test = as.data.frame(stfdf@sp)
test$ID = index(stfdf@sp)
test$weight <-rep(1,dim(test)[1])
strat <- stratfold3d(target.name = c("lon", "lat"), # stavi prcp
                     other.names = c("lon", "lat"),
                     data = test,
                     num.folds = 5,
                     num.means = 5,
                     seed = 42,
                     cum.prop = 0.9)

# set.seed(42)
# indices <- CreateSpacetimeFolds(test ,spacevar="ID", k=5)
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


###### TMAX TMIN IMERG ######
###### create STRK model ###################################################################

load(paste(v, "_", min(years), "_", max(years), '.rda', sep=""))
stfdf$tdif = NULL

temp_df = as.data.frame(stfdf)
temp_df <- temp_df[, c(1:2, 4, 11:14)]

temp_df$completed = complete.cases(temp_df)
### log transformation ###
# temp_df$log <- log1p(temp_df$prcp)
# temp_df$prcp <- log10(temp_df$prcp + 1)
# temp_df$imerg <- log10(temp_df$imerg + 1)
# temp_df$prcp[5]
# log1p(temp_df$prcp[5])
# expm1(log1p(temp_df$prcp[5]))
# 
# log10(temp_df$prcp[5]+1)
# 10^(log10(temp_df$prcp[5]+1)) - 1

summary(temp_df)
summary(complete.cases(temp_df))

### multiple linear regression - trend ###
temp_df.dev = temp_df[complete.cases(temp_df), ]

set.seed(42)
lm = lm(prcp ~ imerg + tmax + tmin, temp_df.dev)
coefficients(lm) # model coefficients, added to tregcoef$tmeanHR
summary(lm)

# confint(lm, level=0.95) # CIs for model parameters 
# anova(lm) # anova table 
# vcov(lm) # covariance matrix for model parameters

# save(lm, file='../models/STRK_lm.rda')
load(file='../models/STRK_lm.rda')

rmse = sqrt(sum((temp_df.dev$prcp - lm$fitted.values)^2)/(length(temp_df.dev$prcp)))
# 5.288771
tss = t(temp_df.dev$prcp - mean(temp_df.dev$prcp)) %*% (temp_df.dev$prcp - mean(temp_df.dev$prcp))
ess = t(temp_df.dev$prcp - lm$fitted.values) %*% (temp_df.dev$prcp - lm$fitted.values)
r2 = (tss-ess)/tss
# 0.4090632

lm_trend = c()
br = 1
for (i in 1:length(temp_df$completed)) {
  if (temp_df$completed[i]){
    lm_trend[i] = lm$fitted.values[br]
    br = br + 1
  }
}
stfdf$lm_trend = lm_trend

### scaterplot ###
plot(stfdf@data$prcp ~ stfdf@data$lm_trend)
abline(v=0, col="red")

### STRK - fitting pooled variogram ###

stfdf$lm_res = stfdf$prcp - stfdf$lm_trend
# save(stfdf, file = "stfdf_temp.rda")
load(file = "stfdf_temp.rda")

library(e1071)
skewness(stfdf$prcp[!is.na(stfdf$prcp)]) # 7.027869
skewness(log10(stfdf$prcp[!is.na(stfdf$prcp)]+1)) # 2.29061
skewness(stfdf$lm_res[!is.na(stfdf$lm_res)]) # 5.669737 # log 0.9175563

summary(stfdf$lm_res)
length(stfdf$lm_res[!is.na(stfdf$lm_res)]) # 92320
length(stfdf$lm_res[abs(stfdf$lm_res)>20 & !is.na(stfdf$lm_res)]) # 1273
length(stfdf$lm_res[stfdf$lm_res>20 & !is.na(stfdf$lm_res)]) # 1103
length(stfdf$lm_res[stfdf$lm_res<(-20) & !is.na(stfdf$lm_res)]) # 170

### Fig10 ###
# tiff("../plot/hist_res.tiff", width = 80, height = 75, units = 'mm', res = 600)
jpeg("../plot/hist_res.jpeg", width = 80, height = 75, units = 'mm', res = 600)
par(cex = 0.7, mar=c(4,4,0.5,0.5))
res_hist = hist(stfdf$lm_res, main=NULL, xlab="Residuals", breaks=300, xlim=c(-20,20))
# multiplier <- res_hist$counts / res_hist$density
# hist_density <- density(stfdf$lm_res, na.rm=T)
# hist_density$y <- hist_density$y * multiplier[1]
# lines(hist_density, col="red", lwd=3)
dev.off()

summary(stfdf$lm_res[stfdf$prcp==0])
# jpeg("../plot/hist_res0.jpeg", width = 80, height = 75, units = 'mm', res = 600)
par(cex = 0.7, mar=c(4,4,1,0.5))
res_hist = hist(stfdf$lm_res[stfdf$prcp==0], main=NULL, xlab="Residuals", breaks=120) #, xlim=c(-20,20))
# dev.off()

var = variogramST(lm_res ~ 1, stfdf) # , tlags = 0:5, cutoff = 300, width = 10, na.omit=T) # tunit="days"
# save(var, file=paste('../models/STRK_vgm', '.rda', sep=""))
load(file=paste('../models/STRK_vgm', '.rda', sep=""))

plot(var, map = F)#, main="2d sample variogram")

plot(var, wireframe=T, zlim=c(0,35),
     zlab=NULL,
     xlab=list("distance (km)", rot=30, cex = 0.7),
     ylab=list("time lag (days)", rot=-35, cex = 0.7),
     scales=list(arrows=F, z = list(distance = 5), cex = 0.5)) #, main="3d sample variogram")
dev.off()

var$dist = var$dist/10
var$spacelag = var$spacelag/10
var$avgDist = var$avgDist/10
attr(var, "boundaries") = attr(var, "boundaries") / 10
# var$timelag = var$timelag/24
# attr(var$timelag, "units") = "days"

estiStAni(var, c(0, 30), "metric",
          vgm(psill=12,"Sph", range=3, nugget=0),
          vgm(psill=25,"Sph", range=2, nugget=0) )

pars.l <- c(sill.s = 0, range.s = 1, nugget.s = 0,
            sill.t = 0, range.t = 1, nugget.t = 0,
            sill.st = 0, range.st = 1, nugget.st = 0,
            anis = 0)

sumMetric <- vgmST("sumMetric",
                   space = vgm(psill=12,"Sph", range=3, nugget=5),
                   time = vgm(psill=25,"Sph", range=2, nugget=5),
                   joint = vgm(psill=20,"Sph", range=3, nugget=5),
                   stAni=12)
# sumMetric <- vgmST("sumMetric",
#                    space = vgm(psill=12,"Exp", range=3, nugget=0),
#                    time = vgm(psill=25,"Sph", range=5, nugget=0),
#                    joint = vgm(psill=20,"Exp", range=5, nugget=0),
#                    stAni=12)

sumMetric_Vgm <- fit.StVariogram(var, sumMetric, method="L-BFGS-B",lower=pars.l)
attr(sumMetric_Vgm, "MSE") # 1.491033
sumMetric_Vgm

var$dist <- var$dist*10
var$spacelag <- var$spacelag*10
var$avgDist = var$avgDist*10
attr(var, "boundaries") = attr(var, "boundaries") * 10

sumMetric_Vgm$space$range <- sumMetric_Vgm$space$range*10
sumMetric_Vgm$joint$range <- sumMetric_Vgm$joint$range*10
sumMetric_Vgm$stAni <- sumMetric_Vgm$stAni*10

# plot(var, sumMetric_Vgm,map=F)#, main="2d fitted sum-metric variogram")

# save(sumMetric_Vgm, file=paste('../models/STRK_fit_vgm', '.rda', sep=""))
load(file=paste('../models/STRK_vgm', '.rda', sep=""))
load(file=paste('../models/STRK_fit_vgm', '.rda', sep=""))

### Fig11 ###
# tiff("../plot/fitted3d_variogram.tiff", width = 130, height = 70, units = 'mm', res = 600)
jpeg("../plot/fitted3d_variogram.jpeg", width = 130, height = 70, units = 'mm', res = 600)
par(cex = 0.7, mar=c(0.5,0.5,0.5,0.5))
plot(var, list(sumMetric_Vgm),all=T, wireframe=T, zlim=c(0,35),
     zlab=NULL,
     xlab=list("distance (km)", rot=30, cex = 0.7),
     ylab=list("time lag (days)", rot=-35, cex = 0.7),
     scales=list(arrows=F, z = list(distance = 5), cex = 0.5))#, main="3d fitted sum-metric variogram")
dev.off()

###### STRK - 5-fold CV ###################################################################

# UK fit_vgm_auto <- autofitVariogram(prcp ~ imerg + tmax + tmin + time, df1, model = c("Sph", "Exp"), dX = 0)
load(file = "stfdf_temp.rda")
load(file="folds.rda")
load(file='../models/STRK_fit_vgm.rda')

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

STRK_5f_obs <- c()
STRK_5f_pred <- c()

for(val_fold in 1:5){
  print(paste("Fold ", val_fold, sep = ""))
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  dev = stfdf[stfdf@sp$ind==1, ]
  val = stfdf[stfdf@sp$ind==2, ]
  
  ### fiting STRK model on dev data ###
  temp_df = as.data.frame(dev)
  temp_df <- temp_df[, c(1:2, 4, 11:15)]
  temp_df$completed = complete.cases(temp_df)
  temp_df.dev = temp_df[complete.cases(temp_df), ]
  set.seed(42)
  lm = lm(prcp ~ imerg + tmax + tmin, temp_df.dev)
  # coefficients(lm) # model coefficients, added to tregcoef$tmeanHR
  # summary(lm)
  lm_trend = c()
  br = 1
  for (i in 1:length(temp_df$completed)) {
    if (temp_df$completed[i]){
      lm_trend[i] = lm$fitted.values[br]
      br = br + 1
    } else {
      lm_trend[i] = NA
    }
  }
  dev$lm_trend = lm_trend
  dev$lm_res = dev$prcp - dev$lm_trend
  df1 <- as.data.frame(dev)
  df1 <- df1[, c(1:2, 4, 17)]
  df1 <- df1[complete.cases(df1), ]
  coordinates(df1) <- ~ lon + lat
  df1@proj4string = stfdf@sp@proj4string
  
  var = variogramST(lm_res ~ 1, dev) # , tlags = 0:5, cutoff = 300, width = 10, na.omit=T) # tunit="days"
  
  var$dist = var$dist/10
  var$spacelag = var$spacelag/10
  var$avgDist = var$avgDist/10
  attr(var, "boundaries") = attr(var, "boundaries") / 10
  
  pars.l <- c(sill.s = 0, range.s = 1, nugget.s = 0,
              sill.t = 0, range.t = 1, nugget.t = 0,
              sill.st = 0, range.st = 1, nugget.st = 0,
              anis = 0)
  
  sumMetric <- vgmST("sumMetric",
                     space = vgm(psill=12,"Sph", range=3, nugget=5),
                     time = vgm(psill=25,"Sph", range=2, nugget=5),
                     joint = vgm(psill=20,"Sph", range=3, nugget=5),
                     stAni=12)
  
  sumMetric_Vgm <- fit.StVariogram(var, sumMetric, method="L-BFGS-B",lower=pars.l)
  attr(sumMetric_Vgm, "MSE")
  sumMetric_Vgm
  
  # var$dist <- var$dist*10
  # var$spacelag <- var$spacelag*10
  # var$avgDist = var$avgDist*10
  # attr(var, "boundaries") = attr(var, "boundaries") * 10
  
  sumMetric_Vgm$space$range <- sumMetric_Vgm$space$range*10
  sumMetric_Vgm$joint$range <- sumMetric_Vgm$joint$range*10
  sumMetric_Vgm$stAni <- sumMetric_Vgm$stAni*10
  
  # plot(var, list(sumMetric_Vgm),all=T, wireframe=T, zlim=c(0,35),
  #      zlab=NULL,
  #      xlab=list("distance (km)", rot=30, cex = 0.7),
  #      ylab=list("time lag (days)", rot=-35, cex = 0.7),
  #      scales=list(arrows=F, z = list(distance = 5), cex = 0.5))#, main="3d fitted sum-metric variogram")
  
  ###
  
  val$lm_trend <- predict(lm, val@data)
  val$lm_res <- val$prcp - val$lm_trend
  
  sp.nmax = 10
  i_1=c(rep(1,1),1:(daysNum -1))
  ip1=c(1:daysNum)
  vario = sumMetric_Vgm
  temp = val
  
  N_POINTS <- length(temp@sp@coords[,1])
  
  registerDoParallel(cores=detectCores()-1)
  cv <- foreach(i = 1:N_POINTS, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
    st= dev@sp
    st$dist=spDists(dev@sp,temp@sp[i,])
    tmp_st<-st[ order(st$'dist') , ]
    # do not remove target station (now is added because of independent CV)
    local_t= row.names(tmp_st[1:sp.nmax,] )
    
    xxx = as.list ( rep(NA, length(time) ) )
    for( ii in 1:length(time) ) {
      data=dev[local_t, i_1[ii]:ip1[ii],'lm_res',drop=F]
      nrowsp <- length(data@sp)
      # count NAs per stations
      numNA <- apply(matrix(data@data[,'lm_res'],
                            nrow=nrowsp,byrow=F), MARGIN=1,
                     FUN=function(x) sum(is.na(x)))
      # Remove stations out of covariates
      rem <- !numNA > 0
      data <-  data[rem,drop=F]
      
      # If there are less than 5 observations
      if (length(data)<5){
        for (m in 2:5){
          local_t1= row.names(tmp_st[1:(sp.nmax*m),] )
          data=dev[local_t1, i_1[ii]:ip1[ii],'lm_res']
          nrowsp <- length(data@sp)
          # count NAs per stations
          numNA <- apply(matrix(data@data[,'lm_res'],
                                nrow=nrowsp,byrow=F), MARGIN=1,
                         FUN=function(x) sum(is.na(x)))
          # Remove stations out of covariates
          rem <- !numNA > 0
          data <-  data[rem,drop=F]
          if (length(data)>5) break
        }
      }
      
      if (length(data)==0){
        xxx[[ii]] <- NA
      } else {
        xxx[[ii]]=krigeST(as.formula("lm_res~1"),
                          data=data, 
                          newdata=STF(as(temp@sp[i,],"SpatialPoints"),
                                      temp@time[ii],  
                                      temp@endTime[ii]),     
                          modelList=vario,
                          computeVar=FALSE)@data[,1]
      }
      
    } # end of  for
    ret=unlist(xxx) 
  }
  stopImplicitCluster()
  
  cv <- do.call(rbind,cv)
  cv <- as.vector(cv)
  cv.temp <- temp
  
  cv.temp$pred.cv <- cv + cv.temp$lm_trend
  cv.temp$resid.cv <- cv.temp$prcp  - cv.temp@data$pred.cv
  
  STRK_5f_obs <- c(STRK_5f_obs, cv.temp$prcp)
  STRK_5f_pred <- c(STRK_5f_pred, cv.temp$pred.cv)
  
}

# save(STRK_5f_pred, file = "STRK_5f_pred.rda")
# save(STRK_5f_obs, file = "STRK_5f_obs.rda")
load(file = "STRK_5f_pred.rda")
load(file = "STRK_5f_obs.rda")

summary(STRK_5f_pred)
STRK_5f_pred <- ifelse(STRK_5f_pred<0, 0, STRK_5f_pred)

STRK_5f_pred <- STRK_5f_pred[!is.na(STRK_5f_obs)]
STRK_5f_obs <- STRK_5f_obs[!is.na(STRK_5f_obs)]

## RMSE
sqrt(mean((STRK_5f_obs - STRK_5f_pred)^2, na.rm = T))
# 3.922399
## CCC
DescTools::CCC(STRK_5f_obs, STRK_5f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.8145016
## MAE
mean(abs(STRK_5f_obs - STRK_5f_pred), na.rm=TRUE)
# 1.19602
## R2
1 - (t(STRK_5f_obs - STRK_5f_pred) %*% (STRK_5f_obs - STRK_5f_pred)) / (t(STRK_5f_obs - mean(STRK_5f_obs)) %*% (STRK_5f_obs - mean(STRK_5f_obs)))
# 0.6749614
## ME
mean((STRK_5f_obs - STRK_5f_pred), na.rm=TRUE)
# -0.09174077

###### RFSI - 5-fold CV ###################################################################

load(file = "stfdf_temp.rda")
load(file="folds.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

rfsi_5f_obs <- c()
rfsi_5f_pred <- c()

nos <- c(6,8,10,12,14)#6:15#c(3,4,5,6,7,8,9,10,13,15,17,20)
ntree <- 250 # 500

n_obs = max(nos)

for(val_fold in 1:5){
  cat("fold ")
  cat(val_fold)
  cat("\n")
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  temp_df <- temp_df[, c(1:4, 11:15)]
  temp_df = temp_df[complete.cases(temp_df), ]
  
  # time <- unique(temp_df$time)
  cpus <- detectCores()-1
  registerDoParallel(cores=cpus)
  nearest_obs <- foreach (t = time) %dopar% {
    
    day_df <- temp_df[temp_df$time==t, c("lon", "lat", "prcp", "ind")]
    
    if (nrow(day_df)==0) {
      return 
    }
    return(near.obs(
      locations = day_df,
      observations = day_df[day_df$ind==1, ],
      zcol = "prcp",
      n.obs = n_obs
    ))
    
  }
  stopImplicitCluster()
  
  nearest_obs <- do.call("rbind", nearest_obs)
  temp_df.dev <- cbind(temp_df, nearest_obs)
  
  temp_df.dev = temp_df.dev[complete.cases(temp_df.dev), ]
  
  # dev = stfdf[stfdf@sp$ind==1, ]
  # val = stfdf[stfdf@sp$ind==2, ]
  dev_df <- temp_df.dev[temp_df.dev$ind==1, -5]
  val_df <- temp_df.dev[temp_df.dev$ind==2, -5]
  
  ### Tune RFSI ###
  indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                  k=5, seed = 42)
  trControl=trainControl(method="cv", number = 5,
                                     verboseIter = T, #classProbs = T,
                                     index = indices$index,
                         search = "random")
  train_rmse <- 100
  br <- 0
  
  for (no in nos){
    cat(no)
    cat("\n")
    dev_df1 <- dev_df[, c(5:(9+no-1), (9+n_obs):(9+n_obs+no-1))]
    # tgrid <- expand.grid(
    #   .mtry = 10, #c(9:11),
    #   .splitrule = c("variance"), #, "extratrees"), #, "maxstat"),
    #   .min.node.size = c(15, 20, 25)
    # )
    model_caret <- caret::train(prcp ~ ., data = dev_df1,
                         method = "ranger",
                         metric = "RMSE",
                         trControl = trControl,
                         # tuneGrid = tgrid,
                         num.trees = ntree,
                         importance = "impurity",
                         tuneLength = 20,
                         seed = 42)
    current_rmse <- model_caret$results[as.integer(row.names(model_caret$bestTune)), "RMSE"]
    # print(model_caret$results[row.names(model_caret$results) == row.names(model_caret$bestTune), ])
    # cat("\n")
    
    if (current_rmse < train_rmse){
      # dev_model <- model_caret$finalModel
      dev_parameters <- model_caret$bestTune
      dev_obs <- no
      train_rmse <- current_rmse
      br <- 0
    } else {
      br <- br + 1
      if (br==3) break # if 3 times in a row RMSE is not lower
    }
  }
  
  print(dev_parameters)
  cat("\n")
  
  dev_df1 <- dev_df1 <- dev_df[, c(5:(9+dev_obs-1), (9+n_obs):(9+n_obs+dev_obs-1))]
  dev_model <- ranger(prcp ~ ., data = dev_df1, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = dev_parameters$mtry,
                      splitrule = dev_parameters$splitrule,
                      min.node.size = dev_parameters$min.node.size)
  
  ### validation ###
  val_df1 <- val_df[, c(6:(9+dev_obs-1), (9+n_obs):(9+n_obs+dev_obs-1))]
  rfsi_pred <- predict(dev_model, val_df1)
  # summary(rfsi_pred$predictions)
  # rfsi_res <- val_df$prcp - rfsi_pred$predictions
  
  rfsi_5f_obs <- c(rfsi_5f_obs, val_df$prcp)
  rfsi_5f_pred <- c(rfsi_5f_pred, rfsi_pred$predictions)
}

# save(rfsi_5f_pred, file="rfsi_5f_pred.rda")
# save(rfsi_5f_obs, file="rfsi_5f_obs.rda")

load(file="rfsi_5f_pred.rda")
load(file="rfsi_5f_obs.rda")

rfsi_5f_pred <- rfsi_5f_pred[!is.na(rfsi_5f_obs)]
rfsi_5f_obs <- rfsi_5f_obs[!is.na(rfsi_5f_obs)]

## RMSE
sqrt(mean((rfsi_5f_obs - rfsi_5f_pred)^2, na.rm = T))
# 3.802413
## CCC
DescTools::CCC(rfsi_5f_obs, rfsi_5f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.8178055
## MAE
mean(abs(rfsi_5f_obs - rfsi_5f_pred), na.rm=TRUE)
# 1.120374
## R2
1 - (t(rfsi_5f_obs - rfsi_5f_pred) %*% (rfsi_5f_obs - rfsi_5f_pred)) / (t(rfsi_5f_obs - mean(rfsi_5f_obs)) %*% (rfsi_5f_obs - mean(rfsi_5f_obs)))
# 0.6945429
## ME
mean((rfsi_5f_obs - rfsi_5f_pred), na.rm=TRUE)
# -0.04959614

###### RF - 5-fold CV ###################################################################

load(file = "stfdf_temp.rda")
load(file="folds.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

rf_5f_obs <- c()
rf_5f_pred <- c()

ntree <- 250 # 500

for(val_fold in 1:5){
  cat("fold ")
  cat(val_fold)
  cat("\n")
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  temp_df <- temp_df[, c(1:4, 11:15)]
  temp_df.dev = temp_df[complete.cases(temp_df), ]
  
  dev_df <- temp_df.dev[temp_df.dev$ind==1, -5]
  val_df <- temp_df.dev[temp_df.dev$ind==2, -5]
  
  indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                  k=5, seed = 42)
  
  trControl=trainControl(method="cv", number = 5,
                         verboseIter = T, #classProbs = T,
                         index = indices$index,
                         search = "random")
  
    # tgrid <- expand.grid(
    #   .mtry = 2:4,
    #   .splitrule = "gini",
    #   .min.node.size = c(10, 20)
    #   # .sample.fraction
    # )
    
    dev_df1 <- dev_df[, c(5:length(dev_df))]
    
    model_caret <- caret::train(prcp ~ ., data = dev_df1,
                                method = "ranger",
                                metric = "RMSE",
                                trControl = trControl,
                                # tuneGrid = tgrid,
                                num.trees = ntree,
                                importance = "impurity",
                                tuneLength = 10,
                                seed = 42)
    dev_parameters <- model_caret$bestTune
    
    dev_model <- ranger(prcp ~ ., data = dev_df1, importance = "impurity", seed = 42,
                  num.trees = ntree, mtry = dev_parameters$mtry,
                  splitrule = dev_parameters$splitrule,
                  min.node.size = dev_parameters$min.node.size)
  
  ### validation ###
  val_df1 <- val_df[, c(5:length(val_df))]
  rf_pred <- predict(dev_model, val_df1)
  # summary(rfsi_pred$predictions)
  # rf_res <- val_df1$prcp - rf_pred$predictions
  # summary(rf_res)
  # summary(rf_pred$predictions)
  
  rf_5f_obs <- c(rf_5f_obs, val_df$prcp)
  rf_5f_pred <- c(rf_5f_pred, rf_pred$predictions)
  
}

# save(rf_5f_pred, file="rf_5f_pred.rda")
# save(rf_5f_obs, file="rf_5f_obs.rda")

load(file="rf_5f_pred.rda")
load(file="rf_5f_obs.rda")

rf_5f_pred <- rf_5f_pred[!is.na(rf_5f_obs)]
rf_5f_obs <- rf_5f_obs[!is.na(rf_5f_obs)]

## RMSE
sqrt(mean((rf_5f_obs - rf_5f_pred)^2, na.rm = T))
# 4.821719
## CCC
DescTools::CCC(rf_5f_obs, rf_5f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.6694631
## MAE
mean(abs(rf_5f_obs - rf_5f_pred), na.rm=TRUE)
# 1.66342
## R2
1 - (t(rf_5f_obs - rf_5f_pred) %*% (rf_5f_obs - rf_5f_pred)) / (t(rf_5f_obs - mean(rf_5f_obs)) %*% (rf_5f_obs - mean(rf_5f_obs)))
# 0.5088261
## ME
mean((rf_5f_obs - rf_5f_pred), na.rm=TRUE)
# 0.01440828

###### RFsp- 5-fold CV ###################################################################

load(file = "stfdf_temp.rda")
load(file="folds.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

### Create BBOX ###
xmin=min(stfdf@sp@coords[, "lon"])-1000; xmax=max(stfdf@sp@coords[, "lon"])+1000;
ymin=min(stfdf@sp@coords[, "lat"])-1000 ; ymax=max(stfdf@sp@coords[, "lat"])+1000 # Serbia (proj bbox is 166021.4431, 833978.5569, 0.0000, 9329005.1825)
bbox = extent(xmin, xmax, ymin, ymax)
bbox <- as(bbox, "SpatialPolygons")
bbox@proj4string = CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
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

rfsp_5f_obs <- c()
rfsp_5f_pred <- c()

## Buffer distances only:
grid.distP <- GSIF::buffer.dist(stfdf@sp["staid"], fishnet[1], as.factor(1:nrow(stfdf@sp)))

for(val_fold in 1:5){
  cat(val_fold)
  val_ids = as.vector(unlist(strat$obs.fold.list[val_fold]))
  dev_ids = as.vector(unlist(strat$obs.fold.list[-val_fold]))
  stfdf@sp$ind = ifelse(index(stfdf@sp) %in% dev_ids, 1, ifelse(index(stfdf@sp) %in% val_ids, 2, NA))
  
  temp_df = as.data.frame(stfdf)
  
  ov_toma <- do.call(cbind, list(stfdf@sp@data, over(stfdf@sp, grid.distP)))
  temp_df <- plyr::join(temp_df, ov_toma, by="staid")
  
  temp_df <- temp_df[, c(1:4, 11:15, 22:length(temp_df))]
  
  temp_df.dev = temp_df[complete.cases(temp_df), ]
  
  dev_df <- temp_df.dev[temp_df.dev$ind==1, -5]
  val_df <- temp_df.dev[temp_df.dev$ind==2, -5]
  
  # RFsp model
  
  # indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
  #                                 k=5, seed = 42)
  # 
  # trControl=trainControl(method="cv", number = 5,
  #                        verboseIter = T, #classProbs = T,
  #                        index = indices$index) #,
  #                        # search = "random")
  # 
  dev_df1 <- dev_df[, c(5:length(dev_df))]
  dev_df1 <- dev_df1[, c(names(dev_df1)[1:4], paste("layer.", unlist(strat$obs.fold.list[-val_fold]), sep=""))]
  # 
  # tgrid <- expand.grid(
  #   .mtry = length(dev_df1)-11,
  #   .splitrule = c("variance", "extratrees"),
  #   .min.node.size = c(2,5,10,15,20)
  #   # .sample.fraction
  # )
  # 
  # splitrule   min.node.size  RMSE      Rsquared   MAE     
  # variance     2             4.592372  0.5387220  1.599585
  # variance     5             4.597378  0.5375311  1.602821
  # variance    10             4.601772  0.5363256  1.601028
  # variance    15             4.612040  0.5341427  1.601251
  # variance    20             4.622727  0.5319226  1.602650
  # extratrees   2             4.617584  0.5329733  1.628651
  # extratrees   5             4.640896  0.5283165  1.640885
  # extratrees  10             4.665543  0.5230956  1.646836
  # extratrees  15             4.690567  0.5181235  1.654341
  # extratrees  20             4.720987  0.5118995  1.660470
  # # Fitting mtry = 58, splitrule = variance, min.node.size = 2 on full training set
  # 
  # model_caret <- caret::train(prcp ~ ., data = dev_df1,
  #                             method = "ranger",
  #                             metric = "RMSE",
  #                             trControl = trControl,
  #                             tuneGrid = tgrid,
  #                             num.trees = ntree,
  #                             importance = "impurity",
  #                             # tuneLength = 10,
  #                             seed = 42)
  # 
  # # dev_model <- model_caret$finalModel
  # dev_parameters <- model_caret$bestTune
  
  dev_model <- ranger(prcp ~ ., data = dev_df1, importance = "impurity", seed = 42,
                      num.trees = ntree, mtry = length(dev_df1)-11,
                      splitrule = "variance",
                      min.node.size = 2)
  
  ### validation ###
  val_df1 <- val_df[, c(5:length(val_df))]
  val_df1 <- val_df1[, c(names(val_df1)[1:4], paste("layer.", unlist(strat$obs.fold.list[-val_fold]), sep=""))]
  
  rfsp_pred <- predict(dev_model, val_df1)
  # summary(rfsi_pred$predictions)
  # rfsp_res <- val_df1$prcp - rfsp_pred$predictions
  # summary(rfsp_res)
  # summary(rfsp_pred$predictions)
  
  # rmse_total = sqrt(sum((rfsp_res)^2, na.rm = T)/(length(rfsp_res[!is.na(rfsp_res)])))
  # # print(paste("Fold ", val_fold, ", RFnl RMSE", sep = ""))
  # print(rmse_total)
  
  rfsp_5f_obs <- c(rfsp_5f_obs, val_df$prcp)
  rfsp_5f_pred <- c(rfsp_5f_pred, rfsp_pred$predictions)
  
}

# save(rfsp_5f_pred, file="rfsp_5f_pred.rda")
# save(rfsp_5f_obs, file="rfsp_5f_obs.rda")
load(file="rfsp_5f_pred.rda")
load(file="rfsp_5f_obs.rda")

rfsp_5f_pred <- rfsp_5f_pred[!is.na(rfsp_5f_obs)]
rfsp_5f_obs <- rfsp_5f_obs[!is.na(rfsp_5f_obs)]

## RMSE
sqrt(mean((rfsp_5f_obs - rfsp_5f_pred)^2, na.rm = T))
# 4.559222
## CCC
DescTools::CCC(rfsp_5f_obs, rfsp_5f_pred, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c
# 0.7194699
## MAE
mean(abs(rfsp_5f_obs - rfsp_5f_pred), na.rm=TRUE)
# 1.555771
## R2
1 - (t(rfsp_5f_obs - rfsp_5f_pred) %*% (rfsp_5f_obs - rfsp_5f_pred)) / (t(rfsp_5f_obs - mean(rfsp_5f_obs)) %*% (rfsp_5f_obs - mean(rfsp_5f_obs)))
# 0.5608499
## ME
mean((rfsp_5f_obs - rfsp_5f_pred), na.rm=TRUE)
# -0.04914544

###### create RFSI model ###################################################################

load(file = "stfdf_temp.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

temp_df = as.data.frame(stfdf)
temp_df <- temp_df[, c(1:4, 11:14)]

dev_df = temp_df[complete.cases(temp_df), ]

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

nos <- c(6,8,10,12,14)
ntree <- 250 # 500

n_obs = max(nos)

cpus <- detectCores()-1
registerDoParallel(cores=cpus)
nearest_obs <- foreach (t = time) %dopar% {
  
  day_df <- dev_df[dev_df$time==t, c("lon", "lat", "prcp")]
  
  return(near.obs(
    locations = day_df,
    observations = day_df,
    zcol = "prcp",
    n.obs = n_obs
  ))
  
}
stopImplicitCluster()

nearest_obs <- do.call("rbind", nearest_obs)
dev_df <- cbind(dev_df, nearest_obs)

indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                k=5, seed = 42)

trControl=trainControl(method="cv", number = 5,
                       verboseIter = T, #classProbs = T,
                       index = indices$index,
                       search = "random")

train_rmse <- 100
br <- 0

for (no in nos){
  cat(no)
  cat("\n")
  
  dev_df1 <- dev_df[, c(5:(9+no-1), (9+n_obs):(9+n_obs+no-1))]
  
  # tgrid <- expand.grid(
  #   .mtry = 10, #c(9:11),
  #   .splitrule = c("variance"), #, "extratrees"), #, "maxstat"),
  #   .min.node.size = c(15, 20, 25)
  # )
  
  model_caret <- caret::train(prcp ~ ., data = dev_df1,
                              method = "ranger",
                              metric = "RMSE",
                              trControl = trControl,
                              # tuneGrid = tgrid,
                              num.trees = ntree,
                              importance = "impurity",
                              tuneLength = 20,
                              seed = 42)
  
  current_rmse <- model_caret$results[as.integer(row.names(model_caret$bestTune)), "RMSE"]
  # print(model_caret$results[row.names(model_caret$results) == row.names(model_caret$bestTune), ])
  # cat("\n")
  
  if (current_rmse < train_rmse){
    dev_parameters <- model_caret$bestTune
    dev_obs <- no
    train_rmse <- current_rmse
    br <- 0
  } else {
    br <- br + 1
    if (br==3) break
  }
  
}

print(dev_parameters)
#      mtry  splitrule   min.node.size   n_obs
# 18   12    extratrees  18              10

dev_df1 <- dev_df[, c(5:(9+dev_obs-1), (9+n_obs):(9+n_obs+dev_obs-1))]
rfsi_model <- ranger(prcp ~ ., data = dev_df1, importance = "impurity", seed = 42,
                    num.trees = ntree, mtry = dev_parameters$mtry,
                    splitrule = dev_parameters$splitrule,
                    min.node.size = dev_parameters$min.node.size,
                    quantreg = TRUE) ### quantreg???

# save(rfsi_model, file = "../models/RFSI.rda")
load(file = "../models/RFSI.rda")

rfsi_model$r.squared
# 0.7288497
sqrt(rfsi_model$prediction.error)
# 3.582545

###### create RF model ###################################################################

load(file = "stfdf_temp.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

temp_df = as.data.frame(stfdf)
temp_df <- temp_df[, c(1:4, 11:14)]

dev_df = temp_df[complete.cases(temp_df), ]

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

ntree <- 250

indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
                                k=5, seed = 42)

trControl=trainControl(method="cv", number = 5,
                       verboseIter = T, #classProbs = T,
                       index = indices$index,
                       search = "random")

dev_df1 <- dev_df[, c(5:length(dev_df))]
  
  # tgrid <- expand.grid(
  #   .mtry = 10, #c(9:11),
  #   .splitrule = c("variance"), #, "extratrees"), #, "maxstat"),
  #   .min.node.size = c(15, 20, 25)
  # )
  
  model_caret <- caret::train(prcp ~ ., data = dev_df1,
                              method = "ranger",
                              metric = "RMSE",
                              trControl = trControl,
                              # tuneGrid = tgrid,
                              num.trees = ntree,
                              importance = "impurity",
                              tuneLength = 10,
                              seed = 42)
  
  dev_parameters <- model_caret$bestTune
print(dev_parameters)
#   mtry splitrule min.node.size
# 1    2   maxstat             5

rf_model <- ranger(prcp ~ ., data = dev_df1, importance = "impurity", seed = 42,
                     num.trees = ntree, mtry = dev_parameters$mtry,
                     splitrule = dev_parameters$splitrule,
                     min.node.size = dev_parameters$min.node.size,
                     quantreg = TRUE) ### quantreg???

# save(rf_model, file = "../models/RF.rda")
load(file = "../models/RF.rda")

rf_model$r.squared
# 0.5097485
sqrt(rf_model$prediction.error)
# 4.817215

###### create RFsp model ###################################################################

load(file = "stfdf_temp.rda")

st_wgs <-stfdf@sp 
stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
st_proj <- stfdf@sp

### Create BBOX ###
xmin=min(stfdf@sp@coords[, "lon"])-1000; xmax=max(stfdf@sp@coords[, "lon"])+1000;
ymin=min(stfdf@sp@coords[, "lat"])-1000 ; ymax=max(stfdf@sp@coords[, "lat"])+1000 # Serbia (proj bbox is 166021.4431, 833978.5569, 0.0000, 9329005.1825)
bbox = extent(xmin, xmax, ymin, ymax)
bbox <- as(bbox, "SpatialPolygons")
bbox@proj4string = CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
fishnet <- raster(bbox)
res(fishnet) <- c(1000, 1000) # 1 x 1 km
fishnet[] <- runif(length(fishnet), -10, 10)
fishnet = as(fishnet, "SpatialPixelsDataFrame")
# plot(fishnet)
st <- as(fishnet, "SpatialPoints")

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

ntree <- 250

## Geographic distances only:
grid.distP <- GSIF::buffer.dist(stfdf@sp["staid"], fishnet[1], as.factor(1:nrow(stfdf@sp)))

temp_df = as.data.frame(stfdf)

ov_toma <- do.call(cbind, list(stfdf@sp@data, over(stfdf@sp, grid.distP)))
temp_df <- plyr::join(temp_df, ov_toma, by="staid")

temp_df <- temp_df[, c(1:4, 11:14, 20:length(temp_df))]

dev_df = temp_df[complete.cases(temp_df), ]

time=zoo::index(stfdf@time)
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

# indices <- CreateSpacetimeFolds(dev_df,spacevar = "sp.ID",
#                                 k=5, seed = 42)
# 
# trControl=trainControl(method="cv", number = 5,
#                        verboseIter = T, #classProbs = T,
#                        index = indices$index,
#                        search = "random")

dev_df1 <- dev_df[, c(5:length(dev_df))]

# tgrid <- expand.grid(
#   .mtry = 60,
#   .splitrule = c("variance"), #, "extratrees"), #, "maxstat"),
#   .min.node.size = c(15, 20, 25)
# )
#   
#   model_caret <- caret::train(prcp ~ ., data = dev_df1,
#                               method = "ranger",
#                               metric = "RMSE",
#                               trControl = trControl,
#                               # tuneGrid = tgrid,
#                               num.trees = ntree,
#                               importance = "impurity",
#                               tuneLength = 10,
#                               seed = 42)
#   
#   dev_parameters <- model_caret$bestTune
# print(dev_parameters)
# #   mtry  splitrule min.node.size
# # 1    60  variance             2

rfsp_model <- ranger(prcp ~ ., data = dev_df1, importance = "impurity", seed = 42,
                     num.trees = ntree, mtry = length(dev_df1)-11,
                     splitrule = "variance",
                     min.node.size = 2,
                     quantreg = TRUE) ### quantreg???

# save(rfsp_model, file = "../models/RFsp.rda")
load(file = "../models/RFsp.rda")

rfsp_model$r.squared
# 0.5298143
sqrt(rfsp_model$prediction.error)
# 4.717602

###### importance plot ###################################################################

load(file = "../models/RFSI.rda")

xlP.g <- as.list(round(rfsi_model$variable.importance))
df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:15,]/100000
print(df)
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
rfsi_importance <- 
  ggplot(df, aes(x=covariate, y=importance)) +
  geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
  geom_point( color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
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
df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:3,]/100000
print(df)
pr <- 3:1
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
  geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
  geom_point( color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
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
  geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
  geom_point( color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
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
ggarrange(rf_importance, rfsi_importance, rfsp_importance, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

###### CV error histogram, scatter plot ###################################################################

load(file = "STRK_5f_pred.rda")
load(file = "STRK_5f_obs.rda")
summary(STRK_5f_pred)
STRK_5f_pred <- ifelse(STRK_5f_pred<0, 0, STRK_5f_pred)

load(file = "rfsi_5f_obs.rda")
load(file = "rfsi_5f_pred.rda")

load(file = "rfsp_5f_obs.rda")
load(file = "rfsp_5f_pred.rda")

load(file = "rf_5f_obs.rda")
load(file = "rf_5f_pred.rda")

# # tiff("../plot/scaterplot.tiff", width = 174, height = 65, units = 'mm', res = 600, compression = "lzw")
# jpeg("../plot/scaterplot.jpeg", width = 174, height = 65, units = 'mm', res = 600)
# par(mfrow=c(1,4), cex = 0.7, mar=c(4.5,4.5,1,0.5))
# plot(STRK_5f_obs ~ STRK_5f_pred, main="STRK", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "Observations [mm]", cex = 0.5)
# abline(a=0, b=1, col="red")
# plot(rf_5f_obs ~ rf_5f_pred, main="RF", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "", cex = 0.5)
# abline(a=0, b=1, col="red")
# plot(rfsi_5f_obs ~ rfsi_5f_pred, main="RFSI", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "", cex = 0.5)
# abline(a=0, b=1, col="red")
# plot(rfsp_5f_obs ~ rfsp_5f_pred, main="RFsp", xlim=c(0, 120), xlab = "Predictions [mm]", ylab = "", cex = 0.5)
# abline(a=0, b=1, col="red")
# dev.off()

STRK_5f_pred <- STRK_5f_pred[!is.na(STRK_5f_obs)]
STRK_5f_obs <- STRK_5f_obs[!is.na(STRK_5f_obs)]
STRK_res <- STRK_5f_obs-STRK_5f_pred
# STRK_res <- STRK_res[abs(STRK_res)>50]
rf_5f_res <- rf_5f_obs-rf_5f_pred
rfsi_5f_res <- rfsi_5f_obs-rfsi_5f_pred
rfsp_5f_res <- rfsp_5f_obs-rfsp_5f_pred

summary(STRK_res)
summary(rf_5f_res)
summary(rfsi_5f_res)
summary(rfsp_5f_res)

STRK_res[which.max(STRK_5f_obs)]
rf_5f_res[which.max(rf_5f_obs)]
rfsi_5f_res[which.max(rfsi_5f_obs)]
rfsp_5f_res[which.max(rfsp_5f_obs)]

STRK_5f_obs[which.min(STRK_res)]
STRK_5f_pred[which.min(STRK_res)]
rf_5f_pred[which.min(STRK_res)]
rfsi_5f_pred[which.min(STRK_res)]
rfsp_5f_pred[which.min(STRK_res)]

summary(STRK_res[which(STRK_5f_obs>100)])
summary(rf_5f_res[which(rf_5f_obs>100)])
summary(rfsi_5f_res[which(rfsi_5f_obs>100)])
summary(rfsp_5f_res[which(rfsp_5f_obs>100)])

summary(STRK_res[which(STRK_5f_obs<10)])
summary(rf_5f_res[which(rf_5f_obs<10)])
summary(rfsi_5f_res[which(rfsi_5f_obs<10)])
summary(rfsp_5f_res[which(rfsp_5f_obs<10)])

summary(STRK_res[which(STRK_5f_obs>10 & STRK_5f_obs<100)])
summary(rf_5f_res[which(rf_5f_obs>10 & rf_5f_obs<100)])
summary(rfsi_5f_res[which(rfsi_5f_obs>10 & rfsi_5f_obs<100)])
summary(rfsp_5f_res[which(rfsp_5f_obs>10 & rfsp_5f_obs<100)])

# scaterplots
my_colors=colorRampPalette(rev(bpy.colors()))(50)
data <- as.data.frame(cbind(STRK_5f_obs, STRK_5f_pred))
p1 <- ggplot(data, aes(x=STRK_5f_pred, y=STRK_5f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 5, 10, 50, 100, 500, 35000, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('5', '10', '50', '100', '500', '35000', '70000+'))+
  theme(plot.title = element_text(hjust = 8),
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
  labs(x = "Residuals [mm]", y = "Observations [mm]", title = "STRK") + coord_fixed()+
  xlim(0, 150) +
  geom_abline(slope=1, intercept=0, size = 0.1)
p1

data <- as.data.frame(cbind(rf_5f_obs, rf_5f_pred))
p2 <- ggplot(data, aes(x=rf_5f_pred, y=rf_5f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 5, 10, 50, 100, 500, 35000, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('5', '10', '50', '100', '500', '35000', '70000+'))+
  theme(plot.title = element_text(hjust = 8),
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
  xlim(0, 150) +
  geom_abline(slope=1, intercept=0, size = 0.1)

data <- as.data.frame(cbind(rfsi_5f_obs, rfsi_5f_pred))
p3 <- ggplot(data, aes(x=rfsi_5f_pred, y=rfsi_5f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 5, 10, 50, 100, 500, 35000, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('5', '10', '50', '100', '500', '35000', '70000+'))+
  theme(plot.title = element_text(hjust = 8),
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
  xlim(0, 150) +
  geom_abline(slope=1, intercept=0, size = 0.1)

data <- as.data.frame(cbind(rfsp_5f_obs, rfsp_5f_pred))
p4 <- ggplot(data, aes(x=rfsp_5f_pred, y=rfsp_5f_obs) ) +
  geom_hex(bins = 50, aes(fill = stat(cut(log(count), breaks = log(c(0, 5, 10, 50, 100, 500, 35000, Inf)), labels = F, right = T, include.lowest = T)))) +
  scale_fill_gradientn(colours = my_colors, name = 'Count', labels = c('5', '10', '50', '100', '500', '35000', '70000+'))+
  theme(plot.title = element_text(hjust = 8),
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
  xlim(0, 150) +
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
# # res_hist2 = histogram(rf_5f_res, main="", xlim=c(-20,20), ylim=c(0,80000), xlab="Residuals [mm]", ylab="", breaks=300, type = "count", col="white",
# #                       scales=(list(cex=0.5)))
# 
# res_hist2 = ggplot(as.data.frame(rf_5f_res), aes(x=rf_5f_res)) +
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
# res_hist3 = ggplot(as.data.frame(rfsi_5f_res), aes(x=rfsi_5f_res)) +
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
# res_hist4 = ggplot(as.data.frame(rfsp_5f_res), aes(x=rfsp_5f_res)) +
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
ggarrange(p1, p2, p3, p4, ncol=4, nrow=1, common.legend = TRUE, legend="bottom") #,
          # ggarrange(res_hist1, res_hist2, res_hist3, res_hist4, ncol=4, nrow=1),
          # ncol=1, nrow=2)
dev.off()

par(mfrow=c(1,1))
