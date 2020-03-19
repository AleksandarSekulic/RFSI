library(spacetime)
library(gstat)
library(rgdal)
library(plotGoogleMaps)
library(hexbin)
library(zoo)
library(snowfall)
library(RSAGA)
library(raster)
library(meteo)
library(doParallel)
library(randomForest)
library(caret)
library(tidyverse)
library(plyr)
library(beepr)
library(hexbin)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
# dir.create("temp_data")

border <- readOGR("catalonia_border/union_of_selected_boundaries_AL4-AL4.shp")

years = 2016:2018
time<-seq(as.Date(paste(as.numeric(min(years)),"-01-01", sep="")), as.Date(paste(max(years),"-12-31", sep="")), by="day")
days<-gsub("-","",time,fixed=TRUE)
daysNum <- length(time)

var="min"
sp.nmax = 100
computeVar = FALSE

load(file = paste("temp_data/", var, "_", min(years), "_", max(years), '.rda', sep=""))
st=stfdf_temp@sp
st@proj4string=CRS("+proj=longlat +datum=WGS84")

local_t= row.names(st)

stfdf_temp.local<-stfdf_temp[local_t,,'tres'] # tres are residuals from global model
stfdf_temp.local@sp@proj4string = CRS("+proj=longlat +datum=WGS84")

data(tvgms)
data(tregcoef)

if(var=="max"){
  coef = as.vector(tregcoef[[8]])
  vario= tvgms[[8]]
} else if (var=="min"){
  coef = as.vector(tregcoef[[5]])
  vario= tvgms[[5]]
} else {
  coef = as.vector(tregcoef[[2]])
  vario= tvgms[[2]]
}

i_1=c(1,1:(daysNum -1))
ip1=c(2,2,3:daysNum)

dir.create(var)

dem <- readGDAL(paste('dem_twi/dem_cat.tif', sep = ""))
r = raster(dem)
names(dem) = 'dem'
dem$twi<- readGDAL(paste('dem_twi/twi_cat.tif', sep = ""))$band1
gg <- as(dem, "SpatialPixelsDataFrame")
gg@proj4string = CRS("+proj=longlat +datum=WGS84")
rm(dem) 

gtt_pred <- tgeom2STFDF(gg, time = time, variable = var)

gc();gc()

tlm  = coef[1] + coef[2]*as.numeric(gtt_pred$temp_geo)+coef[3]*rep(as.numeric(gg$dem),daysNum) + coef[4]*rep(as.numeric(gg$twi),daysNum) 
tlm=matrix(tlm,ncol=daysNum)

### predictions ###
Mpoint=data.frame(x=mean(gg@coords[,1]),y=mean(gg@coords[,2]) )
coordinates(Mpoint)=~x+y
Mpoint@proj4string=CRS("+proj=longlat +datum=WGS84")
st=stfdf_temp@sp

# Find nearest 100 stations
st$dist=spDists(st,Mpoint, longlat =TRUE)
tmp_st<-st[ order(st$'dist') ,]
local_t= row.names(tmp_st[1:sp.nmax,] )

stfdf_temp.local<-stfdf_temp[local_t,,'tres']

gc();gc()

# xxx<- lapply(1:daysNum, function(i) {
registerDoParallel(cores=detectCores()-1)
xxx <- foreach(i = 1:daysNum, .packages = c("raster","spacetime","gstat","rgdal","raster","doParallel","snowfall","meteo")) %dopar% {
  obs=as(stfdf_temp.local[,i_1[i]:ip1[i],'tres'],"STSDF")
  
  if (length(obs@data$tres)<5){
    for (m in 2:10){
      local_t= row.names(tmp_st[1:(sp.max*m),] )
      stfdf_temp.local<-stfdf_temp[local_t,,'tres']
      obs=as(stfdf_temp.local[,i_1[i]:ip1[i],'tres'],"STSDF")
      if (length(obs@data$tres)>5) break
    }
  }
  
  krigeST(as.formula("tres~1"),
          data=obs,
          newdata=STF(as(gg,"SpatialPoints"),
                      stfdf_temp.local@time[i],
                      stfdf_temp.local@endTime[i]),  
          modelList=vario,
          computeVar=FALSE)$var1.pred
  
}
stopImplicitCluster()

res=do.call(cbind,xxx)
temp= tlm + res
row.names(temp)<-1:nrow(gg)
rm(tlm, res)

# temp=as.numeric(temp)
temp=round(temp*10)

for(j in 1:daysNum){
  pre = temp[, j]
  gg@data$pred=pre
  p <- try( raster( as(gg["pred"], "SpatialPixelsDataFrame")) )
  if(inherits(p, "try-error")) {
    p <- rasterize( gg, r, "pred")
  }
  writeRaster(p, paste(var,"/",days[j], sep = ""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
  rm(pre, p)
}
