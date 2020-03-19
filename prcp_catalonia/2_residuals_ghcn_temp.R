library(sp)
library(meteo) # rforge
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

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

border <- readOGR("catalonia_border/union_of_selected_boundaries_AL4-AL4.shp")

dir.create("temp_data")
setwd(paste(wd, "/temp_data/", sep = ""))

var = "max"
el = 'TMAX'
years = 2016:2018

names <- c('staid','date','temp')

data(tvgms)
data(tregcoef)
if(var=="max"){
    coef = as.vector(tregcoef[[8]])
} else if (var=="min"){
    coef = as.vector(tregcoef[[5]])
} else {
    coef = as.vector(tregcoef[[2]])
}

ob <- c()

for (year in years){  
  time<-seq(as.Date(paste(as.numeric(year),"-01-01", sep="")), as.Date(paste(year,"-12-31", sep="")), by="day")
  days<-gsub("-","",time,fixed=TRUE)
  daysNum <- length(time)
  
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

# save(ob, file = paste("ob_", var, ".rda", sep = ""))
load(file = paste("ob_", var, ".rda", sep = ""))

### GHCND stations ###
url <- "ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"
destfile <- paste("ghcnd_stations.txt", sep = "")
download.file(url, destfile)

sta_ghcnd <- read.fwf("ghcnd_stations.txt", width = c(11,9,10,7), header=F) # 11 ILI 12 !!!!!!!!!!!!!!!!!!!!!!!!!!
names(sta_ghcnd) <- c('staid', 'lat', 'lon', 'h')
sta_ghcnd$staid <- as.character(sta_ghcnd$staid)
sta_ghcnd$h <- ifelse(sta_ghcnd$h==-999.9, NA, sta_ghcnd$h)
sta_ghcnd = sta_ghcnd[!is.na(sta_ghcnd$lon),]

# ob <- ob[ob$staid %in% sta$staid, ]

stfdf_temp <- meteo2STFDF ( obs      = ob,
                       stations = sta_ghcnd,
                       crs      = CRS("+proj=longlat +datum=WGS84"),
                       obs.staid.time=c(1,2),
                       stations.staid.lon.lat=c(1,3,2)
)

stfdf_temp = rm.dupl(stfdf_temp, zcol = 1, zero.tol = 0.6) # 1 stations removed

nrowsp <- length(stfdf_temp@sp)

numNA <- apply(matrix(stfdf_temp@data[,"temp"],
                      nrow=nrowsp,byrow=F), MARGIN=1,
               FUN=function(x) sum(is.na(x)))

# Remove stations out of covariates
rem <- numNA != daysNum
stfdf_temp <-  stfdf_temp[rem,drop=F] # 368

# plot(border)
# plot(stfdf_temp@sp, add=T)

### SP OVERLAY ###

r <- raster("../../../meteo/dem_twi/dem.sdat") # DEM of whole world
e <- raster::extract(r,stfdf_temp@sp)
stfdf_temp@sp$dem=e

r <- raster("../../../meteo/dem_twi/twi.sdat") # TWI of whole world  
e <- raster::extract(r, stfdf_temp@sp) 
stfdf_temp@sp$twi=e

stfdf_temp <- stfdf_temp[!is.na(stfdf_temp@sp$dem),] # 1 station removed

time<-seq(as.Date(paste(as.numeric(min(years)),"-01-01", sep="")), as.Date(paste(max(years),"-12-31", sep="")), by="day")
days<-gsub("-","",time,fixed=TRUE)
daysNum <- length(time)

tg <- tgeom2STFDF(stfdf_temp@sp, time = time, variable = var )

identical(tg@sp, stfdf_temp@sp)

stfdf_temp@data$tres  =stfdf_temp@data$temp -  ( coef[1] + coef[2]*as.numeric(tg@data$temp_geo)+coef[3]*rep(as.numeric(stfdf_temp@sp$dem),daysNum) + coef[4]*rep(as.numeric(stfdf_temp@sp$twi),daysNum) )

save(stfdf_temp, file=paste(var, "_", min(years), "_", max(years), '.rda', sep=""))
    
rm(ob, stfdf_temp, tg, e, e1)
    
