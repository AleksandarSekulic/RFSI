library(R.utils)
library(doParallel)
#library(RCurl)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
dir.create("imerg")
setwd("imerg")
years <- 2016:2018

# ftp://jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/2019/02/
# 3B-DAY-L.GIS.IMERG.20190201.V05B.zip
  
time=seq(as.Date(paste(years[1],"-01-01", sep="")), as.Date(paste(years[length(years)],"-12-31", sep="")), by="day")
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

version <- c('06B', '05B', '04B')

# ftp://arthurhou.pps.eosdis.nasa.gov/gpmdata/2016/01/01/
# gis/3B-DAY-GIS.MS.MRG.3IMERG.20160101-S000000-E235959.0000.V06A.zip
# ftp://arthurhou.pps.eosdis.nasa.gov/gpmallversions/V05/2016/01/01/gis/
# 3B-DAY-GIS.MS.MRG.3IMERG.20160101-S000000-E235959.0000.V05B.zip

registerDoParallel(cores=15)
foreach(i = 1:daysNum, .packages = c('R.utils')) %dopar% {
# for(i in c(1:daysNum)) {
    url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/'
              ,substr(days[i], 1, 4), '/', substr(days[i], 5, 6), '/', 
              '3B-HHR-L.MS.MRG.3IMERG.', days[i], '-S233000-E235959.1410.V06A.1day.zip', sep = '')
    destfile = paste(days[i], '.zip', sep = '')
    d <- try(download.file(url, destfile))
    if(inherits(d, "try-error")) {
      cat(day)
      for (v in 1:length(version)) {
        url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/'
                    ,substr(days[i], 1, 4), '/', substr(days[i], 5, 6), '/',
                    '3B-HHR-L.MS.MRG.3IMERG.', days[i], '-S233000-E235959.1410.V', version[v], '.1day.zip', sep = '')
        d <- try(download.file(url, destfile))
        if(inherits(d, "try-error")) {
          next
        } else {
          break
        }
      }
    }
    unzip(destfile)
    unlink(destfile)
}
stopImplicitCluster()

files_remove <- list.files(pattern = "1day.liquid")
unlink(files_remove)
files_remove <- list.files(pattern = "1day.ice")
unlink(files_remove)

### 2017-10 missing data ###
year <- 2017
time=seq(as.Date(paste(year,"-10-01", sep="")), as.Date(paste(year,"-10-31", sep="")), by="day")
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

registerDoParallel(cores=15)
foreach(i = 1:daysNum, .packages = c('R.utils')) %dopar% {
  url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/'
              ,year, '/', substr(days[i], 5, 6), '/V04/', 
              '3B-DAY-L.GIS.IMERG.', days[i], '.V04B.zip', sep = '')
  destfile = paste(days[i], '.zip', sep = '')
  download.file(url, destfile)
  unzip(destfile)
  unlink(destfile)
}
stopImplicitCluster()

### + october V05B ###
year=2017

time=seq(as.Date(paste(year,"-10-01", sep="")), as.Date(paste(year,"-10-31", sep="")), by="day")
daysNum = length(time)
days=gsub("-","",time,fixed=TRUE)

registerDoParallel(cores=15)
foreach(i = 1:daysNum, .packages = c('R.utils')) %dopar% {
  url = paste('ftp://aleksandarsale.sekulic%40gmail.com:aleksandarsale.sekulic%40gmail.com@jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/'
              ,year, '/', substr(days[i], 5, 6), '/V04/', 
              '3B-HHR-L.MS.MRG.3IMERG.', days[i], '-S053000-E055959.0330.V04B.1day.zip', sep = '')
  destfile = paste(days[i], '.zip', sep = '')
  download.file(url, destfile)
  unzip(destfile)
  unlink(destfile)
}
stopImplicitCluster()