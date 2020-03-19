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
library(zoo)
library(ggpubr)

wd="/media/geocomp/060C0FE30C0FCC9D/RFnl/final_catalonia/"
setwd(wd)

# lonmin=3.3 ;lonmax=7.3 ; latmin=50.7 ;latmax=53.6 # Holland
# nld_border <- getData("GADM", country="NLD", level=0)
# lonmin=0 ;lonmax=4 ; latmin=40 ;latmax=43 # Ctalonia
border <- readOGR("catalonia_border/union_of_selected_boundaries_AL4-AL4.shp")

# dir.create("temp_data")
setwd(paste(wd, "temp_data/", sep = ""))

v = "prcp"
# /media/sekulic/e42d0594-47d8-4b72-99f7-bf3809dbeb4f/__DOKTORAT/MEDCLIVAR 2018/stations/1all_data.R
# years = 2009:2018
# years = 2016:2018
# time1<-seq(as.Date(paste(as.numeric(min(years)),"-01-01", sep="")), as.Date(paste(max(years),"-12-31", sep="")), by="day")
# days<-gsub("-","",time,fixed=TRUE)
# daysNum <- length(time)

el = 'PRCP'
names <- c('staid','date','prcp')

# load(paste(v, "_", min(years), "_", max(years), '.rda', sep=""))
load('stfdf_temp.rda')
time=zoo::index(stfdf@time)
days<-gsub("-","",time,fixed=TRUE)
daysNum <- length(time)

st=1

load('PRK_obs.rda')
load('PRK_pred.rda')
load('rf_5f_pred.rda')
load('rfsi_5f_pred.rda')
load('rfsp_5f_pred.rda')
PRK_pred <- PRK_pred[!is.na(PRK_obs)]

temp_df <- as.data.frame(stfdf)
temp_df <- temp_df[complete.cases(temp_df), ]
temp_df <- cbind(temp_df[, c(3,4,9,11)], PRK_pred, rf_5f_pred, rfsi_5f_pred, rfsp_5f_pred)

### plot predictions and observations at 2 stations ###

year = "2018"

gt = F
# promeni geom_line !!!

for (year in c("2016", "2017", "2018")){

a1 <- temp_df[as.numeric(as.character(temp_df$dem))==826, c("prcp", "time")]
n2 <- as.character(stfdf@sp[stfdf@sp$dem==826, ]@data$staid)
b1 <- temp_df[as.numeric(as.character(temp_df$dem))==826, c("PRK_pred", "time")]
c1 <- temp_df[as.numeric(as.character(temp_df$dem))==826, c("rf_5f_pred", "time")]
d1 <- temp_df[as.numeric(as.character(temp_df$dem))==826, c("rfsi_5f_pred", "time")]
e1 <- temp_df[as.numeric(as.character(temp_df$dem))==826, c("rfsp_5f_pred", "time")]
time=zoo::index(stfdf@time)

# summary(a-b)
# summary(a1-b1)

a1 = a1[format(as.Date(time, format="%d/%m/%Y"),"%Y")==year, "prcp"]
b1 = b1[format(as.Date(time, format="%d/%m/%Y"),"%Y")==year, "PRK_pred"]
c1 = c1[format(as.Date(time, format="%d/%m/%Y"),"%Y")==year, "rf_5f_pred"]
d1 = d1[format(as.Date(time, format="%d/%m/%Y"),"%Y")==year, "rfsi_5f_pred"]
e1 = e1[format(as.Date(time, format="%d/%m/%Y"),"%Y")==year, "rfsp_5f_pred"]
time=time[format(as.Date(time, format="%d/%m/%Y"),"%Y")==year]

if (gt) {
  b1 = b1[a1>=1]
  c1 = c1[a1>=1]
  d1 = d1[a1>=1]
  e1 = e1[a1>=1]
  time=time[a1>=1]
  a1 = a1[a1>=1]
  time = as.factor(as.numeric(strftime(time, format = "%j")))
} else {
  b1 = b1[a1<1]
  c1 = c1[a1<1]
  d1 = d1[a1<1]
  e1 = e1[a1<1]
  time=time[a1<1]
  a1 = a1[a1<1]
  time = as.Date(time)
}

# b1 = b1[a1==0]
# c1 = c1[a1==0]
# d1 = d1[a1==0]
# time=time[a1==0]
# a1 = a1[a1==0]
# time = as.Date(time)

# Library
library(tidyverse)

# Create data
data=data.frame(x=time, prcp=a1, pooled_rk=b1, rf=c1, rfsi=d1, rfsp=e1 )
data$max <- apply(data[, 2:6], 1, function(x) as.numeric(max(x)))
data$min <- apply(data[, 2:6], 1, function(x) as.numeric(min(x)))
summary(data)

Sys.setlocale("LC_TIME", "C")
p1 <- ggplot(data) +
  geom_line( aes(x=x, y=pooled_rk, color="green", group = 1), alpha = 0.5) +
  geom_line( aes(x=x, y=rf, color="mediumblue", group = 1), alpha = 0.5) +
  geom_line( aes(x=x, y=rfsi, color="red", group = 1), alpha = 0.5) +
  geom_line( aes(x=x, y=rfsp, color="orange", group = 1), alpha = 0.5) +
  # geom_segment( aes(x=x, xend=x, y=0, yend=max), color="grey", size=0.3, alpha=0.7) +
  
  geom_line( aes(x=x, y=prcp, color="black", group = 1), alpha = 0.5) +
  
  # geom_point( aes(x=x, y=pooled_ok, fill="green"), size=1.5, alpha=0.5, shape=21, color="green" ) +
  # geom_point( aes(x=x, y=rf_nl, fill="mediumblue"), size=1.5, alpha=0.5, shape=21, color="mediumblue" ) +
  # geom_point( aes(x=x, y=rf_avg, fill="orange"), size=1.5, alpha=0.5, shape=21, color="orange" ) +
  # geom_point( aes(x=x, y=prcp, fill="darkorange", color="black"), size=1.5, alpha=0.5, shape=21) +
  
  # geom_line( aes(x=x, y=prcp, color="black", group = 1)) +
  
  # coord_flip()+
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
    # legend.position = "none",
    # panel.border = element_blank()
  ) +
  ggtitle(paste("Year ", year)) +
  
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 7),
        text = element_text(size = 7),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size= unit(0.2, "cm"),
        # legend.margin = unit(0, "cm"),
        legend.title = element_text(size=7, face="bold")) + 
  # scale_fill_manual(values=c("blue", "cyan4")) + 
  # scale_colour_manual(values=c("white", "black"))
  # scale_fill_manual(name = "Precipitation",
  #                     labels = c("Observations", "Pooled OK", "RFnl", "RFavg"),
  #                     values = c("darkorange"="darkorange", "green"="green", "mediumblue"="mediumblue", "red"="red")) +
  # # scale_shape_manual(name = "Precipitation",
  # #                    labels = c("Observations", "Pooled OK", "RFnl", "RFavg"),
  # #                    values = c(21, 21, 21, 21)) +
  # scale_colour_manual(name = "Precipitation",
  #                     labels = c("Observations", "Pooled OK", "RFnl", "RFavg"),
  #                     values = c("black"="black", "green"="green", "mediumblue"="mediumblue", "red"="red"))
  # scale_fill_manual(name = "Predictions",
  #                 labels = c("Pooled OK", "RFnl", "RFavg"),
  #                 values = c("black"="black", "green"="green", "mediumblue"="mediumblue", "orange"="orange"),
  #                 limits = c("green", "mediumblue", "orange")) +
  # scale_shape_manual(name = "Precipitation",
  #                    labels = c("Observations", "Pooled OK", "RFnl", "RFavg"),
  #                    values = c(21, 21, 21, 21)) +
  scale_colour_manual(name = "Precipitation",
                      labels = c("Observations", "Pooled RK", "RF", "RFSI", "RFsp"),
                      values = c("black"="black",
                                 "green"="green",
                                 "mediumblue"="mediumblue",
                                 "red"="red",
                                 "orange"="orange"),
                      limits = c("black", "green", "mediumblue", "red", "orange"))
print(p1)

if (gt) {
  p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

if (year == "2017"){
  p1 = p1 + ylab("Precipitation [mm]")
} else {
  p1 = p1 + ylab("")
}

if (year == "2018"){
  if (gt){
    p1 = p1 + xlab("Day of the year")
  } else {
    p1 = p1 + xlab("Month")
  }
} else {
  p1 = p1 + xlab("")
}

print(p1)

if (gt) {
  a="gt"
} else {
  a="lt"
}

if (year == "2016"){
  p2016 = p1
} else if (year == "2017"){
  p2017 = p1
} else {
  p2018 = p1
}

}

# tiff(paste("../plot/ts_", a, "1.tiff", sep=""), width = 174, height = 210, units = 'mm', res = 600, compression = "lzw")
jpeg(paste("../plot/ts_", a, "1.jpeg", sep=""), width = 174, height = 210, units = 'mm', res = 600)
par(mfrow=c(1,1))
ggarrange(p2016, p2017, p2018, ncol=1, nrow=3, common.legend = TRUE, legend="bottom")
dev.off()

par(mfrow=c(1,1))

# time = as.Date(time)
# Sys.setlocale("LC_TIME", "C")
# plot(d1 ~ time , type="b" , bty="l" , xlab="month" , ylab="Precipitation [mm]" , col=alpha("blue", 0.5) , lwd=0.5 , pch=20)# , ylim=c(1,5) )
# lines(a1 ~ time , col= alpha("orange", 0.5), lwd=0.5 , pch=20 , type="b" )
# lines(c1 ~ time , col=alpha("red", 0.5) , lwd=0.5 , pch=20 , type="b" )
# lines(b1 ~ time , col=alpha("chartreuse3", 0.5) , lwd=0.5 , pch=20 , type="b" )





Sys.setlocale("LC_TIME", "C")
# tiff("plot/pred_per_station.tiff", width = 25, height = 25, units = 'cm', res = 300)
# jpeg("plot/pred_per_station.jpeg", width = 25, height = 25, units = 'cm', res = 300)
par(mfrow=c(1,1))
# plot.zoo(cbind(a1[1:365], b1[1:365], c1[1:365], d1[1:365]),
plot.zoo(cbind(a1, b1, c1, d1),
         plot.type = "single",
         col = c("orange", "chartreuse3", "red", "blue"),
         ylab = "Predictions [mm]", xlab = "Month", main = paste(n2, " (H = 145 m)"),
         xlim = NULL)#c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# abline(h=0, lty=2)
legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Observation", "Pooled OK", "RFnl", "RFavg"),  lty = 1, bty = "n", col = c("black", "chartreuse3", "red", "blue"), cex = 1)

plot.zoo(cbind(a, b, c, d),
         plot.type = "single",
         col = c("black", "chartreuse3", "red", "blue"),
         ylab = "Predictions [mm]", xlab = "", main = paste(n1, " (H = 2523 m)"),
         xlim = NULL)#c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
abline(h=0, lty=2)
legend("topright", inset=c(0,0), y.intersp = 1, legend = c("Observation", "Pooled OK", "RFnl", "RFavg"),  lty = 1, bty = "n", col = c("black", "chartreuse3", "red", "blue"), cex = 1)
# dev.off()
par(mfrow=c(1,1))










