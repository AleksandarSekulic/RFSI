# RFSI

This repository contains scripts and data used in the manuscript:

Sekulić, A., Kilibarda, M., Heuvelink, G. B., Nikolić, M. & Bajat, B. Random Forest Spatial Interpolation.Remote. Sens. 12, 1687, https://doi.org/10.3390/rs12101687 (2020).

## Methodology

**Random Forest Spatial Interpolation (RFSI)** is a novel methodology for spatial interpolation using machine learning, i.e. random forest (RF) (Breiman 2001). The main novelty is that it uses observations at n nearest locations and distances from these locations to the prediction location as spatial covariates to improve accuracy of RF.

***Note that Out-of-bag (OOB) error statistics from RFSI model are biased and should not be considered as accuracy metrics (they do not show spatial accuracy)! The proper way to assess accuaracy of the RFSI model is by using the nested k-fold cross-validation (`cv.rfsi` function of the R package [meteo](https://github.com/AleksandarSekulic/Rmeteo), Sekulić et al. 2020b).***

## Case studies

### Synthetic case study

For this case study, multiple realisations of a spatially autocorrelated random field are simulated using a known semivariogram. The performance of RFSI is compared with the performance of ordinary kriging (OK), RFsp (Hengl et al. 2018), Inverse distance Weighting (IDW), nearest neighbour (NN), and trend surface (TS). 

The scripts and data for the synthetic case study are in the [simulation](simulation/) folder, and are the following:
- [accuracy_results.rda](simulation/accuracy_results.rda) - R data frame with accuracy metrics from evaluation,
- [boxplots.R](simulation/boxplots.R) - a script for making the boxplots from accuracy metrics (accuracy_results.rda),
- [sim_500.R](simulation/sim_500.R) - a script for number of nearest locations sensitivity study,
- [simulation.R](simulation/simulation.R) - a script for simulation, making the models and prediction maps, and calculation of the accuracy metrics for OK, RFSI, RFsp, IDW, NN, and TS.

### Precipitation case study

RFSI is applied to a daily precipitation dataset for Catalonia for the years 2016–2018. Its performance is compared to space-time regression kriging (STRK), IDW, RFsp, and regular RF using nested leave-location out 5-fold cross-validation (LLOCV).

The scripts and data for the precipitation case study are in the [prcp_catalonia](prcp_catalonia/) folder, and are the following:
- [catalonia_border](prcp_catalonia/catalonia_border/) - a folder containing Catalonian border (ESRI Shapefile),
- [dem_twi](prcp_catalonia/dem_twi/) - a folder containing digital elevation model (DEM) and topographic wetness index (TWI) for Catalonia (GeoTIFF),
- [imerg](prcp_catalonia/imerg/) - a folder containing IMERG (Integrated Multi-satellitE Retrievals for GPM, Huffman et al. 2014) maps of daily precipitation estimates for 1–4 January 2016 (GeoTIFF),
- [max](prcp_catalonia/max/) - a folder containing maximum daily temperature for 1–4 January 2016 estimated with models proposed by Kilibarda et al. (2014) (GeoTIFF),
- [min](prcp_catalonia/min/) - a folder containing minimum daily temperature for 1–4 January 2016 estimated with models proposed by Kilibarda et al. (2014) (GeoTIFF),
- [models](prcp_catalonia/models/) - a folder containing RF, RFsp, RFSI models, STRK trend (multiple linear regresion) model, STRK sample and fitted semivariogram (Rdata),
- [temp_data](prcp_catalonia/temp_data/) - a folder containing some of the data relevant to script 4 and 5 (Rdata),
- [1_download_gpmdata.R](prcp_catalonia/1_download_gpmdata.R) - a script for downloading the IMERG data,
- [2_residuals_ghcn_temp.R](prcp_catalonia/2_residuals_ghcn_temp.R) - a script for calculation of the residuals for maximum and minimum daily tempreature,
- [3_prediction_temperature.R](prcp_catalonia/3_prediction_temperature.R) - a script for prediction of maximum and minimum daily tempreature,
- [4_prcp_case_study_catalonia_RK.R](prcp_catalonia/4_prcp_case_study_catalonia_RK.R) - a script for making the STRK, IDW, RF, RFsp, and RFSI models and LLOCV,
- [5_prcp_prediction.R](prcp_catalonia/5_prcp_prediction.R) - a script for making prediction and interquartile range (IQR) maps from STRK, IDW, RF, RFsp, and RFSI models,
- [stratfolds.R](prcp_catalonia/stratfolds.R) - a script for creation of spatially stratified folds.

*Note that the number in the script name refers to the order in which the scripts should be run. Also, scripts 4 and 5 can be run with the data in folder temp_data.*

### Temperature case study

RFSI is applied to a daily mean temperature dataset for Croatia for the year 2018. Its performance is compared to space-time regression kriging (STRK), IDW, RFsp, and regular RF using nested leave-location out 10-fold cross-validation (LLOCV).

The scripts and data for the precipitation case study are in the [temp_croatia](temp_croatia/) folder, and are the following:
- [borders](temp_croatia/borders/) - a folder containing Croatian border (ESRI Shapefile),
- [models](temp_croatia/models/) - a folder containing RF, RFsp (over 100 MB), RFSI models (Rdata),
- [temp_data](temp_croatia/temp_data/) - a folder containing some of the data relevant to script 1 and 2 (Rdata),
- [1_models.R](temp_croatia/1_models.R) - a script for making the IDW, RF, RFsp, and RFSI models and LLOCV,
- [2_predictions.R](temp_croatia/2_predictions.R) - a script for making prediction maps from IDW, RF, RFsp, and RFSI models,
- [dem.tif](temp_croatia/dem.tif) - digital elevation model (DEM) for Croatia (GeoTIFF),
- [stratfolds.R](temp_croatia/stratfolds.R) - a script for creation of spatially stratified folds.

## How to make an RFSI model

Complete RFSI examples (including tune.rfsi and cv.rfsi) can be found in the R package [meteo](https://github.com/AleksandarSekulic/Rmeteo), in the [demo](https://github.com/AleksandarSekulic/Rmeteo/tree/master/demo) folder.
```
library(meteo)
library(sp)
library(spacetime)
library(gstat)
library(plyr)
library(xts)
library(snowfall)
library(doParallel)
library(CAST)
library(ranger)

# preparing data
demo(meuse, echo=FALSE)
meuse <- meuse[complete.cases(meuse@data),]

#################### rfsi ####################

fm.RFSI <- as.formula("zinc ~ dist + soil + ffreq")

rfsi_model <- rfsi(formula = fm.RFSI,
                   data = meuse,
                   zero.tol = 0,
                   n.obs = 5, # number of nearest observations
                   s.crs = NA, # or meuse@proj4string # nedded only if in lon/lat (WGS84)
                   t.crs = NA, # or meuse@proj4string # nedded only if in lon/lat (WGS84)
                   cpus = detectCores()-1,
                   progress = TRUE,
                   # ranger parameters
                   importance = "impurity",
                   seed = 42,
                   num.trees = 250,
                   mtry = 5,
                   splitrule = "variance",
                   min.node.size = 5,
                   sample.fraction = 0.95,
                   quantreg = FALSE)
rfsi_model

# Note that OOB error statistics are biased and should not be considered as accuracy metrics
# (they do not show spatial accuracy)!
# The proper way to assess accuaracy of the RFSI model is by using the nested k-fold
# cross-validation (cv.rfsi function)

sort(rfsi_model$variable.importance)

#################### pred.rfsi ####################

rfsi_prediction <- pred.rfsi(model = rfsi_model,
                             data = meuse, # meuse.df (use data.staid.x.y.time)
                             zcol = "zinc",
                             newdata = meuse.grid, # meuse.grid.df (use newdata.staid.x.y.time)
                             output.format = "SpatialPixelsDataFrame",
                             zero.tol = 0,
                             s.crs = meuse@proj4string, # NA
                             newdata.s.crs = meuse@proj4string, # NA
                             t.crs = meuse@proj4string, # NA
                             cpus = detectCores()-1,
                             progress = TRUE
)

spplot(rfsi_prediction['pred'])
```

## References

Breiman, L. (2001). Random Forests. Machine Learning, 45(1), 5–32. https://doi.org/10.1023/A:1010933404324

Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B. M., & Gräler, B. (2018). Random forest as a generic framework for predictive modeling of spatial and spatio-temporal variables. PeerJ, 6, e5518. https://doi.org/10.7717/peerj.5518

Huffman, G., Bolvin, D., Braithwaite, D., Hsu, K., Joyce, R., Xie, P. (2014). Integrated Multi-satellitE Retrievals for GPM (IMERG), version 4.4. NASA's Precipitation Processing Center, accessed 31 March, 2019, ftp://jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/

Kilibarda, M., Hengl, T., Heuvelink, G. B. M., Gräler, B., Pebesma, E., Perčec Tadic, M.,  Bajat, B. (2014). Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution. Journal of Geophysical Research Atmospheres, 119(5), 2294–2313. https://doi.org/10.1002/2013JD020803