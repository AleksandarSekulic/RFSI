# RFSI

This repository contains scripts and data used in the manuscript "Random forest spatial interpolation" by Sekulić A., Kilibarda M., Heuvelink G.B.M., Nikolić M. and Bajat B., that was submitted to the Remote Sensing journal in March 2020.

## Methodology

**Random Forest Spatial Interpolation (RFSI)** is a novel methodology for spatial interpolation using machine learning, i.e. random forest (RF). The main novelty is that it uses observations at n nearest locations and distances from these locations to the prediction location as spatial covariates to improve accuracy of RF.

## Case studies

RFSI is evaluated on two case studies, synthetic and precipitation.

### Synthetic case study

For this case study, multiple realisations of a spatially autocorrelated random field are simulated using a known semivariogram. The performance of RFSi is compared with the performance of ordinary kriging (OK), RFsp (Hengl et al. 2018), Inverse distance Weighting (IDW), nearest neighbour (NN), and trend surface (TS). 

The scripts and data for the synthetic case study are in the *simulation* folder, and are the following:
- *accuracy_results.rda* - R data frame with accuracy metrics from evaluation
- *boxplots.R* - a script for making the boxplots from accuracy metrics (accuracy_results.rda)
- *simulation.R* - a script for simulation, making the models and prediction maps, and calculation of the accuracy metrics for OK, RFSI, RFsp, IDW, NN, and TS

### Precipitation case study (prcp_catalonia)

RFSI is applied to a daily precipitation dataset for Catalonia for the years 2016–2018. Its performance is compared to space-time regression kriging (STRK), RFsp, and regular RF using nested leave-location out 5-fold cross-validation (LLOCV).

The scripts and data for the precipitation case study are in the *prcp_catalonia* folder, and are the following:
- *catalonia_border* - a folder containing Catalonian border (ESRI Shapefile)
- *dem_twi* - a folder containing digital elevation model (DEM) and topographic wetness index (TWI) for Catalonia (GeoTIFF)
- *imerg* - a folder containing IMERG (Integrated Multi-satellitE Retrievals for GPM, Huffman et al. 2014) maps of daily precipitation estimates for 1–4 January 2016 (GeoTIFF)
- *max* - a folder containing maximum daily temperature for 1–4 January 2016 estimated with models proposed by Kilibarda et al. (2014) (GeoTIFF)
- *min* - a folder containing minimum daily temperature for 1–4 January 2016 estimated with models proposed by Kilibarda et al. (2014) (GeoTIFF)
- *models* - a folder containing RF, RFsp, RFSI models, STRK trend (multiple linear regresion) model, STRK sample and fitted semivariogram (Rdata)
- *temp_data* - a folder containing some of the data relevant to script 4 and 5 (Rdata)
- *1_download_gpmdata.R* - a script for downloading the IMERG data
- *2_residuals_ghcn_temp.R* - a script for calculation of the residuals for maximum and minimum daily tempreature
- *3_prediction_temperature.R* - a script for prediction of maximum and minimum daily tempreature
- *4_prcp_case_study_catalonia_RK.R* - a script for making the STRK, RF, RFsp, and RFSI models and LLOCV
- *5_prcp_prediction.R* - a script for making prediction and interquartile range (IQR) maps from STRK, RF, RFsp, and RFSI models
- *stratfolds.R* - a script for creation of spatially stratified folds
*Note that the number in the script name refers to the order in which the scripts should be run. Also, scripts 4 and 5 cam be run with the data in folder temp_data.*

## How to make an RFSI model

TODO

## References

Hengl, T., Nussbaum, M., Wright, M. N., Heuvelink, G. B. M., & Gräler, B. (2018). Random forest as a generic framework for predictive modeling of spatial and spatio-temporal variables. PeerJ, 6, e5518. https://doi.org/10.7717/peerj.5518

Huffman, G., Bolvin, D., Braithwaite, D., Hsu, K., Joyce, R., Xie, P. (2014). Integrated Multi-satellitE Retrievals for GPM (IMERG), version 4.4. NASA's Precipitation Processing Center, accessed 31 March, 2019, ftp://jsimpson.pps.eosdis.nasa.gov/data/imerg/gis/

Kilibarda, M., Hengl, T., Heuvelink, G. B. M., Gräler, B., Pebesma, E., Perčec Tadic, M.,  Bajat, B. (2014). Spatio-temporal interpolation of daily temperatures for global land areas at 1 km resolution. Journal of Geophysical Research Atmospheres, 119(5), 2294–2313. https://doi.org/10.1002/2013JD020803