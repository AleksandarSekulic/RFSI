# ### install latest version of R meteo package ###
# install.packages("meteo", repos="http://R-Forge.R-project.org")
# 
# # build from source
# library(devtools)
# build(pkg = "/home/sekulic/meteo/packageWV/pkg/")
# install(pkg = "/home/sekulic/meteo/packageWV/pkg/")

### load packages ###
library(raster)
library(gstat)
library(doParallel)
library(ranger)
library(RColorBrewer)
library(ggplot2)
library(rgeos)
library(nabor)
library(spatial)
library(fields)
library(splines)
library(GSIF)
library(DescTools)
library(dplyr)
library(ggpubr)
library(grid)
library(meteo)

### Number of nearest observations (stations) ###
n_st <- 25 # 10
#################################################

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

accuracy_results <- c()
# load(file = "accuracy_results.rda") # simulation results
# summary(accuracy_results)

### Create BBOX ###
xmin=0; xmax=500; ymin=0; ymax=500
bbox = extent(xmin, xmax, ymin, ymax)
bbox <- as(bbox, "SpatialPolygons")
fishnet <- raster(bbox)
res(fishnet) <- c(1, 1) # 1 x 1

st <- as(fishnet, "SpatialPoints")

### seeds ###
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(42)
seeds <- round(runif(100, 1, 1000))

### check only first 5 seeds ###
# seeds <- seeds[1:5]
################################

cpus <- detectCores()-1 # more thatn 3 will full the RAM memory

### TO CHECK ONLY ONE NUGGET/RANGE/NUMBER OF TRAINING POINTS/SEED ###
### NUGGET ###
nug <- 2.5
##############

### RANGE ###
# range <- 50
##############
ranges <- c(25,50,75,100,200)# c(50, 100, 150, 200)

### N POINTS ###
n_po <- 500
##############

### SEED ###
seed <- seeds[1]
##############


        
seed_result <- c()

### randomly selection of training points ###
set.seed(seed)
fold_i <- sample.int(n = length(st), size = n_po)
head(fold_i)
fold <- rep(FALSE, length(st))
fold[fold_i] <- TRUE
# summary(fold)
values(fishnet) <- fold
# plot(fishnet)

##########################################################################################
### Simulation ###########################################################################
##########################################################################################

for (range in ranges) {
  ### variogram ###
  vmf <- vgm((sill-nug), "Sph", range, nug)
  
  ### simulation ###
  set.seed(seed)
  sim_ok <- krige(z ~ 1, loc = NULL, newdata = st,
                  model = vmf, nsim = 1, nmax = 50, # nmax - number of points
                  beta = 20, dummy = T)
  # summary(sim_ok)
  values(fishnet) <- sim_ok$sim1
  # plot(fishnet)
  sim_ok$fold <- fold
  
  test <- sim_ok[!sim_ok$fold, ] # test points
  train <- sim_ok[sim_ok$fold, ] # training points
  
  ##########################################################################################
  ### Ordinary Kriging validation ##########################################################
  ##########################################################################################
  
  ok_plot <- krige(sim1 ~ 1, train, as(fishnet, "SpatialPixelsDataFrame"), model = vmf, nmax = n_st) #, beta=0) - beta is for SK
  # ok_plot <- krige(sim1 ~ 1, train, as(fishnet, "SpatialPixelsDataFrame"), model = auto_vmf$var_model, nmax = n_st) #, beta=0)
  ok_raster <- fishnet
  values(ok_raster) <- ok_plot$var1.pred

  ##########################################################################################
  ### RFSI #################################################################################
  ##########################################################################################
  
  dir.create(paste(n_po, "_", range, sep=""))
  
  ### observations at nearest stations and deistances to them - calculation ###
  st_num <- length(sim_ok)
  sim_df <- as.data.frame(sim_ok)
  
  nos <- 1:40# c(2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50)
  n_obs = max(nos)
  
  nearest_obs <- near.obs(
    locations = sim_ok,
    observations = train,
    zcol = "sim1",
    n.obs = n_obs
  )
  sim_df <- cbind(sim_df, nearest_obs)
  sim_df.dev <- sim_df[sim_df$fold,]
  
  rmse_n <- c()
  
  for (n in nos) {
    fm.obs = paste(names(nearest_obs)[c(1:n, (n_obs+1):(n_obs+n))], collapse="+")
    fm.RFSI <- as.formula(paste("sim1 ~ ", fm.obs))
    print(fm.RFSI)
    
    ### RFSI model ###
    set.seed(seed)
    rf_RFSI <- ranger(fm.RFSI, data = sim_df.dev, importance = "impurity", num.trees = 250) #, mtry = n_st) # mtry - kriging range
    # rf_RFSI
    
    ### validation ###
    rfsi_plot <- predict(rf_RFSI, sim_df) # doesn't take training stations for their predictions!
    rfsi_raster <- fishnet
    values(rfsi_raster) <- rfsi_plot$predictions
    
    pred_test_rfsi <- raster::extract(rfsi_raster, test)
    rmse_n = c(rmse_n, sqrt(mean((test$sim1 - pred_test_rfsi)^2, na.rm=TRUE)))
    
    rfsi_raster <- rfsi_raster * 100
    writeRaster(rfsi_raster, paste(n_po, "_", range, "/", n, ".tif", sep=""), "GTiff",NAflag= -32767, datatype='INT2S', overwrite=T)
    # plot(rfsi_raster)
    
  }
  
  assign(paste("rmse_", n_po, "_", range, sep=""), rmse_n)
  
}


### Maps with different n ###
r <- raster("500_50/2.tif")
r <- r/100
all_plot <- as(r, "SpatialPixelsDataFrame")
# names(all_plot@data)[1] <- "n = 2"

for (n in c(5,10,20,30,40)) {
  r <- raster(paste("500_50/", n, ".tif", sep=""))
  r <- r/100
  # all_plot@data <- cbind(all_plot@data, values(r))
  all_plot[[paste("X", n, sep="")]] <- values(r)
}

# tiff("plot/n_comparison.tiff", width = 140, height = 93, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/n_comparison.jpeg", width = 140, height = 93, units = 'mm', res = 600) # res = 1200
print(spplot(all_plot, col.regions=terrain.colors, par.settings=list(fontsize=list(text=7)),
             # sp.layout = list(pts, which=1, col="black", cex=0.1, alpha=0.5),
             # scales = list(draw = TRUE),
             as.table=TRUE,
             names.attr = paste("n = ", c(2,5,10,20,30,40), sep="")
)) #, main = "Comparison of predictions"
dev.off()

### RMSE ~ n ###

plot(rmse_500_25~nos)
plot(rmse_500_50~nos)
plot(rmse_500_75~nos)
plot(rmse_500_100~nos)
plot(rmse_500_200~nos)

# plot(rmse_500_50~nos)
# plot(rmse_500_100~nos)
# plot(rmse_500_150~nos)
# plot(rmse_500_200~nos)
# sve na jedan plot i obeleze najmanje rmse
# rmse_500_50 <- round(rmse_500_50, 2)
# rmse_500_100 <- round(rmse_500_100, 2)
# rmse_500_150 <- round(rmse_500_150, 2)
# rmse_500_200 <- round(rmse_500_200, 2)

theme = theme_set(theme_minimal())
rmse_n_plot <-
  ggplot() + #data = my_sum, aes(x=Npoints, y=mean, group=c(Method))) +
    geom_line(aes(x=nos, y=rmse_500_25, color="a"), size=0.5) +
    geom_point(aes(x=nos, y=rmse_500_25, color="a"), size=0.8) +
    geom_point(aes(x=nos[which.min(rmse_500_25)], y=min(rmse_500_25), color="a"), size=2.5) +
    
    geom_line(aes(x=nos, y=rmse_500_50, color="b"), size=0.5) +
    geom_point(aes(x=nos, y=rmse_500_50, color="b"), size=0.8) +
    geom_point(aes(x=nos[which.min(rmse_500_50)], y=min(rmse_500_50), color="b"), size=2.5) +
    
    geom_line(aes(x=nos, y=rmse_500_75, color="c"), size=0.5) +
    geom_point(aes(x=nos, y=rmse_500_75, color="c"), size=0.8) +
    geom_point(aes(x=nos[which.min(rmse_500_75)], y=min(rmse_500_75), color="c"), size=2.5) +
    
    geom_line(aes(x=nos, y=rmse_500_100, color="d"), size=0.5) +
    geom_point(aes(x=nos, y=rmse_500_100, color="d"), size=0.8) +
    geom_point(aes(x=nos[which.min(rmse_500_100)], y=min(rmse_500_100), color="d"), size=2.5) +
    
    geom_line(aes(x=nos, y=rmse_500_200, color="e"), size=0.5) +
    geom_point(aes(x=nos, y=rmse_500_200, color="e"), size=0.8) +
    geom_point(aes(x=nos[which.min(rmse_500_200)], y=min(rmse_500_200), color="e"), size=2.5) +
  # scale_fill_manual(breaks = c("50", "100", "150", "200"),
    #                   values=c("red", "green", "deepskyblue", "magenta")) +
    # scale_color_manual(breaks = c("50", "100", "150", "200"),
    #                    values=c("red", "green", "deepskyblue", "magenta")) +
    # scale_fill_manual(name = 'Ranges', 
    #                   values =c("a"="red", "b"="green", "c"="deepskyblue", "d"="blueviolet", "e"="chocolate1"),
    #                   labels = c("25", "50", "75", "100", "200")) +
    scale_colour_manual(name = 'Ranges', 
                        values =c("a"="red", "b"="green", "c"="deepskyblue", "d"="blueviolet", "e"="chocolate1"),
                        labels = c("25", "50", "75", "100", "200")) +
    xlab("Number of nearest locations") +
    ylab("RMSE") +
    theme(plot.title = element_text(hjust = 8),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          text = element_text(size = 8),
          legend.key.size= unit(0.2, "cm"),
          legend.margin = unit(0, "cm"),
          legend.title = element_text(size=8, face="bold"),
          legend.text=element_text(size=8))

rmse_n_plot

# tiff("plot/rmse_n.tiff", width = 174, height = 50, units = 'mm', res = 600, compression = "lzw")
jpeg("plot/rmse_n.jpeg", width = 174, height = 50, units = 'mm', res = 600) # res = 1200
rmse_n_plot
dev.off()


