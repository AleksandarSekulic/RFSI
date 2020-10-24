### install latest version of R meteo package ###
# install.packages("meteo", repos="http://R-Forge.R-project.org")

# build from source
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

nugs <- c(0, 2.5, 5)
sill <- 10
ranges <- c(50, 200)

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

### number of simulated points ##
n_pois <- c(100, 200, 500, 1000, 2000, 5000)

### seeds ###
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(42)
seeds <- round(runif(100, 1, 1000))

### check only first 5 seeds ###
# seeds <- seeds[1:5]
################################

# if (n_po < 2000) {
#   cpus <- 10
# }
# if (n_po >= 2000) {
#   cpus <- 3
# }
# cpus <- 10
cpus <- detectCores()-1 # more thatn 3 will full the RAM memory

#############
### DO NOT GO THORUGH FOR LOOPS IF YOU WANT TO CHECK ONLY ONE NUGGET/RANGE/NUMBER OF TRAINING POINTS/SEED!!! ###
#############

### TO CHECK ONLY ONE NUGGET/RANGE/NUMBER OF TRAINING POINTS/SEED ###
# ### NUGGET ###
# nug <- 2.5
# ##############
# 
# ### RANGE ###
# range <- 50
# ##############
#
# ### N POINTS ###
# n_po <- 500
# ##############
#
# ### SEED ###
# seed <- seeds[1]
# ##############

### loop through nuggets ###
for (nug in nugs) {
  cat("NUGGET")
  cat(nug)
  cat("\n")
  ### loop through ranges ###
  for (range in ranges) {
    cat("RANGE")
    cat(range)
    cat("\n")

### loop through different number of training points ###
    for(n in 1:length(n_pois)) {
      n_po <- n_pois[n]
      
      cat(n_po)
      cat("\n")
    
      # loop through all seeds (parallel) ###
      registerDoParallel(cores=cpus)
      seed_results <- foreach (see = 1:length(seeds), .packages = c("raster","gstat","ranger","nabor","GSIF","DescTools")) %dopar% {
      # for(see in 1:length(seeds)){
        seed <- seeds[see]
        
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
    
        # auto_vmf <- autofitVariogram(sim1 ~ 1, train, model = "Sph")
        # auto_vmf$var_model
        time1 <- Sys.time()
        ok_plot <- krige(sim1 ~ 1, train, as(fishnet, "SpatialPixelsDataFrame"), model = vmf, nmax = n_st) #, beta=0) - beta is for SK
        # ok_plot <- krige(sim1 ~ 1, train, as(fishnet, "SpatialPixelsDataFrame"), model = auto_vmf$var_model, nmax = n_st) #, beta=0)
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "PT", nug, range, n_po, seed, "OK"))
        ok_raster <- fishnet
        values(ok_raster) <- ok_plot$var1.pred
        # plot(ok_raster)
        pred_test_ok <- raster::extract(ok_raster, test)
    
        # res_test_ok <- test$sim1 - pred_test_ok
        # summary(res_test_ok)
        # hist(res_test_ok)
    
        ME = mean((test$sim1 - pred_test_ok), na.rm=TRUE)
        MAE = mean(abs(test$sim1 - pred_test_ok), na.rm=TRUE)
        RMSE = sqrt(mean((test$sim1 - pred_test_ok)^2, na.rm=TRUE))
        R2 = 1 - (t(test$sim1 - pred_test_ok) %*% (test$sim1 - pred_test_ok)) / (t(test$sim1 - mean(test$sim1)) %*% (test$sim1 - mean(test$sim1)))
        ccc = DescTools::CCC(test$sim1, pred_test_ok, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c$est
    
        # R2; RMSE; MAE; ccc
        seed_result <- rbind(seed_result, c(R2, "R2", nug, range, n_po, seed, "OK"))
        seed_result <- rbind(seed_result, c(RMSE, "RMSE", nug, range, n_po, seed, "OK"))
        seed_result <- rbind(seed_result, c(ccc, "CCC", nug, range, n_po, seed, "OK"))
        seed_result <- rbind(seed_result, c(MAE, "MAE", nug, range, n_po, seed, "OK"))
        seed_result <- rbind(seed_result, c(ME, "ME", nug, range, n_po, seed, "OK"))
    
        ##########################################################################################
        ### RFSI #################################################################################
        ##########################################################################################
    
        ### observations at nearest stations and deistances to them - calculation ###
        st_num <- length(sim_ok)
        sim_df <- as.data.frame(sim_ok)
    
        time1 <- Sys.time()
        nearest_obs <- near.obs(
          locations = sim_ok,
          observations = train,
          zcol = "sim1",
          n.obs = n_st
        )
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "DT", nug, range, n_po, seed, "RFSI"))
    
        sim_df <- cbind(sim_df, nearest_obs)
        sim_df.dev <- sim_df[sim_df$fold,]
    
        fm.obs = paste(names(nearest_obs), collapse="+")
        fm.RFSI <- as.formula(paste("sim1 ~ ", fm.obs))
        
        ### RFSI model ###
        set.seed(seed)
        time1 <- Sys.time()
        rf_RFSI <- ranger(fm.RFSI, data = sim_df.dev, importance = "impurity", num.trees = 250) #, mtry = n_st) # mtry - kriging range
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "MT", nug, range, n_po, seed, "RFSI"))
        rf_RFSI
    
        ### validation ###
        time1 <- Sys.time()
        rfsi_plot <- predict(rf_RFSI, sim_df) # doesn't take training stations for their predictions!
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "PT", nug, range, n_po, seed, "RFSI"))
        rfsi_raster <- fishnet
        values(rfsi_raster) <- rfsi_plot$predictions
        # plot(rfsi_raster)
    
        pred_test_rfsi <- raster::extract(rfsi_raster, test)
        # res_test_rfsi <- test$sim1 - pred_test_rfsi
        # summary(res_test_rfsi)
        # hist(res_test_ok)
        # hist(res_test_rfsi)
    
        ME = mean((test$sim1 - pred_test_rfsi), na.rm=TRUE)
        MAE = mean(abs(test$sim1 - pred_test_rfsi), na.rm=TRUE)
        RMSE = sqrt(mean((test$sim1 - pred_test_rfsi)^2, na.rm=TRUE))
        R2 = 1 - (t(test$sim1 - pred_test_rfsi) %*% (test$sim1 - pred_test_rfsi)) / (t(test$sim1 - mean(test$sim1)) %*% (test$sim1 - mean(test$sim1)))
        ccc = DescTools::CCC(test$sim1, pred_test_rfsi, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c$est
    
        # R2; RMSE; MAE; ccc
        seed_result <- rbind(seed_result, c(R2, "R2", nug, range, n_po, seed, "RFSI"))
        seed_result <- rbind(seed_result, c(RMSE, "RMSE", nug, range, n_po, seed, "RFSI"))
        seed_result <- rbind(seed_result, c(ccc, "CCC", nug, range, n_po, seed, "RFSI"))
        seed_result <- rbind(seed_result, c(MAE, "MAE", nug, range, n_po, seed, "RFSI"))
        seed_result <- rbind(seed_result, c(ME, "ME", nug, range, n_po, seed, "RFSI"))
    
        ##########################################################################################
        ### Nearest neighbour ####################################################################
        ##########################################################################################
    
        gs <- gstat(formula=sim1~1, locations=train, nmax=1, set=list(idp = 0))
        nn_raster <- interpolate(fishnet, gs)
        # plot(nn_raster)
    
        pred_test_nn <- raster::extract(nn_raster, test)
        # res_test_nn <- test$sim1 - pred_test_nn
        # summary(res_test_nn)
    
        ### validation ###
        ME = mean((test$sim1 - pred_test_nn), na.rm=TRUE)
        MAE = mean(abs(test$sim1 - pred_test_nn), na.rm=TRUE)
        RMSE = sqrt(mean((test$sim1 - pred_test_nn)^2, na.rm=TRUE))
        R2 = 1 - (t(test$sim1 - pred_test_nn) %*% (test$sim1 - pred_test_nn)) / (t(test$sim1 - mean(test$sim1)) %*% (test$sim1 - mean(test$sim1)))
        ccc = DescTools::CCC(test$sim1, pred_test_nn, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c$est
        
        # R2; RMSE; MAE; ccc
        seed_result <- rbind(seed_result, c(R2, "R2", nug, range, n_po, seed, "NN"))
        seed_result <- rbind(seed_result, c(RMSE, "RMSE", nug, range, n_po, seed, "NN"))
        seed_result <- rbind(seed_result, c(ccc, "CCC", nug, range, n_po, seed, "NN"))
        seed_result <- rbind(seed_result, c(MAE, "MAE", nug, range, n_po, seed, "NN"))
        seed_result <- rbind(seed_result, c(ME, "ME", nug, range, n_po, seed, "NN"))
    
        ##########################################################################################
        ### IDW ##################################################################################
        ##########################################################################################
    
        # better results with nmax (local), worse with all stations (global)!!!
        idw_plot = gstat::idw(sim1~1, train, st, na.action = na.pass, idp = 2, nmax=n_st) # , idp = 3 # , idp = 4
        idw_raster <- fishnet
        values(idw_raster) <- idw_plot$var1.pred
        # plot(idw_raster)
    
        pred_test_idw <- raster::extract(idw_raster, test)
        # res_test_idw <- test$sim1 - pred_test_idw
        # summary(res_test_idw)
    
        ### validation ###
        ME = mean((test$sim1 - pred_test_idw), na.rm=TRUE)
        MAE = mean(abs(test$sim1 - pred_test_idw), na.rm=TRUE)
        RMSE = sqrt(mean((test$sim1 - pred_test_idw)^2, na.rm=TRUE))
        R2 = 1 - (t(test$sim1 - pred_test_idw) %*% (test$sim1 - pred_test_idw)) / (t(test$sim1 - mean(test$sim1)) %*% (test$sim1 - mean(test$sim1)))
        ccc = DescTools::CCC(test$sim1, pred_test_idw, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c$est
    
        # R2; RMSE; MAE; ccc
        seed_result <- rbind(seed_result, c(R2, "R2", nug, range, n_po, seed, "IDW"))
        seed_result <- rbind(seed_result, c(RMSE, "RMSE", nug, range, n_po, seed, "IDW"))
        seed_result <- rbind(seed_result, c(ccc, "CCC", nug, range, n_po, seed, "IDW"))
        seed_result <- rbind(seed_result, c(MAE, "MAE", nug, range, n_po, seed, "IDW"))
        seed_result <- rbind(seed_result, c(ME, "ME", nug, range, n_po, seed, "IDW"))
    
        ##########################################################################################
        ### Trend surface - 2nd order ############################################################
        ##########################################################################################
    
        sim_df <- as.data.frame(sim_ok)
        sim_df.dev <- sim_df[sim_df$fold,]
        
        m_ts2 = lm(sim1 ~ I(x^2) + I(y^2) + I(x*y) + x + y, sim_df.dev)
    
        ts2_plot <- predict(m_ts2, sim_df)
        ts2_raster <- fishnet
        values(ts2_raster) <- ts2_plot
        # plot(ts2_raster)
    
        pred_test_ts2 <- raster::extract(ts2_raster, test)
        # res_test_ts2 <- test$sim1 - pred_test_ts2
        # summary(res_test_ts2)
    
        ### validation ###
        ME = mean((test$sim1 - pred_test_ts2), na.rm=TRUE)
        MAE = mean(abs(test$sim1 - pred_test_ts2), na.rm=TRUE)
        RMSE = sqrt(mean((test$sim1 - pred_test_ts2)^2, na.rm=TRUE))
        R2 = 1 - (t(test$sim1 - pred_test_ts2) %*% (test$sim1 - pred_test_ts2)) / (t(test$sim1 - mean(test$sim1)) %*% (test$sim1 - mean(test$sim1)))
        ccc = DescTools::CCC(test$sim1, pred_test_ts2, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c$est
    
        # R2; RMSE; MAE; ccc
        seed_result <- rbind(seed_result, c(R2, "R2", nug, range, n_po, seed, "TS2"))
        seed_result <- rbind(seed_result, c(RMSE, "RMSE", nug, range, n_po, seed, "TS2"))
        seed_result <- rbind(seed_result, c(ccc, "CCC", nug, range, n_po, seed, "TS2"))
        seed_result <- rbind(seed_result, c(MAE, "MAE", nug, range, n_po, seed, "TS2"))
        seed_result <- rbind(seed_result, c(ME, "ME", nug, range, n_po, seed, "TS2"))
        # print(R2)
    
        ##########################################################################################
        ### RFsp (Hengl et al. 2018) #############################################################
        ##########################################################################################
    
        co_grids = as(sim_ok, "SpatialPixelsDataFrame")
        # identical(co_grids@data$sim1, sim_df$sim1)
        time1 <- Sys.time()
        grid.distP <- GSIF::buffer.dist(train, co_grids[1], as.factor(1:nrow(train))) # Hengl's buffer distances
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "DT", nug, range, n_po, seed, "RFsp"))
        sim_df <- as.data.frame(sim_ok)
        sim_df <- cbind(sim_df, grid.distP@data)
        
        # remove tmp files of GSIF:buffer.dist function #
        tmp_files <- list.files(dirname(rasterTmpFile()), full.names = T)
        tmp_files_time <- file.info(tmp_files)$ctime
        if (n_po > 2000) {p_sec = 2*7200} else {p_sec = 600}
        tmp_files <- tmp_files[tmp_files_time < (Sys.time()-p_sec)]
        unlink(tmp_files, recursive = TRUE)
    
        sim_df.dev <- sim_df[sim_df$fold,]
    
        fm.dist <- paste(names(grid.distP), collapse="+")
        fm.RFsp <- as.formula(paste("sim1 ~ ", fm.dist))
    
        # RFsp model #
        set.seed(seed)
        time1 <- Sys.time()
        rf_RFsp <- ranger(fm.RFsp, data = sim_df.dev, importance = "impurity", num.trees = 250, mtry = 2*n_po/3)
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "MT", nug, range, n_po, seed, "RFsp"))
        # rf_RFsp
    
        ### validation ###
        time1 <- Sys.time()
        rfsp_plot <- predict(rf_RFsp, sim_df)
        seed_result <- rbind(seed_result, c(as.numeric(Sys.time()-time1,units="secs"), "PT", nug, range, n_po, seed, "RFsp"))
        rfsp_raster <- fishnet
        values(rfsp_raster) <- rfsp_plot$predictions
        # plot(rfsp_raster)
    
        pred_test_rfsp <- raster::extract(rfsp_raster, test)
        # res_test_rfsp <- test$sim1 - pred_test_rfsp
        # summary(res_test_rfsi)
        # hist(res_test_ok)
        # hist(res_test_rfsi)
    
        ME = mean((test$sim1 - pred_test_rfsp), na.rm=TRUE)
        MAE = mean(abs(test$sim1 - pred_test_rfsp), na.rm=TRUE)
        RMSE = sqrt(mean((test$sim1 - pred_test_rfsp)^2, na.rm=TRUE))
        R2 = 1 - (t(test$sim1 - pred_test_rfsp) %*% (test$sim1 - pred_test_rfsp)) / (t(test$sim1 - mean(test$sim1)) %*% (test$sim1 - mean(test$sim1)))
        ccc = DescTools::CCC(test$sim1, pred_test_rfsp, ci = "z-transform", conf.level = 0.95, na.rm=TRUE)$rho.c$est
    
        # R2; RMSE; MAE; ccc
        seed_result <- rbind(seed_result, c(R2, "R2", nug, range, n_po, seed, "RFsp"))
        seed_result <- rbind(seed_result, c(RMSE, "RMSE", nug, range, n_po, seed, "RFsp"))
        seed_result <- rbind(seed_result, c(ccc, "CCC", nug, range, n_po, seed, "RFsp"))
        seed_result <- rbind(seed_result, c(MAE, "MAE", nug, range, n_po, seed, "RFsp"))
        seed_result <- rbind(seed_result, c(ME, "ME", nug, range, n_po, seed, "RFsp"))
    
        gc();gc()
        rm(grid.distP)
        
        # ##########################################################################################
        # ### Predictions ##########################################################################
        # ##########################################################################################
        # 
        # dir.create("plot")
        # 
        # all_plot <- as(ok_raster, "SpatialPixelsDataFrame")
        # names(all_plot) <- "OK"
        # all_plot$SIM <- sim_ok$sim1
        # all_plot@data <- all_plot@data[, 2:1]
        # all_plot$RFSI <- values(rfsi_raster)
        # # all_plot$RFSIavg <- values(rfsi_avg_raster)
        # all_plot$RFsp <- values(rfsp_raster)
        # all_plot$IDW <- values(idw_raster)
        # all_plot$NN <- values(nn_raster)
        # all_plot$TS <- values(ts2_raster)
        # # all_plot$fold <- sim_ok$fold
        # pts <- train
        # 
        # ### Fig5 ###
        # # tiff(paste("plot/pred_comparison_", n_po, "_s", seed, "_nug", nug, ".tif", sep = ""), width = 174, height = 93, units = 'mm', res = 600, compression = "lzw")
        # jpeg(paste("plot/pred_comparison_", n_po, "_s", seed, "_nug", nug, ".jpeg", sep = ""), width = 174, height = 93, units = 'mm', res = 600) # res = 1200
        # print(spplot(all_plot, col.regions=terrain.colors, par.settings=list(fontsize=list(text=7)),
        #              sp.layout = list(pts, which=1, col="black", cex=0.1, alpha=0.5),
        #              # scales = list(draw = TRUE),
        #              as.table=TRUE
        #              )) #, main = "Comparison of predictions"
        # dev.off()
        # 
        # ##########################################################################################
        
        # ##########################################################################################
        # ### plot residuals with all combinations #################################################
        # ##########################################################################################
        #
        # test$OK <- (all_plot$SIM - all_plot$OK)[!sim_ok$fold]
        # test$RFSI <- (all_plot$SIM - all_plot$RFSI)[!sim_ok$fold]
        # # test$RFSIavg <- res_test_rf_avg
        # test$RFsp <- (all_plot$SIM - all_plot$RFsp)[!sim_ok$fold]
        # test$IDW <- (all_plot$SIM - all_plot$IDW)[!sim_ok$fold]
        # test$NN <- (all_plot$SIM - all_plot$NN)[!sim_ok$fold]
        # test$TS2 <- (all_plot$SIM - all_plot$TS2)[!sim_ok$fold]
        # 
        # pal <-  brewer.pal(n = 11, name = "RdYlBu")
        # # pal <- pal[-c(5,7)]
        # pal[6] <- "white"
        # cut <- c(min(test@data[,3:ncol(test)]), -7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7, max(test@data[,3:6]))
        # # tiff(paste("plot/res_comparison_", n_po, "_s", seed, "_nug", nug, ".tif", sep = ""), width = 174, height = 120, units = 'mm', res = 600, compression = "lzw")
        # jpeg(paste("plot/res_comparison_", n_po, "_s", seed, "_nug", nug, ".jpeg", sep = ""), width = 174, height = 120, units = 'mm', res = 600) # res = 1200
        # print(spplot(as(test, "SpatialPixelsDataFrame"), zcol = 3:ncol(test), col.regions = pal, at = cut,#col.regions = bpy.colors(64),
        #              par.settings=list(fontsize=list(text=7)), as.table=TRUE #,
        #      # main = "Comparison of residuals",
        #      # sub = paste("OK: R2 = ", round(r2, 2), ", RMSE = ", round(rmse, 2), ", Lin = ", round(lin, 2),
        #      #        "\nRFSI: R2 = ", round(r2_1, 2), ", RMSE = ", round(rmse_1,2), ", Lin = ", round(lin_1, 2),
        #      #        "\nNN: R2 = ", round(r2_3, 2), ", RMSE = ", round(rmse_3,2), ", Lin = ", round(lin_3, 2),
        #      #        "\nIDW: R2 = ", round(r2_4, 2), ", RMSE = ", round(rmse_4,2), ", Lin = ", round(lin_4, 2),
        #      #        "\nTS2: R2 = ", round(r2_6, 2), ", RMSE = ", round(rmse_6,2), ", Lin = ", round(lin_6, 2),
        #      #        "\nRFsp: R2 = ", round(r2_7, 2), ", RMSE = ", round(rmse_7,2), ", Lin = ", round(lin_7, 2),
        #             # sep = "")
        #      ))
        # dev.off()
        # ##########################################################################################
        
        # ###### importance plot ###################################################################
        # 
        # # RFSI #
        # xlP.g <- as.list(round(rf_RFSI$variable.importance))
        # df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:15,]/100000
        # print(df)
        # pr <- 15:1
        # df <- as.data.frame(cbind(cbind(df, names(df)), pr))
        # df$V2 <- as.character(df$V2)
        # names(df) <- c("importance", "covariate", "order")
        # 
        # # normalize
        # df$importance <- as.numeric(as.character(df$importance))
        # # df$importance <- (df$importance-min(df$importance))/(max(df$importance)-min(df$importance))
        # df$importance <- df$importance / max(df$importance)
        # summary(df$importance)
        # 
        # # Reorder the data
        # df <- df %>%
        #   arrange(importance) %>%
        #   mutate(covariate=factor(covariate,covariate))
        # 
        # theme = theme_set(theme_minimal())
        # 
        # # Plot
        # rfsi_importance <-
        #   ggplot(df, aes(x=covariate, y=importance)) +
        #   geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
        #   geom_point( color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
        #   # theme_ipsum() +
        #   coord_flip() +
        #   theme(
        #     legend.position="none"
        #   ) +
        #   xlab("Covariate") +
        #   ylab("Importance") +
        #   ggtitle("RFSI") +
        #   theme(plot.title = element_text(hjust = 0.5),
        #         axis.text = element_text(size = 5),
        #         axis.title = element_text(size = 7),
        #         text = element_text(size = 7))
        # 
        # # RFsp #
        # xlP.g <- as.list(round(rf_RFsp$variable.importance))
        # df = t(data.frame(xlP.g[order(unlist(xlP.g), decreasing=TRUE)]))[1:15,]/100000
        # print(df)
        # pr <- 15:1
        # df <- as.data.frame(cbind(cbind(df, names(df)), pr))
        # df$V2 <- as.character(df$V2)
        # names(df) <- c("importance", "covariate", "order")
        # 
        # # normalize
        # df$importance <- as.numeric(as.character(df$importance))
        # # df$importance <- (df$importance-min(df$importance))/(max(df$importance)-min(df$importance))
        # df$importance <- df$importance / max(df$importance)
        # summary(df$importance)
        # 
        # # Reorder the data
        # df <- df %>%
        #   arrange(importance) %>%
        #   mutate(covariate=factor(covariate,covariate))
        # 
        # rfsp_importance <-
        #   ggplot(df, aes(x=covariate, y=importance)) +
        #   geom_segment( aes(x=covariate, xend=covariate, y=0, yend=importance ), color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
        #   geom_point( color=ifelse(df$covariate %in% c("tmin","tmax","imerg"), "red", "black"), size=ifelse(df$covariate %in% c("tmin","tmax","imerg"), .5, .25) ) +
        #   # theme_ipsum() +
        #   coord_flip() +
        #   theme(
        #     legend.position="none"
        #   ) +
        #   xlab("Covariate") +
        #   ylab("Importance") +
        #   ggtitle("RFsp") +
        #   theme(plot.title = element_text(hjust = 0.5),
        #         axis.text = element_text(size = 5),
        #         axis.title = element_text(size = 7),
        #         text = element_text(size = 7))
        # 
        # ### Fig6 ###
        # # tiff("plot/importance.tiff", width = 75, height = 70, units = 'mm', res = 1200, compression = "lzw")
        # jpeg("plot/importance.jpeg", width = 75, height = 70, units = 'mm', res = 1200)
        # ggarrange(rfsp_importance, rfsi_importance, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
        # dev.off()
        # ##########################################################################################
        
      return(seed_result)
      }
      stopImplicitCluster()
      
      accuracy_results <- rbind(accuracy_results, do.call("rbind", seed_results))
    }
  }
}

# save(accuracy_results, file = paste("accuracy_results.rda", sep = ""))
summary(accuracy_results)
