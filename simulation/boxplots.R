library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

load("accuracy_results.rda")
accuracy_results <- as.data.frame(accuracy_results)
names(accuracy_results) <- c("value", "parameter", "nugget", "range", "points", "seed", "method")
accuracy_results$value <- as.numeric(as.character(accuracy_results$value))
accuracy_results$nugget <- as.numeric(as.character(accuracy_results$nugget))
accuracy_results$range <- as.integer(as.character(accuracy_results$range))
accuracy_results$seed <- as.integer(as.character(accuracy_results$seed))

# accuracy_results[accuracy_results$nugget == 2.5 &
#                    accuracy_results$range == 50 &
#                    accuracy_results$points == 500 &
#                    accuracy_results$seed == 915 &
#                    accuracy_results$parameter == "R2", ]

n_st <- 25
##################
nugs <- c(0, 2.5, 5)
ranges <- c(50, 200)

unique(accuracy_results$parameter)

n_pois <- c(100, 200, 500, 1000, 2000, 5000)

for (nug in nugs) {
  for (range in ranges) {
    
    # # ### NUGGET ###
    # nug <- 2.5
    # # ##############
    # # 
    # # ### RANGE ###
    # range <- 200
    # # ##############
    
    ### box plot R2 vs number of points ###
    
    r2_all <- accuracy_results[accuracy_results$parameter == "R2" & accuracy_results$nugget == nug & accuracy_results$range == range, c(1,5,6,7)]
    names(r2_all) <- c("R2", "Npoints", "Seed", "Method")
    # r2_all$Point_percent <- as.factor(r2_all$Point_percent)
    r2_all$R2 <- round(r2_all$R2, 5)*100
    
    r2_all$Method = factor(r2_all$Method,levels(r2_all$Method)[c(3,4,5,1,2,6)])
    head(r2_all)
    
    # Calculates mean, sd, se and IC
    my_sum <- r2_all %>%
      group_by(Npoints, Method) %>%
      summarise( 
        n=n(),
        mean=mean(R2),
        median=median(R2),
        sd=sd(R2)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
    
    ymin <- min(my_sum$mean - my_sum$se)
    ymax <- max(my_sum$mean + my_sum$se)
  
    theme = theme_set(theme_minimal())
    assign(paste("r2.", nug, ".", range, sep=""),
    ggplot(data = my_sum, aes(x=Npoints, y=mean, fill=Method)) +
      geom_bar(position="dodge", stat = "identity") +
      geom_errorbar( aes(x=Npoints, ymin=mean-se, ymax=mean+se), width=0.3, colour="black", alpha=0.9, size=0.3, position=position_dodge(.9)) +
      scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
      scale_fill_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                        values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
      scale_color_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                         values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
      geom_vline(xintercept = 1.5, linetype="dashed", size=0.2) +
      geom_vline(xintercept = 2.5, linetype="dashed", size=0.2) +
      geom_vline(xintercept = 3.5, linetype="dashed", size=0.2) +
      geom_vline(xintercept = 4.5, linetype="dashed", size=0.2) +
      geom_vline(xintercept = 5.5, linetype="dashed", size=0.2) +
      geom_vline(xintercept = 6.5, linetype="dashed", size=0.2) +
      # geom_vline(xintercept = 7.5, linetype="dashed", size=0.2) +
      
      xlab("Number of sample locations") +
      ylab(expression(R^{2})) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      
      theme(plot.title = element_text(hjust = 8),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 8),
            text = element_text(size = 8),
            # legend.position = "bottom",
            # legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size=8, face="bold"),
            legend.text=element_text(size=8))
    )
    
    ### box plot RMSE vs number of points ###

    rmse_all <- accuracy_results[accuracy_results$parameter == "RMSE" & accuracy_results$nugget == nug & accuracy_results$range == range, c(1,5,6,7)]
    names(rmse_all) <- c("RMSE", "Npoints", "Seed", "Method")
    # rmse_all$Point_percent <- as.factor(rmse_all$Point_percent)
    rmse_all$RMSE <- round(as.numeric(as.character(rmse_all$RMSE)), 3)

    rmse_all$Method = factor(rmse_all$Method,levels(rmse_all$Method)[c(3,4,5,1,2,6)])
    head(rmse_all)
    
    # Calculates mean, sd, se and IC
    my_sum <- rmse_all %>%
      group_by(Npoints, Method) %>%
      summarise( 
        n=n(),
        mean=mean(RMSE),
        median=median(RMSE),
        sd=sd(RMSE)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
    
    ymin <- min(my_sum$mean - my_sum$se)
    ymax <- max(my_sum$mean + my_sum$se)
    
    theme = theme_set(theme_minimal())
    assign(paste("rmse.", nug, ".", range, sep=""),
           ggplot(data = my_sum, aes(x=Npoints, y=mean, fill=Method)) +
             geom_bar(position="dodge", stat = "identity") +
             geom_errorbar( aes(x=Npoints, ymin=mean-se, ymax=mean+se), width=0.3, colour="black", alpha=0.9, size=0.3, position=position_dodge(.9)) +
             scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
             scale_fill_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                               values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             scale_color_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                                values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             geom_vline(xintercept = 1.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 2.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 3.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 4.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 5.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 6.5, linetype="dashed", size=0.2) +
             # geom_vline(xintercept = 7.5, linetype="dashed", size=0.2) +
             
             xlab("Number of sample locations") +
             ylab("RMSE") +
             coord_cartesian(ylim = c(ymin, ymax)) +
             
             theme(plot.title = element_text(hjust = 8),
                   axis.text = element_text(size = 8),
                   axis.title = element_text(size = 8),
                   text = element_text(size = 8),
                   # legend.position = "bottom",
                   # legend.direction = "horizontal",
                   legend.key.size= unit(0.2, "cm"),
                   legend.margin = unit(0, "cm"),
                   legend.title = element_text(size=8, face="bold"),
                   legend.text=element_text(size=8))
    )

    ### box plot CCC vs number of points ###

    lin_all <- accuracy_results[accuracy_results$parameter == "CCC" & accuracy_results$nugget == nug & accuracy_results$range == range, c(1,5,6,7)]
    names(lin_all) <- c("CCC", "Npoints", "Seed", "Method")
    # lin_all$Point_percent <- as.factor(lin_all$Point_percent)
    lin_all$CCC <- round(as.numeric(as.character(lin_all$CCC)), 3)

    lin_all$Method = factor(lin_all$Method,levels(lin_all$Method)[c(3,4,5,1,2,6)])
    head(lin_all)
    
    my_sum <- lin_all %>%
      group_by(Npoints, Method) %>%
      summarise( 
        n=n(),
        mean=mean(CCC),
        median=median(CCC),
        sd=sd(CCC)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
    
    ymin <- min(my_sum$mean - my_sum$se)
    ymax <- max(my_sum$mean + my_sum$se)
    
    theme = theme_set(theme_minimal())
    assign(paste("lin.", nug, ".", range, sep=""),
           ggplot(data = my_sum, aes(x=Npoints, y=mean, fill=Method)) +
             geom_bar(position="dodge", stat = "identity") +
             geom_errorbar( aes(x=Npoints, ymin=mean-se, ymax=mean+se), width=0.3, colour="black", alpha=0.9, size=0.3, position=position_dodge(.9)) +
             scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
             scale_fill_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                               values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             scale_color_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                                values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             geom_vline(xintercept = 1.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 2.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 3.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 4.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 5.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 6.5, linetype="dashed", size=0.2) +
             # geom_vline(xintercept = 7.5, linetype="dashed", size=0.2) +
             
             xlab("Number of sample locations") +
             ylab("CCC") +
             coord_cartesian(ylim = c(ymin, ymax)) +
             
             theme(plot.title = element_text(hjust = 8),
                   axis.text = element_text(size = 8),
                   axis.title = element_text(size = 8),
                   text = element_text(size = 8),
                   # legend.position = "bottom",
                   # legend.direction = "horizontal",
                   legend.key.size= unit(0.2, "cm"),
                   legend.margin = unit(0, "cm"),
                   legend.title = element_text(size=8, face="bold"),
                   legend.text=element_text(size=8))
    )
    
    ### box plot MAE vs number of points ###
    
    mae_all <- accuracy_results[accuracy_results$parameter == "MAE" & accuracy_results$nugget == nug & accuracy_results$range == range, c(1,5,6,7)]
    names(mae_all) <- c("MAE", "Npoints", "Seed", "Method")
    mae_all$MAE <- round(as.numeric(as.character(mae_all$MAE)), 3)
    
    mae_all$Method = factor(mae_all$Method,levels(mae_all$Method)[c(3,4,5,1,2,6)])
    head(mae_all)
    
    my_sum <- mae_all %>%
      group_by(Npoints, Method) %>%
      summarise( 
        n=n(),
        mean=mean(MAE),
        median=median(MAE),
        sd=sd(MAE)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
    
    ymin <- min(my_sum$mean - my_sum$se)
    ymax <- max(my_sum$mean + my_sum$se)
    
    theme = theme_set(theme_minimal())
    assign(paste("mae.", nug, ".", range, sep=""),
           ggplot(data = my_sum, aes(x=Npoints, y=mean, fill=Method)) +
             geom_bar(position="dodge", stat = "identity") +
             geom_errorbar( aes(x=Npoints, ymin=mean-se, ymax=mean+se), width=0.3, colour="black", alpha=0.9, size=0.3, position=position_dodge(.9)) +
             scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
             scale_fill_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                               values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             scale_color_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                                values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             geom_vline(xintercept = 1.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 2.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 3.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 4.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 5.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 6.5, linetype="dashed", size=0.2) +
             # geom_vline(xintercept = 7.5, linetype="dashed", size=0.2) +
             
             xlab("Number of sample locations") +
             ylab("MAE") +
             coord_cartesian(ylim = c(ymin, ymax)) +
             
             theme(plot.title = element_text(hjust = 8),
                   axis.text = element_text(size = 8),
                   axis.title = element_text(size = 8),
                   text = element_text(size = 8),
                   # legend.position = "bottom",
                   # legend.direction = "horizontal",
                   legend.key.size= unit(0.2, "cm"),
                   legend.margin = unit(0, "cm"),
                   legend.title = element_text(size=8, face="bold"),
                   legend.text=element_text(size=8))
    )
    
    ### box plot ME vs number of points ###
    
    me_all <- accuracy_results[accuracy_results$parameter == "ME" & accuracy_results$nugget == nug & accuracy_results$range == range, c(1,5,6,7)]
    names(me_all) <- c("ME", "Npoints", "Seed", "Method")
    me_all$ME <- round(as.numeric(as.character(me_all$ME)), 3)
    
    me_all$Method = factor(me_all$Method,levels(me_all$Method)[c(3,4,5,1,2,6)])
    head(me_all)
    
    my_sum <- me_all %>%
      group_by(Npoints, Method) %>%
      summarise( 
        n=n(),
        mean=mean(ME),
        median=median(ME),
        sd=sd(ME)
      ) %>%
      mutate( se=sd/sqrt(n))  %>%
      mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
    
    ymin <- min(my_sum$mean - my_sum$se)
    ymax <- max(my_sum$mean + my_sum$se)
    
    theme = theme_set(theme_minimal())
    assign(paste("me.", nug, ".", range, sep=""),
           ggplot(data = my_sum, aes(x=Npoints, y=mean, fill=Method)) +
             geom_bar(position="dodge", stat = "identity") +
             geom_errorbar( aes(x=Npoints, ymin=mean-se, ymax=mean+se), width=0.3, colour="black", alpha=0.9, size=0.3, position=position_dodge(.9)) +
             scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
             scale_fill_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                               values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             scale_color_manual(breaks = c("OK", "RFSI", "RFsp", "IDW", "NN", "TS2"),
                                values=c("red", "green", "deepskyblue", "blueviolet", "chocolate1", "bisque3"))+#, "gold")) +
             geom_vline(xintercept = 1.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 2.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 3.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 4.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 5.5, linetype="dashed", size=0.2) +
             geom_vline(xintercept = 6.5, linetype="dashed", size=0.2) +
             # geom_vline(xintercept = 7.5, linetype="dashed", size=0.2) +
             
             xlab("Number of sample locations") +
             ylab("ME") +
             coord_cartesian(ylim = c(ymin, ymax)) +
             
             theme(plot.title = element_text(hjust = 8),
                   axis.text = element_text(size = 8),
                   axis.title = element_text(size = 8),
                   text = element_text(size = 8),
                   # legend.position = "bottom",
                   # legend.direction = "horizontal",
                   legend.key.size= unit(0.2, "cm"),
                   legend.margin = unit(0, "cm"),
                   legend.title = element_text(size=8, face="bold"),
                   legend.text=element_text(size=8))
    )
    
    }
}

dir.create("boxplots")
# # tiff("boxplots/r2.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/r2.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
# ggarrange(r2.0.50, r2.0.200, r2.2.5.50, r2.2.5.200, r2.5.50, r2.5.200, ncol=2, nrow=3, common.legend = TRUE, legend="right")
# dev.off()
# 
# # tiff("boxplots/rmse.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/rmse.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
# ggarrange(rmse.0.50, rmse.0.200, rmse.2.5.50, rmse.2.5.200, rmse.5.50, rmse.5.200, ncol=2, nrow=3, common.legend = TRUE, legend="right")
# dev.off()
# 
# # tiff("boxplots/lin.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/lin.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
# ggarrange(lin.0.50, lin.0.200, lin.2.5.50, lin.2.5.200, lin.5.50, lin.5.200, ncol=2, nrow=3, common.legend = TRUE, legend="right")
# dev.off()

### Fig2 ###
# tiff("boxplots/mae.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
jpeg("boxplots/mae.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
f1=ggarrange(mae.0.50, mae.0.200, mae.2.5.50, mae.2.5.200, mae.5.50, mae.5.200, ncol=2, nrow=3,
          common.legend = TRUE, legend="right")
annotate_figure(f1,
                top = text_grob("Range 50                                                               Range 200", face = "bold", size = 10),
                left = text_grob("  Nugget-to-sill ratio 0.50                 Nugget-to-sill ratio 0.25                  Nugget-to-sill ratio 0.00", face = "bold", rot = 90, size = 10)
)
dev.off()

# # tiff("boxplots/me.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/me.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
# ggarrange(me.0.50, me.0.200, me.2.5.50, me.2.5.200, me.5.50, me.5.200, ncol=2, nrow=3, common.legend = TRUE, legend="right")
# dev.off()

### Fig3 ###
# tiff("boxplots/2.5.200.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
jpeg("boxplots/2.5.200.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
ggarrange(r2.2.5.200, lin.2.5.200, mae.2.5.200, rmse.2.5.200, ncol=2, nrow=2, common.legend = TRUE, legend="right")
dev.off()


### box plot modeling time (DT) vs number of points ###

dt_all <- accuracy_results[accuracy_results$parameter == "DT", c(1,5,6,7)]
names(dt_all) <- c("DT", "Npoints", "Seed", "Method")
dt_all$DT <- round(as.numeric(as.character(dt_all$DT)), 3)

dt_all$Method = factor(dt_all$Method,levels(dt_all$Method)[c(3,4,5,1,2,6)])
head(dt_all)

my_sum <- dt_all %>%
  group_by(Npoints, Method) %>%
  summarise( 
    n=n(),
    mean=mean(DT),
    median=median(DT),
    sd=sd(DT)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

theme = theme_set(theme_minimal())
assign(paste("dt", sep=""),
         ggplot(data = my_sum, aes(x=Npoints, y=mean, group=c(Method))) +
         geom_line(aes(color=Method), size=0.5) +
         geom_point(aes(color=Method), size=0.8) +
         geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),
                       alpha=0.2) +
         geom_label_repel(aes(label=round(mean, 2), fill = Method), size=2.5, colour = "white", fontface = "bold") + 
         scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
         scale_fill_manual(breaks = c("RFSI", "RFsp"),
                           values=c("chartreuse4", "deepskyblue")) +
         scale_color_manual(breaks = c("RFSI", "RFsp"),
                            values=c("chartreuse4", "deepskyblue")) +
         xlab("Number of sample locations") +
         ylab("Distance calculation time [sec]") +
         theme(plot.title = element_text(hjust = 8),
               axis.text = element_text(size = 8),
               axis.title = element_text(size = 8),
               text = element_text(size = 8),
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(size=8, face="bold"),
               legend.text=element_text(size=8))
)

dt
# # tiff("boxplots/dt.tiff", width = 174, height = 87, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/dt.jpeg", width = 174, height = 87, units = 'mm', res = 600) # res = 1200
# dt
# dev.off()


# my_sum1 <- my_sum[my_sum$Method == "RFsp", ]
my_sum1 <- my_sum[my_sum$Method == "RFSI", ]
my_sum1$Npoints <- as.numeric(as.character(my_sum1$Npoints))
my_sum1 <- my_sum1[order(my_sum1$Npoints), ]

model_exp <- lm(log(mean)~Npoints, my_sum1)
summary(model_exp)
# RFsp 0.8806 # 0.626
# RFSI 0.5035 # 0.05933
sum(model_exp$residuals^2)
# RFsp 1.567725 # RFSI 0.01407862

model_lm <- lm(mean~Npoints, my_sum1)
summary(model_lm)
coef(model_lm)
# RFsp 0.9092 # 465.4
# (Intercept)      Npoints 
# -357.2686107    0.7974501 
# RFSI 0.538 # 0.08974
# (Intercept)      Npoints 
# 1.516409e+00 5.620022e-05
sum(model_lm$residuals^2)
# RFsp 866570 # RFSI 0.03221077

# Plot fitted curve
plot(my_sum1$Npoints, my_sum1$mean)
lines(my_sum1$Npoints, predict(model_lm, list(Npoints = my_sum1$Npoints)), col = 'skyblue', lwd = 3)
lines(my_sum1$Npoints, predict(model_lm, list(Npoints = my_sum1$Npoints)), col = 'red', lwd = 3)

### box plot modeling time (MT) vs number of points ###

mt_all <- accuracy_results[accuracy_results$parameter == "MT", c(1,5,6,7)]
names(mt_all) <- c("MT", "Npoints", "Seed", "Method")
mt_all$MT <- round(as.numeric(as.character(mt_all$MT)), 3)

mt_all$Method = factor(mt_all$Method,levels(mt_all$Method)[c(3,4,5,1,2,6)])
head(mt_all)

my_sum <- mt_all %>%
  group_by(Npoints, Method) %>%
  summarise( 
    n=n(),
    mean=mean(MT),
    median=median(MT),
    sd=sd(MT)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

theme = theme_set(theme_minimal())
assign(paste("mt", sep=""),
       ggplot(data = my_sum, aes(x=Npoints, y=mean, group=c(Method))) +
         geom_line(aes(color=Method), size=0.5) +
         geom_point(aes(color=Method), size=0.8) +
         geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),
                     alpha=0.2) +
         geom_label_repel(aes(label=round(mean, 2), fill = Method), size=2.5, colour = "white", fontface = "bold") + 
         scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
         scale_fill_manual(breaks = c("RFSI", "RFsp"),
                           values=c("chartreuse4", "deepskyblue")) +
         scale_color_manual(breaks = c("RFSI", "RFsp"),
                            values=c("chartreuse4", "deepskyblue")) +
         xlab("Number of sample locations") +
         ylab("Modelling time [sec]") +
         theme(plot.title = element_text(hjust = 8),
               axis.text = element_text(size = 8),
               axis.title = element_text(size = 8),
               text = element_text(size = 8),
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(size=8, face="bold"),
               legend.text=element_text(size=8))
)

mt
# # tiff("boxplots/mt.tiff", width = 174, height = 87, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/mt.jpeg", width = 174, height = 87, units = 'mm', res = 600) # res = 1200
# mt
# dev.off()

### box plot modeling time (PT) vs number of points ###

pt_all <- accuracy_results[accuracy_results$parameter == "PT", c(1,5,6,7)]
names(pt_all) <- c("PT", "Npoints", "Seed", "Method")
pt_all$PT <- round(as.numeric(as.character(pt_all$PT)), 3)

pt_all$Method = factor(pt_all$Method,levels(pt_all$Method)[c(3,4,5,1,2,6)])
head(pt_all)

my_sum <- pt_all %>%
  group_by(Npoints, Method) %>%
  summarise( 
    n=n(),
    mean=mean(PT),
    median=median(PT),
    sd=sd(PT)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

theme = theme_set(theme_minimal())
assign(paste("pt", sep=""),
       ggplot(data = my_sum, aes(x=Npoints, y=mean, group=c(Method))) +
         geom_line(aes(color=Method), size=0.5) +
         geom_point(aes(color=Method), size=0.8) +
         geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),
                     alpha=0.2) +
         geom_label_repel(aes(label=round(mean, 2), fill = Method), size=2.5, colour = "white", fontface = "bold") + 
         scale_x_discrete(limits=as.character(n_pois), labels=as.character(n_pois)) +
         scale_fill_manual(breaks = c("OK", "RFSI", "RFsp"),
                           values=c("red", "chartreuse4", "deepskyblue")) +
         scale_color_manual(breaks = c("OK", "RFSI", "RFsp"),
                            values=c("red", "chartreuse4", "deepskyblue")) +
         xlab("Number of sample locations") +
         ylab("Prediction time [sec]") +
         theme(plot.title = element_text(hjust = 8),
               axis.text = element_text(size = 8),
               axis.title = element_text(size = 8),
               text = element_text(size = 8),
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(size=8, face="bold"),
               legend.text=element_text(size=8))
)

pt
# # tiff("boxplots/pt.tiff", width = 174, height = 87, units = 'mm', res = 600, compression = "lzw")
# jpeg("boxplots/pt.jpeg", width = 174, height = 87, units = 'mm', res = 600) # res = 1200
# pt
# dev.off()

### Fig4 ###
# tiff("boxplots/time.tiff", width = 174, height = 174, units = 'mm', res = 600, compression = "lzw")
jpeg("boxplots/time.jpeg", width = 174, height = 174, units = 'mm', res = 600) # res = 1200
ggarrange(dt, mt, pt, ncol=1, nrow=3) #, common.legend = TRUE, legend="right")
dev.off()
