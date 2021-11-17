#### hawaii coral bleaching analysis: spatial variability ####
## written by: morgan winston

## this script uses cluster-level data described & linked to here: {https://www.fisheries.noaa.gov/inport/item/64324}
## code performs the following: investigates differences b/w locations (zones in MHI, islands in NWHI); creates figures

#### initialization ####
# load packages
library(ggplot2) # for plotting figures
library(plotrix) # for calculating standard error
library(lme4) # for running linear mixed models
library(emmeans) # for tukey post hoc tests
library(multcomp) # for assigning letters of pair wise comparisons
library(data.table)
library(gridExtra)

# set working directory & load data
setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/For InPort/Environmental Drivers")
hcbc <- read.csv("HCBC_2019_ClusteredData.csv")
head(hcbc)

#### data prep ####
hcbc_clust <- hcbc[ which(hcbc$Obs_Year == "2019"),] # only looking @ 2019 variability across space here
hcbc_clust$PctBleached_mean_sqrt <- sqrt(hcbc_clust$CoralBleached_Perc_mn) # square-root transform response variable

# assign regions
hcbc_clust$Region_Name <- NA
for(i in c(1:nrow(hcbc_clust))){
  if(hcbc_clust$Island_Name[i] %in% c("French Frigate Shoals", "Lisianski", "Pearl and Hermes", "Kure")){
    hcbc_clust$Region_Name[i] <- "NWHI" 
  }
  else{
    hcbc_clust$Region_Name[i] <- "MHI"
  }
}

#### data analysis ####
# compare response across finest spatial scale possible (zones in MHI, isl/atolls in NWHI)
mod_sp <- aov(PctBleached_mean_sqrt ~ ZoneName, data = hcbc_clust, weights = Weights_DriversAndSpatialAnalysis)
summary(mod_sp) # significant difference p < 0.0001 ***

hcbc_clust$predict_sp <- predict(mod_sp, hcbc_clust, interval = "confidence") # get model predictions per zone
hcbc_clust$predict_ble <- hcbc_clust$predict_sp[,1]
hcbc_clust$error_sp <- hcbc_clust$predict_sp[,1]- hcbc_clust$predict_sp[,2]

# post hoc tests 
tuk <- TukeyHSD(mod_sp, "ZoneName") # tukey post hoc tests 
tuk <- as.data.frame(tuk$ZoneName) # for p-values
tuk[ which(tuk$`p adj` < 0.05),]
marginal = emmeans(mod_sp, ~ ZoneName) # determine the groups for labels on figure
cld <- data.frame(cld(marginal, Letters=letters))

# order of ZoneNames for plotting
ZoneName_ord <- c("Kure", "Pearl and Hermes", "Lisianski", "French Frigate Shoals", "Kauai_N", "Kauai_S", 
                  "Oahu_N", "Oahu_E", "Oahu_S", "Oahu_W", "Molokai_SE", "Molokai_S", "Lanai_NE", "Lanai_S", "Maui_N", "Maui_S", "Maui_SE", "Maui_W", "Maui_WNW", "Maui_NW",
                  "Hawaii_NW", "Hawaii_SW", "Hawaii_E")
IslandName_ord <- c("Kure", "Pearl and Hermes", "Lisianski", "French Frigate Shoals", "Kauai", "Oahu", "Molokai", "Lanai", "Maui", "Hawaii") # create order from northwest to southeast

hcbc_clust <- merge(hcbc_clust, cld, all.x = TRUE)
hcbc_clust <- hcbc_clust[order(factor(hcbc_clust$ZoneName, levels = ZoneName_ord)),] # order by zone
hcbc_clust$ZoneName <- factor(hcbc_clust$ZoneName, levels = ZoneName_ord)
hcbc_clust <- hcbc_clust[order(factor(hcbc_clust$Island_Name, levels = IslandName_ord)),] # order by island
hcbc_clust$Island_Name <- factor(hcbc_clust$Island_Name, levels = IslandName_ord,
                                 labels =  c("Kure", "PHR", "LIS", "FFS", "Kauai", "Oahu", "Molokai", "Lanai", "Maui", "Hawaii") # create order from northwest to southeast
)
hcbc_clust$ZoneLabel <- NA
for(i in c(1:nrow(hcbc_clust))){
  if(hcbc_clust$Region_Name[i] == "MHI"){
    hcbc_clust$ZoneLabel[i] <- substring(as.character(hcbc_clust$ZoneName[i]), regexpr("_", as.character(hcbc_clust$ZoneName[i])) + 1, nchar(as.character(hcbc_clust$ZoneName[i])))
  }  
}

blues_6 = c("#084594", "#3182BD", "#6BAED6", "#9ECAE1", "#C6DBEF", "#EFF3FF")

# plot raw data as a boxplot with raw data points sized by weight
zone_raw_plot_nwhi <- ggplot(hcbc_clust[ which(hcbc_clust$Region_Name == "NWHI"),], aes(x = Island_Name, y = CoralBleached_Perc_mn)) +
  geom_jitter(aes(x = Island_Name, y = CoralBleached_Perc_mn, size = Weights_DriversAndSpatialAnalysis), width = 0.2) +
  geom_boxplot(aes(fill = Island_Name), outlier.shape = NA, alpha = 0.6) +
  facet_grid(~Island_Name, scales = "free_x", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 18),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        text = element_text(size = 18),
        axis.text.y = element_text(colour="black"),
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank()
  ) +
  scale_fill_manual(values = c("#08306B","#08306B","#08306B","#08306B")) +
  xlab("") +
  ylab("% Bleached\n") +
  scale_x_discrete(labels = c("","","","")) 

zone_raw_plot_mhi <- ggplot(hcbc_clust[ which(hcbc_clust$Region_Name == "MHI"),], aes(x = ZoneName, y = CoralBleached_Perc_mn)) +
  geom_jitter(aes(x = ZoneName, y = CoralBleached_Perc_mn, size = Weights_DriversAndSpatialAnalysis), width = 0.2) +
  geom_boxplot(aes(fill = Island_Name), outlier.shape = NA, alpha = 0.6) +
  facet_grid(~Island_Name, scales = "free_x", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 18),
        legend.position = c(0.945,0.9),
        axis.line = element_line(color = "black"),
        text = element_text(size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y =element_blank()
  ) +
  scale_fill_manual(values = blues_6, guide = F) +
  scale_size_continuous(name = "Weight") +
  xlab("") +
  ylab("% Bleached\n") +
  scale_x_discrete(breaks = hcbc_clust$ZoneName, labels =  hcbc_clust$ZoneLabel) 

grid.arrange(zone_raw_plot_nwhi, zone_raw_plot_mhi, nrow = 1, widths = c(1.1,3))

# plot predictions with 95% CI as error bars and sig letters
zone_pred_plot_nwhi <- ggplot(hcbc_clust[ which(hcbc_clust$Region_Name == "NWHI"),], aes(x = Island_Name, y = predict_ble, fill = Island_Name)) +
  geom_bar(stat = "identity", position = position_dodge(), size = 0.3, color = "black") +
  facet_grid(~Island_Name, scales = "free_x", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        text = element_text(size = 8),
        axis.text.y = element_text(colour="black"),
        axis.ticks.x = element_blank()
        #axis.text.x = element_blank()
  ) +
  scale_fill_manual(values = c("#08306B","#08306B","#08306B","#08306B")) +
  xlab("") +
  ylab("Predicted % Bleached\n") +
  scale_y_continuous(limits = c(-1.1,10),breaks = c(0,2.5, 5, 7.5, 10), labels = c(0,5,25,56.25,100)) +
  geom_errorbar(aes(ymin = predict_ble - error_sp, 
                    ymax = predict_ble + error_sp),
                width=.2, size = 0.3,
                position=position_dodge(.9)) +
  scale_x_discrete(labels = c("","","","")) +
  geom_text(aes(x=Island_Name,y=error_sp+predict_ble,label=.group),vjust=-0.5, size = 2)

zone_pred_plot_mhi <- ggplot(hcbc_clust[ which(hcbc_clust$Region_Name == "MHI"),], aes(x = ZoneName, y = predict_ble, fill = Island_Name)) +
  geom_bar(stat = "identity", position = position_dodge(), size = 0.3, color = "black") +
  facet_grid(~Island_Name, scales = "free_x", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        text = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y =element_blank()
  ) +
  scale_fill_manual(values = blues_6) +
  xlab("") +
  ylab("Predicted % Bleached\n") +
  scale_y_continuous(limits = c(-1.10,10),breaks = c(0,2.5, 5, 7.5, 10), labels = c(0,5,25,56.25,100)) +
  geom_errorbar(aes(ymin = predict_ble - error_sp, 
                    ymax = predict_ble + error_sp),
                width=.2, size = 0.3,
                position=position_dodge(.9)) +
  scale_x_discrete(breaks = hcbc_clust$ZoneName, labels =  hcbc_clust$ZoneLabel) +
  geom_text(aes(x=ZoneName,y=error_sp+predict_ble,label=.group),vjust=-0.5, size = 2)

setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Drafts/For Submission - PLoS One/Figures")
tiff("Fig4.tiff", width = 2250, height = 1000, units = "px", res=300)
grid.arrange(zone_pred_plot_nwhi, zone_pred_plot_mhi, nrow = 1, widths = c(1.1,3))
dev.off()
