#### hawaii coral bleaching analysis: taxonomic variability ####
## written by: morgan winston

## this script uses survey-level data described & linked to here: {InPort record}
## code performs the following:
### i) adds depth bins, and transforms the data from long to wide
### ii) explores taxa-level differences in bleaching during the 2019 marine heatwave event in hawaii using nmds plots & basic bar graphs

#### initialization #### 

# load packages & useful functions
library(dplyr)
library(tidyr) # to run complete() which adds missing combinations of data to a data frame
library(plyr) # to use ddply for applying function to each subset of a data frame
library(ggplot2) # for creating plots
library(vegan) # for running metaMDS(), to create nmds plots
library(ade4)
library(gclus)
library(FD)
library(plotrix) # for std.error()
library(stringr) # for working with strings using function word()
library(gridExtra) # to create plots with multiple panels
library(ggpubr)
'%!in%' <- function(x,y)!('%in%'(x,y)) # the opposite of %in%
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)} # function to extract legend from plot

# set working directory: setwd("")
# import data: hcbc <- read.csv("") # data accessible at: 
setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/For InPort/Version III")
hcbc <- read.csv("HCBC_2019_Observations_Update_3.9.21.csv")


#### prep data ####

# rename columns 
colnames(hcbc)[16] <- "LiveCoralCover_Perc"
colnames(hcbc)[17] <- "LiveCoralCover_Bin"
colnames(hcbc)[18] <- "CoralBleached_Perc"
colnames(hcbc)[19] <- "CoralBleached_Bin"

# add median depth to observations that only have bins 
unique(hcbc[ which(is.na(hcbc$Depth_ft)),]$Depth_bin) # mid and shallow
for(i in c(1:nrow(hcbc))){
  if(is.na(hcbc$Depth_ft[i])){
    if(hcbc$Depth_bin[i] == "Mid"){ # mid: 19.8-59'
      hcbc$Depth_ft[i] <- 39
    }
    if(hcbc$Depth_bin[i] == "Shallow"){ ## shallow: 0-19.7'
      hcbc$Depth_ft[i] <- 10
    }
  }
}

# add depth bins to all
hcbc$Depth_bin <- as.character(hcbc$Depth_bin)
for(i in c(1:nrow(hcbc))){
  if(hcbc$Depth_bin[i] == ""){
    if(hcbc$Depth_ft[i] <= 19.7){
      hcbc$Depth_bin[i] <- "Shallow"
    }
    if(hcbc$Depth_ft[i] > 19.7 & hcbc$Depth_ft[i] <=59){
      hcbc$Depth_bin[i] <- "Mid"
    }
    if(hcbc$Depth_ft[i] > 59){
      hcbc$Depth_bin[i] <- "Deep"
    }
  }
}

# add median bin value of bleaching and live coral cover to observations missing measurement (i.e. only recorded bins)
# read in supplemental data: "BinConversion.csv" found in data > supplemental within repo 
# BinLU=read.csv("BinConversion.csv") - this dataframe gives min, max, and median of percentages per bin level 
hcbc[ which(is.na(hcbc$CoralBleached_Perc)),]$CoralBleached_Perc <- BinLU$Median[match(hcbc[ which(is.na(hcbc$CoralBleached_Perc)),]$CoralBleached_Bin, BinLU$Bin)] # bleaching
hcbc[ which(is.na(hcbc$LiveCoralCover_Perc)),]$LiveCoralCover_Perc <- BinLU$Median[match(hcbc[ which(is.na(hcbc$LiveCoralCover_Perc)),]$LiveCoralCover_Bin, BinLU$Bin)] # live coral cover

# remove columns not needed in analysis
hcbc$Comments_Observation <- NULL
hcbc$AveSeverity <- NULL
hcbc$Sector_Name <- NULL
hcbc$Survey_ID_InternalToYourProject <- NULL

# remove DAR AIS and ESD PQ from this analysis -- we  use rapid visual observations from same sites instead
hcbc <- hcbc[ which(hcbc$Dataset_Name != "ESD_PQ"),] 
hcbc <- hcbc[ which(hcbc$Dataset_Name != "DAR_Oahu_PQ"),]

# transform from wide to long format
names(hcbc)[19:58]<-c("Taxon1.Code", "Taxon1.Cover", "Taxon1.Bleached", "Taxon1.AveSev", "Taxon1.MaxSev", 
                      "Taxon2.Code", "Taxon2.Cover", "Taxon2.Bleached", "Taxon2.AveSev", "Taxon2.MaxSev",
                      "Taxon3.Code", "Taxon3.Cover", "Taxon3.Bleached", "Taxon3.AveSev", "Taxon3.MaxSev",
                      "Taxon4.Code", "Taxon4.Cover", "Taxon4.Bleached", "Taxon4.AveSev", "Taxon4.MaxSev",
                      "Taxon5.Code", "Taxon5.Cover", "Taxon5.Bleached", "Taxon5.AveSev", "Taxon5.MaxSev",
                      "Taxon6.Code", "Taxon6.Cover", "Taxon6.Bleached", "Taxon6.AveSev", "Taxon6.MaxSev",
                      "Taxon7.Code", "Taxon7.Cover", "Taxon7.Bleached", "Taxon7.AveSev", "Taxon7.MaxSev",
                      "Taxon8.Code", "Taxon8.Cover", "Taxon8.Bleached", "Taxon8.AveSev", "Taxon8.MaxSev") # rename columns so they work with reshape()

hcbc_taxa <- hcbc[c(7,19:58)] # subset to just taxa level metrics (the data that needs to be taken from wide to long)
hcbc_taxa <- reshape(hcbc_taxa, direction = 'long', # turn from wide to long format using reshape()
                            varying=list(grep("Code", colnames(hcbc_taxa), value=T), grep("Cover", colnames(hcbc_taxa), value=T), grep("Bleached", colnames(hcbc_taxa), value=T),
                                         grep("AveSev", colnames(hcbc_taxa), value=T), grep("MaxSev", colnames(hcbc_taxa), value=T)),
                            timevar="Taxon_Rank",
                            times=c("Taxon1", "Taxon2", "Taxon3", "Taxon4", "Taxon5", "Taxon6", "Taxon7", "Taxon8"),
                            v.names=c("TaxonCode", "Taxon_CoralCover_Perc", "Taxon_CoralBleached_Perc", "Taxon_AveSev", "Taxon_MaxSev"),
                            idvar="Survey_Name")

hcbc_taxa <- merge(hcbc_taxa, hcbc[c(1:18)]) # add survey level data back in
hcbc_taxa <- hcbc_taxa[ which(hcbc_taxa$TaxonCode != ""),] # remove rows where taxoncode is empty because conversion from wide to long creates extra rows if not all 8 dominant taxa recorded per site
hcbc_taxa <- hcbc_taxa[ which(!is.na(hcbc_taxa$Taxon_CoralBleached_Perc)),] # only use taxa-level data where % bleaching was recorded

# add columns for full species and genera names
# read in supplemental data: "Genus_lookup.csv" found in data > supplemental within repo 
# lookup <- read.csv("Genus_lookup.csv") # look-up table for taxa code to full name
colnames(lookup)[2] <- "TaxonCode"
hcbc_taxa <- merge(hcbc_taxa, lookup, by = "TaxonCode", all.x = TRUE)
hcbc_taxa$GENUSNAME <- paste(word(hcbc_taxa$TAXONNAME, 1), "spp.", sep = " ")
hcbc_taxa$TAXON_N <- paste(substr(hcbc_taxa$TAXONNAME, 1, 1), ". ", word(hcbc_taxa$TAXONNAME, 2), sep = "")
colnames(hcbc_taxa)

# add in site susceptibility data
s <- read.csv("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/Current Database/For Analysis/Supplemental Data/HCBC_2019_SiteSusceptibilityScores.csv")
hcbc_taxa$SiteSuscp <- s$SiteSuscep[match(hcbc_taxa$Survey_Name, s$Survey_Name)] # add scores into working data
head(hcbc_taxa)

# calculate the relative cover of each taxa at each site
hcbc_taxa$Taxon_CoralCover_Perc_rel <- (hcbc_taxa$Taxon_CoralCover_Perc*100)/hcbc_taxa$LiveCoralCover_Perc
hcbc_taxa$Taxon_CoralCover_Prop_rel <- hcbc_taxa$Taxon_CoralCover_Perc_rel/100


#### NMDS PLOTS ####
# based on previous attempts, site HAW-TNC-5C is messing up the NMDS
hcbc2 <- hcbc_taxa[ which(hcbc_taxa$Survey_Name != "HAW-TNC-5C"),]

# taxa matrix set-up: need to transform the data frame - long to wide at the survey level
taxa.livecov <- reshape(hcbc2[,c(1,2,31)], idvar = "Survey_Name", timevar = "TaxonCode", direction = "wide") ## we want to create the matrix with relative cover of each coral taxa
names(taxa.livecov) <- gsub("Taxon_CoralCover_Prop_rel.", "", names(taxa.livecov))
taxa.livecov[is.na(taxa.livecov)] <- 0

# cluster info table - island, survey, survey level suscp, survey level % bleaching
cluster.info <- unique(hcbc2[,c(2,12,22,29)])

# convert to matrices
taxa.livecov_m <- as.matrix(taxa.livecov[,-1])

# heatmaps
heatmap(taxa.livecov_m) 
heatmap(sqrt(taxa.livecov_m)) # should use sqrt transformed data

# visualize difference between islands
cov_NMDS <- metaMDS(sqrt(taxa.livecov_m), autotransform = FALSE, distance = "bray", k = 2, maxit = 999, trymax = 200) # stress = 0.15
stressplot(cov_NMDS) 

# for plotting in ggplot, must extract the scores (x and y coords of each cluster (row) and species and add the grouping variables (island, % bleaching)
data.scores <- as.data.frame(scores(cov_NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$isl <- cluster.info$Island_Name  #  add the island variable
data.scores$survey <- cluster.info$Survey_Name
data.scores$ble <- sqrt(cluster.info$CoralBleached_Perc) #  add the bleaching per cluster variable
data.scores$sus <- log(cluster.info$SiteSuscp) # add suscp
head(data.scores)  # look at the data

species.scores <- as.data.frame(scores(cov_NMDS, "species"))  # using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  # look at the data

# color by suscp plot
suscp <- ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=sus,shape=isl),size=4) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  scale_color_viridis_c(name = "b) Taxonomic Susceptibility Score", limits = c(-0.01, log(4)), breaks = c(log(1), log(2), log(3), log(4)), labels = c(1,2,3,4)) +
  scale_fill_discrete(name = "Island") +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "top", legend.box="vertical", legend.margin=margin(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))+
  scale_shape_manual(name = "Island", values = c(15,16,17,18), guide = F) +
  guides(color=guide_colourbar(barwidth=10))


# color by bleaching plot
ble <- ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=ble,shape=isl),size=4) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  scale_color_viridis_c(name = "a) % of Live Coral Cover Bleached", breaks = c(0,2.5,5,7.5,10), labels = c(0,(2.5^2),(5^2),(7.5^2),(10^2))) +
  scale_fill_discrete(name = "Island") +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.position = "top", legend.box="vertical", legend.margin=margin(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))+
  scale_shape_manual(name = "Island", values = c(15,16,17,18), guide = F) +
  guides(color=guide_colourbar(barwidth=11))


shapes <- ggplot() + # just make this plot to get the legend
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=isl),size=4) + # add the point markers
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +
  scale_shape_manual(name = "Island", values = c(15,16,17,18)) +
  theme(
    legend.box="vertical", legend.margin=margin(),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

shape_leg <- g_legend(shapes)

setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Taxa")
png(width = 1150, height = 950, filename = "NMDS_plots_new.png")
grid.arrange(ble, arrangeGrob(suscp, shape_leg,nrow = 1, widths = c(10,2)), nrow = 1, widths = c(2,2.4))
dev.off()


#### BAR PLOTS ####
# show mean % bleaching per species per island, colored by mean absolute % cover of the species per island 
sp_isl_count <- aggregate(hcbc_taxa$Survey_Name, by = list(isl = hcbc_taxa$Island_Name, taxon = hcbc_taxa$TaxonCode), function(x) length(unique(x)))
sp_isl_count <- sp_isl_count[ which(sp_isl_count$x < 6),] # species must occur within at least 6 different surveys per island to be included (in that island and also domain estimates)
drop_sp_isl <- paste(sp_isl_count$isl, sp_isl_count$taxon, sep = "_")
hcbc_taxa$drop <- paste(hcbc_taxa$Island_Name, hcbc_taxa$TaxonCode, sep = "_")
hcbc_taxa <- hcbc_taxa[ which(hcbc_taxa$drop %!in% drop_sp_isl),]
hcbc_taxa$drop <- NULL

# domain level taxa data for plotting
isl_sp <- ddply(hcbc_taxa,.(GENUS_CODE, GENUSNAME, TaxonCode, TAXONNAME),summarize, # calculate mean and std er of bleaching for each coral taxa per island
                TaxonPctBleached_se=std.error(Taxon_CoralBleached_Perc),
                TaxonPctBleached_mean=mean(Taxon_CoralBleached_Perc),
                TaxonLiveCoralCover_se_rel=std.error(Taxon_CoralCover_Perc_rel),
                TaxonLiveCoralCover_mean_rel=mean(Taxon_CoralCover_Perc_rel),
                TaxonLiveCoralCover_se_abs=std.error(Taxon_CoralCover_Perc),
                TaxonLiveCoralCover_mean_abs=mean(Taxon_CoralCover_Perc))
# order from most to least susceptible
taxa_ble <- ddply(hcbc_taxa,.(TAXONNAME),summarize, # calculate mean and std er of bleaching for each coral taxa overall
                  TaxonPctBleached_mean=mean(Taxon_CoralBleached_Perc))
taxa_ble <- taxa_ble[order(taxa_ble$TaxonPctBleached_mean),]
sp_ord <- taxa_ble$TAXONNAME # create the order
isl_sp$TAXONNAME <- factor(isl_sp$TAXONNAME, levels = sp_ord)
isl_sp <- isl_sp[order(factor(isl_sp$TAXONNAME, levels = sp_ord)),]

# island level taxa  data for plotting
hcbc_sp=ddply(hcbc_taxa,.(Island_Name, GENUS_CODE, GENUSNAME, TaxonCode, TAXONNAME),summarize, # calculate mean and std er of bleaching for each coral taxa per island
              TaxonPctBleached_se=std.error(Taxon_CoralBleached_Perc),
              TaxonPctBleached_mean=mean(Taxon_CoralBleached_Perc),
              TaxonLiveCoralCover_se_rel=std.error(Taxon_CoralCover_Perc_rel),
              TaxonLiveCoralCover_mean_rel=mean(Taxon_CoralCover_Perc_rel),
              TaxonLiveCoralCover_se_abs=std.error(Taxon_CoralCover_Perc),
              TaxonLiveCoralCover_mean_abs=mean(Taxon_CoralCover_Perc))
# put in same order as domain level
hcbc_sp$TAXONNAME <- factor(hcbc_sp$TAXONNAME, levels = sp_ord)
hcbc_sp <- hcbc_sp[order(factor(hcbc_sp$TAXONNAME, levels = sp_ord)),]
tax_names <- hcbc_sp$TAXONNAME
isl <- data.frame("TAXONNAME" = rep(unique(tax_names), 4), # adds a row for every taxa per island for plotting purposes
                  "Island_Name" = rep(c("Hawaii", "Lanai", "Oahu", "Maui"), each = length(unique(tax_names))))
hcbc_sp.2<- merge(hcbc_sp, isl, all = TRUE)

# create plots
## domain level
dom.tax <- ggplot(isl_sp, aes(x = TAXONNAME, y = TaxonPctBleached_mean, fill =  TaxonLiveCoralCover_mean_abs)) +
  geom_bar(stat = 'identity', position = position_dodge(), color = "black") +
  geom_errorbar(aes(ymin=TaxonPctBleached_mean-TaxonPctBleached_se, ymax=TaxonPctBleached_mean+TaxonPctBleached_se),  color = "black", width=.2,
                position=position_dodge(.9)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0, "lines"), 
    panel.background = element_rect(colour = "black", fill = "white"),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.y = element_text(hjust = 1, face = "italic"),
    legend.position = c(0.8,0.1),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 15),
    text = element_text(size = 18),
    title = element_text(size = 15)
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(-1,105)) +
  scale_fill_viridis_c(name = "Absolute Live Coral Cover (%)", guide = guide_colourbar(title.position = "top", barwidth = 10)) +
  coord_flip() +
  ggtitle("a)") 

## island level
isl.tax <- ggplot(hcbc_sp.2, aes(x = TAXONNAME, y = TaxonPctBleached_mean, fill =  TaxonLiveCoralCover_mean_abs)) +
  geom_bar(stat = 'identity', position = position_dodge(), color = "black") +
  facet_wrap(~ Island_Name, scales = 'free_y') +
  geom_errorbar(aes(ymin=TaxonPctBleached_mean-TaxonPctBleached_se, ymax=TaxonPctBleached_mean+TaxonPctBleached_se),  color = "black", width=.2,
                position=position_dodge(.9)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0, "lines"), 
    panel.background = element_rect(colour = "black", fill = "white"),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text.y = element_text(hjust = 1, face = "italic"),
    legend.position = "none",
    #legend.title = element_blank(),
    #legend.text = element_text(face = "italic"),
    text = element_text(size = 18),
    title = element_text(size = 15)
  ) +
  xlab("") +
  ylab("") +
  scale_y_continuous(expand = c(0,0), limits = c(-1,105)) +
  scale_fill_viridis_c(name = "Absolute Live Coral Cover (%)") +
  coord_flip() +
  ggtitle("b)")

setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/2019 Manuscript/Figures/Taxa")
png(width = 1500, height = 650, filename = "Taxa_Bar.png")
grid.arrange(dom.tax, isl.tax, nrow = 1, bottom = text_grob("% of Live Coral Cover Bleached", size = 18))
dev.off()