# load libraries
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(ggmap)
library(sf)
library(ggnewscale)
library(scales)
library(gridExtra)
library(plyr)
library(grid)

# load extra functions
g_legend<-function(a.gplot){ # function to extract legend from a plot
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# load googlemaps API Key
ggmapAPI = readChar("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Projects/Tom's Mapping Example/morgan.winston_googlemapsAPIkey.txt",nchars = 39)
register_google(key = ggmapAPI)

# set working directory
setwd("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching")
# import zone shapefile - readOGR() converts the shapefile into a Spatial Polygons Data Frame
## shapefile available for download here: https://www.nodc.noaa.gov/archive/arc0186/0239862/1.1/data/0-data/
zones = readOGR("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/Current Database/For Analysis/Shapefiles/HCBC_Zones.shp")
# import 2019 coral bleaching observations at the cluster level
## cluster level bleaching data & metadata described and linked to here: https://www.fisheries.noaa.gov/inport/item/64324
h_19 = read.csv("C:/Users/Morgan.Winston/Desktop/MHI NWHI 2019 Coral Bleaching/Data/Bleaching Assessments/Combined/Current Database/For Analysis/Analysis Ready/HCBC_2019_SpatialAnalysis.csv")

# convert the 2019 bleaching observations to a spatial points data frame for mapping
h.sp = SpatialPointsDataFrame(coords=h_19[,c("Longitude_DD_mn","Latitude_DD_mn")], data=h_19)
crs(h.sp) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
scale = 1.3 # set the scale by which to expand the map extent (generaly 1.3 works for all islands with exceptions)

### MHI MAP ####

## ## ## ## ##
#### MAUI ####
## ## ## ## ##
maui <- h.sp[ which(h.sp$Island_Name == "Maui"),] # subset data to only include clusters from Maui for mapping (except I ended up using the full dataset in most of these maps for geom_point to get the full depth bin legend)
isl_bb = scale*extent(maui) # set the extent (xmin, xmax, ymin, ymax) for the Maui base map

# convert zone shapefile to simple feature object (sf) for stable plotting in ggplot
ISL_zone <- st_as_sf(zones)
ISL_zone <- ISL_zone[ which(ISL_zone$ISLAND == "Maui"),] # only include Maui zones for this map

# download base map and check zoom level
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=10, maptype="satellite") # zoom 10 works for nui
#ggmap(isl_base) # quick look at the map

# create a dataframe that has coordinates for the labels of each of Maui's zones
zone_points <- sf::st_point_on_surface(ISL_zone)
zone_coords <- as.data.frame(sf::st_coordinates(zone_points))
zone_coords$NAME <- ISL_zone$ZoneName
zone_coords$NAME <- c("South", "West", "West\nNorthwest", "North", "Northwest", "Southeast") # i know from looking at zone_points that this is the order of zone names
zone_coords$X <- c(-156.4, -156.55, -156.72, -156.42, -156.63, -156.12) # manually set the coordinates for where the zone names are located on the map
zone_coords$Y <- c(20.5, 20.69615, 20.75, 21, 21.11, 20.58)

# run final plot
maui_map=ggmap(isl_base) + # plots base map
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones, transparent with white border
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(h.sp), size = 5) + # plots clusters (% bleaching) - color = PctBleached_mean tells R to use a continuous gradient to color by
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) + # sets the palette for coloring the cluster points, and limits = c(0,100) makes sure this is consistent between islands (have to include in all island maps)
  geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 8)+ # adds the text 
  scale_shape_manual(name = "", values = c(21,24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    title = element_text(size = 20),
    legend.position = "none"
  ) +
  # crop the map to remove extra space and manually set the breaks/labels
  coord_sf(expand = F) + # if you don't include this then gray space remains after you crop the map
  scale_x_continuous(limits = c(-156.8, -156),breaks = c(-156.7, -156.5, -156.3, -156.1), labels = c(-156.7, -156.5, -156.3, -156.1)) + 
  scale_y_continuous(limits = c(20.47,21.13), breaks = c(20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1), labels = c(20.5, 20.6, 20.7, 20.8, 20.9, 21.0, 21.1)) +
  ggtitle("Maui")

## ## ## ## ## ##
#### MOLOKAI ####
## ## ## ## ## ##
# for molokai create an inset map because there were really only surveys done along a limited section of coastline

mol <- h.sp[ which(h.sp$Island_Name == "Molokai"),]
isl_bb = 2.5*extent(mol)

# load zones shapefile, convert to SF for stable plotting in ggplot
ISL_zone <- st_as_sf(zones)
ISL_zone <- ISL_zone[ which(ISL_zone$ISLAND == "Molokai"),]

# download base map and check zoom level
isl_base_zoom = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=11, maptype="satellite")
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=10, maptype="satellite")
zone_points <- sf::st_point_on_surface(ISL_zone)
zone_coords <- as.data.frame(sf::st_coordinates(zone_points))
zone_coords$NAME <- ISL_zone$ZoneName
zone_coords$NAME <- c("Southeast", "South")
zone_coords$Y <- c(21.005, 21.035)

# create zoomed in map
mol_zoom=ggmap(isl_base_zoom) +
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(h.sp), size = 5) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(21,24,25)) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  coord_sf(expand = F) +
  geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 8)+ 
  theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank(),
    panel.border = element_rect(colour = "white", fill=NA, size=2)
  ) +
  scale_y_continuous(limits = c(20.99,21.12))

# save the zoomed in map as a grob
mol_zoom_grob <- ggplotGrob(mol_zoom)

# create overall map that nests the zoomed in map inside of it
mol_map = ggmap(isl_base) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  coord_sf(expand = F) +
  ggtitle("Molokai") +
  geom_segment(aes(x = -157.2, y = 20.99, xend = -157.4, yend = 20.94), colour="white") + # draw lines 
  geom_segment(aes(x = -156.8, y = 20.99, xend = -156.59, yend = 20.95), colour="white") +
  annotation_custom(grob = mol_zoom_grob, xmin = -157.42, xmax = -156.55, # add the zoomed in map
                    ymin = 20.68, ymax = 20.98) +
  geom_rect(aes(xmin = -157.2, xmax = -156.8, ymin = 20.99, ymax =  21.12), # draw rectangle showing extent of zoomed in map
            fill = "transparent", color = "white", size = 1) +
  scale_y_continuous(limits = c(20.69, 21.32), breaks = c(20.7, 20.8, 20.9, 21, 21.1, 21.2, 21.3), labels = c(20.7, 20.8, 20.9, 21, 21.1, 21.2, 21.3)) +
  scale_x_continuous(breaks = c(-157.4, -157.2, -157, -156.8, -156.6), labels = c(-157.4, -157.2, -157, -156.8, -156.6)) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    title=element_text(size = 20),
    axis.text = element_text(size = 20, colour = "black")
  )

## ## ## ## ##
### LANAI ####
## ## ## ## ##

lan <- h.sp[ which(h.sp$Island_Name == "Lanai"),]
isl_bb = scale*extent(lan)

# load zones shapefile, convert to SF for stable plotting in ggplot
ISL_zone <- st_as_sf(zones)

# download base map and check zoom level
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=11, maptype="satellite")
zone_points <- sf::st_point_on_surface(ISL_zone)
zone_points <- sf::st_point_on_surface(ISL_zone[ which(ISL_zone$ISLAND == "Lanai"),])
zone_coords <- as.data.frame(sf::st_coordinates(zone_points))
zone_coords$NAME <- c("South", "Northeast")
zone_coords$Y <- c(20.675, 20.95)

# run final plot
lan_map=ggmap(isl_base) +
  geom_sf(data = ISL_zone[ which(ISL_zone$ISLAND == "Lanai"),],fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(h.sp), size = 6) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(21,24,25)) +
  coord_sf(expand = F) +
  geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 8)+ 
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    legend.position = "none",
    title=element_text(size = 20)
  ) +
  scale_x_continuous(limits = c(-157.1, -156.75), breaks = c(-157, -156.9, -156.8), labels = c(-157, -156.9, -156.8)) +
  scale_y_continuous(limits = c(20.665, 20.97), breaks = c(20.7, 20.75, 20.8, 20.85, 20.9, 20.95), labels = c(20.7, 20.75, 20.8, 20.85, 20.9, 20.95)) +
  ggtitle("Lanai")

## ## ## ## ##
## HAWAII ####
## ## ## ## ##
## since Hawaii is so big and only want to show data from three zones, create inset maps for each zone

haw <- h.sp[ which(h.sp$Island_Name == "Hawaii"),]
isl_bb = scale*extent(haw)

# load zones shapefile, convert to SF for stable plotting in ggplot
ISL_zone <- st_as_sf(zones)

# download base map and check zoom level
ISL_zone <- st_as_sf(crop(x=gBuffer(zones,byid = T,width=0),y=isl_bb))
isl_base = get_map(location = c(-156, mean(isl_bb[3:4])), zoom=7, maptype="satellite")
ggmap(isl_base)

# create overall island map
hawaii=ggmap(isl_base) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  coord_sf(expand = F) +
  # geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 7)+ # labels, remove for NW and SW Hawaii
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 18, colour = "black"),
    legend.position = "none"
  ) +
  scale_x_discrete(breaks = c(-159, -158, -157, -156, -155, -154, -153), labels = c(-159, -158, -157, -156, -155, -154, -153)) +
  scale_y_discrete(breaks = c(17, 18, 19, 20, 21, 22), labels = c(17, 18, 19, 20, 21, 22))

## HAWAII NW ###
# create inset map for NW zone
haw_nw <- h.sp[ which(h.sp$ZoneName == "Hawaii_NW"),]
isl_bb = scale*extent(haw_nw)

ISL_zone <- st_as_sf(zones[ which(zones$ZoneName == "Hawaii_NW"),])
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=9, maptype="satellite")
ggmap(isl_base)

hawaii_nw=ggmap(isl_base) +
  geom_sf(data = ISL_zone[ which(ISL_zone$ZoneName == "Hawaii_NW"),],fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(haw_nw), size = 6) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-156.1,-155.82), breaks = seq(-156.1, -155.8, len = 4)) +
  scale_y_continuous(limits = c(19.46, 20.18))

## HAWAII SW ###
# create inset map for SW zone
haw_sw <- h.sp[ which(h.sp$ZoneName == "Hawaii_SW"),]
isl_bb = scale*extent(haw_sw)

ISL_zone <- st_as_sf(zones[ which(zones$ZoneName == "Hawaii_SW"),])
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=9, maptype="satellite")
ggmap(isl_base)

hawaii_sw=ggmap(isl_base) +
  geom_sf(data = ISL_zone[which(ISL_zone$ZoneName == "Hawaii_SW"),],fill=NA,color="white",alpha=.5, inherit.aes = FALSE) +
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(haw_sw), size = 6) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-156.035,-155.85)) +
  scale_y_continuous(limits = c(19, 19.51)) 

## HAWAII E ###
# create inset map for E zone
haw_e <- h.sp[ which(h.sp$ZoneName == "Hawaii_E"),]
isl_bb = scale*extent(haw_e)

ISL_zone <- st_as_sf(zones[ which(zones$ZoneName == "Hawaii_E"),])
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=9, maptype="satellite")
ggmap(isl_base)

hawaii_e=ggmap(isl_base) +
  geom_sf(data = ISL_zone[which(ISL_zone$ZoneName == "Hawaii_E"),],fill=NA,color="white",alpha=.5, inherit.aes = FALSE) +
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(haw_e), size = 7) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    legend.position = "none"
  ) +
  scale_x_continuous(limits = c(-155.12,-154.96)) +
  scale_y_continuous(limits = c(19.72, 19.82)) 

## COMBINED HAWAII ###
# now start to combine the overall map with the inset maps
remove <- theme(axis.line=element_blank(),axis.text.x=element_blank(),
                axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),legend.position="none",
                panel.background=element_blank(),panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),plot.background=element_blank(),
                panel.border = element_rect(colour = "white", fill=NA, size=2))
h_nw <- hawaii_nw + remove + ggtitle("  Northwest") + theme(plot.title = element_text(size = 18, colour = "white", margin = margin(t = 4, b = -20)))
h_sw <- hawaii_sw + remove + ggtitle("  Southwest") + theme(plot.title = element_text(size = 18, colour = "white", margin = margin(t = 4, b = -20)))
h_e <- hawaii_e + remove + ggtitle("  East") + theme(plot.title = element_text(size = 18, colour = "white", margin = margin(t = 4, b = -20)))
nw <- ggplotGrob(h_nw)
sw <- ggplotGrob(h_sw)
e <- ggplotGrob(h_e)

# create final map
hawaii_map <- hawaii + 
  geom_segment(aes(x = -156.14, y = 19.42, xend = -157, yend = 16.75), colour="white") +
  annotation_custom(grob = nw, xmin = -159.8, xmax = -156.5,
                    ymin = 16.5, ymax = 22.7) +
  geom_rect(aes(xmin = -156.14, xmax = -155.79, ymin = 19.42, ymax =  20.25), # nw
            fill = "transparent", color = "white", size = 1) +
  
  geom_segment(aes(x = -155.83, y = 18.96, xend = -154.3, yend = 16.9), colour="white") +
  annotation_custom(grob = sw, xmin = -154.6, xmax = -152.6,
                    ymin = 16, ymax = 22) +
  geom_rect(aes(xmin = -156.035, xmax = -155.83, ymin = 18.96, ymax =  19.54), # sw
            fill = "transparent", color = "white", size = 1) +
  
  geom_segment(aes(x = -154.96, y = 19.82, xend = -154.86, yend = 22), colour="white") +
  annotation_custom(grob = e, xmin = -157.15, xmax = -154.35,
                    ymin = 21, ymax = 22.5) +
  geom_rect(aes(xmin = -155.12, xmax = -154.96, ymin = 19.72, ymax =  19.82), # east
            fill = "transparent", color = "white", size = 0.5) +
  ggtitle("Hawaii") + theme(title = element_text(size = 20))

## ## ## ## ##
### KAUAI ####
## ## ## ## ##
# for kauai we also create inset maps of north and south zones

kau <- h.sp[ which(h.sp$Island_Name == "Kauai"),]
isl_bb = scale*extent(kau)

ISL_zone <- st_as_sf(crop(x=gBuffer(zones,byid = T,width=0),y=isl_bb))
isl_base_south = get_map(location = c(mean(isl_bb[1:2]), 21.85), zoom=12, maptype="satellite")
isl_base_north = get_map(location = c(mean(isl_bb[1:2]), 22.25), zoom=12, maptype="satellite")
isl_base = get_map(location = c(-159.67, mean(isl_bb[3:4])), zoom=10, maptype="satellite")

zone_points <- sf::st_point_on_surface(ISL_zone)
zone_coords <- as.data.frame(sf::st_coordinates(zone_points))
zone_coords$NAME <- ISL_zone$ZoneName

zone_coords$NAME <- c("South", "North")
zone_coords$Y <- c(21.83, 22.28)

# RASTER PLOTS
# run final plot
k_n = ggmap(isl_base_north) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(h.sp), size = 7)+
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(23,24,25)) +
  coord_sf(expand = F) +
  geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 8) +
  theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank(),
    panel.border = element_rect(colour = "white", fill=NA, size=2)
  )+
  scale_y_continuous(limits = c(22.19, 22.29))+
  scale_x_continuous(limits = c(-159.58, -159.42))

k_s = ggmap(isl_base_south) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(h.sp), size = 7) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(21,24,25)) +
  coord_sf(expand = F) +
  geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 8) +
  theme(
    axis.line=element_blank(),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position="none",
    panel.background=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),plot.background=element_blank(),
    panel.border = element_rect(colour = "white", fill=NA, size=2)
  )+
  scale_y_continuous(limits = c(21.82, 21.9)) +
  scale_x_continuous(limits = c(-159.59, -159.41))

k_north <- ggplotGrob(k_n)
k_south <- ggplotGrob(k_s)

# here we create the final map
kauai_map=ggmap(isl_base) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  coord_sf(expand = F) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    legend.position = "none",
    title = element_text(size = 20)
  ) +
  ggtitle("Kauai") +
  geom_segment(aes(x = -159.59, y = 22.28, xend = -159.63, yend = 22.44), colour="white") +
  annotation_custom(grob = k_north, xmin = -160.1, xmax = -159.61,
                    ymin = 22.075, ymax = 22.52) +
  geom_rect(aes(xmin = -159.59, xmax = -159.41, ymin = 22.19, ymax =  22.29),
            fill = "transparent", color = "white", size = 1) +
  
  geom_segment(aes(x = -159.59, y = 21.82, xend = -159.63, yend = 21.78), colour="white") +
  annotation_custom(grob = k_south, xmin = -160.1, xmax = -159.61,
                    ymin = 21.75, ymax = 22) +
  geom_rect(aes(xmin = -159.59, xmax = -159.41, ymin = 21.82, ymax =  21.9),
            fill = "transparent", color = "white", size = 1) +
  scale_x_continuous(breaks = c(-160, -159.8, -159.6, -159.4), labels = c(-160, -159.8, -159.6, -159.4)) +
  scale_y_continuous(limits = c(21.76, 22.46), breaks = c(21.8, 21.9, 22, 22.1, 22.2, 22.3, 22.4), labels = c(21.8, 21.9, 22, 22.1, 22.2, 22.3, 22.4))

## ## ## ## ##
#### OAHU ####
## ## ## ## ##
oah <- h.sp[ which(h.sp$Island_Name == "Oahu"),]
isl_bb = scale*extent(oah)

# load zones shapefile, convert to SF for stable plotting in ggplot
ISL_zone <- st_as_sf(zones)
ISL_zone <- ISL_zone[ which(ISL_zone$ISLAND == "Oahu"),] # only include Maui zones for this map

# download base map and check zoom level
isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=10, maptype="satellite") # zoom 9 for hawaii, 10 for oahu and maui nui
zone_points <- sf::st_point_on_surface(ISL_zone)
zone_coords <- as.data.frame(sf::st_coordinates(zone_points))
zone_coords$NAME <- ISL_zone$ZoneName
zone_coords$NAME <- c("South","West","East","North")
zone_coords$X <- c(-157.9699, -158.290, -157.63, -157.9038)
zone_coords$Y <- c(21.227, 21.41524, 21.44492, 21.77)

# RASTER PLOTS
# run final plot
oahu_map=ggmap(isl_base) +
  geom_sf(data = ISL_zone,fill=NA,color="white",alpha=.5, inherit.aes = FALSE) + # plots zones
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(h.sp), size = 5) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(21,24,25)) +
  coord_sf(expand = F) +
  geom_text(data = zone_coords, aes(X, Y, label = NAME), colour = "white", size = 8) + # labels, remove for NW and SW Hawaii
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 20, colour = "black"),
    legend.position = "none",
    title = element_text(size = 20),
    legend.key.width = unit(1, "inch") 
  ) +
  ggtitle("Oahu") +
  scale_x_continuous(limits = c(-158.33, -157.59), breaks = c(-158.3, -158.1, -157.9, -157.7), labels = c(-158.3, -158.1, -157.9, -157.7)) +
  scale_y_continuous(limits = c(21.18, 21.82), breaks = c(21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8), labels = c(21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8))

## > ALL MHI COMBINED < ####
## extract legend -- this is just a map that we create for solely the purpose of generating a legend that is then extracted 
legend_creator=ggmap(isl_base) +
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, fill = CoralBleached_Perc_mn, shape = Depth_bin), data = as.data.frame(oah), size = 8) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(21,24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 22, colour = "black"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 22),
    title = element_text(size = 22),
    legend.key.width = unit(1, "inch") 
  ) 

mylegend<-g_legend(legend_creator) # run the function to get legend from plot

## create the multi-panel plot
# 2 x 2 layout
grid.arrange(arrangeGrob(kauai_map, oahu_map, mol_map, maui_map, lan_map, hawaii_map, nrow = 3, widths = c(3,3), left = textGrob("Latitude",gp=gpar(fontsize=22), rot = 90), bottom = textGrob("Longitude",gp=gpar(fontsize=22))), mylegend, nrow = 2, heights = c(9.5,0.5))

## ## ## ## ## ##
### NWHI MAP ####
## ## ## ## ## ##

## KURE ####
kur <- h.sp[ which(h.sp$Island_Name == "Kure"),]
isl_bb = scale*extent(kur)

isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=12, maptype="satellite")
# ggmap(isl_base)

kure.p=ggmap(isl_base) + # plots base map
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = CoralBleached_Perc_mn), data = as.data.frame(kur), size = 4) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(23,24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none"
  ) +
  ggtitle("Kure") +
  scale_x_continuous(limits = c(-178.41,-178.25)) +
  scale_y_continuous(limits = c(28.35, 28.49))


## PHR ####
phr <- h.sp[ which(h.sp$Island_Name == "Pearl and Hermes"),]
isl_bb = scale*extent(phr)

isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=11, maptype="satellite")
ggmap(isl_base)

phr.p=ggmap(isl_base) + # plots base map
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = CoralBleached_Perc_mn), data = as.data.frame(phr), size = 4) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none"
  ) +
  ggtitle("Pearl and Hermes Atoll") +
  scale_x_continuous(limits = c(-176.05,-175.7)) +
  scale_y_continuous(limits = c(27.73, 28), breaks = c(27.75, 27.8, 27.85, 27.9, 27.95, 28), labels = c(27.75, 27.8, 27.85, 27.9, 27.95, 28))

## LIS ####
lis <- h.sp[ which(h.sp$Island_Name == "Lisianski"),]
isl_bb = scale*extent(lis)

isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=11, maptype="satellite")
ggmap(isl_base)

lis.p=ggmap(isl_base) + # plots base map
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = CoralBleached_Perc_mn), data = as.data.frame(lis), size = 4) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(23,24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none"
  ) +
  ggtitle("Lisianski") +
  scale_x_continuous(limits = c(-174.09,-173.8)) +
  scale_y_continuous(limits = c(25.9, 26.2))

## FFS ####
ffs <- h.sp[ which(h.sp$Island_Name == "French Frigate Shoals"),]
isl_bb = scale*extent(ffs)

isl_base = get_map(location = c(mean(isl_bb[1:2]), mean(isl_bb[3:4])), zoom=11, maptype="satellite")
ggmap(isl_base)

ffs.p=ggmap(isl_base) + # plots base map
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = CoralBleached_Perc_mn), data = as.data.frame(ffs), size = 4) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "none"
  ) +
  ggtitle("French Frigate Shoals") + 
  scale_y_continuous(limits = c(23.7,23.94)) +
  scale_x_continuous(limits = c(-166.4, -166.08))

## > ALL NWHI COMBINED < ####
## extract legend
map_leg=ggmap(isl_base) + # plots base map
  geom_point(aes(x = Longitude_DD_mn, y = Latitude_DD_mn, shape = Depth_bin, fill = CoralBleached_Perc_mn), data = as.data.frame(kur), size = 4) +
  scale_fill_distiller(palette = "Spectral", name = "Percent of Live Coral Bleached (%)", limits = c(0,100)) +
  scale_shape_manual(name = "", values = c(23,24,25)) +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  scale_y_continuous(limits = c(23.7,23.94)) +
  scale_x_continuous(limits = c(-166.4, -166.08))

mylegend<-g_legend(map_leg)

# 2 x 2 layout
grid.arrange(arrangeGrob(kure.p, phr.p, lis.p, ffs.p, nrow = 2, left = "Latitude", bottom = "Longitude"), mylegend, nrow = 2, heights = c(8,1))
