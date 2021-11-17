### island-level time series of degree heating weeks (DHW) and sea surface temperature (SST) ###

#### load packages ####
library(sf)
library(maps) 
library(mapdata) 
library(maptools) 
library(rerddap) 
library(rerddapXtracto)
library(ggplot2)
library(dplyr)
library(ncdf4)
library(httr)
library(rgdal)
library(plyr)

#### import shapefile ####
## shapefile available for download here: https://www.nodc.noaa.gov/archive/arc0186/0239862/1.1/data/0-data/
# setwd("~/Morgan Winston/Projects/Satellite Course/Shapefiles")
# isl=st_read('haw_isl_6km_buffer.shp') # load island shapefiles (6km buffer around each island of interest)
# ## this shapefile has an attribute table that has the name of each island
names <- unique(isl$ISLAND_CD) # unique island codes (note, i subsetted the shapefile to only include the islands i want for my manuscript)
names <- names[c(2:10)]
years <- seq(from = 1985, to = 1993) # sequence of years i want to download data for

#### extract data ####
# datasets to download:
## SST ## 
# parameter name: SST_CRW_Daily
# sensor dataset: Coral Reef Watch
# frequency: daily
# URL: https://oceanwatch.pifsc.noaa.gov/erddap/
# dataset ID: CRW_sst_v1_0
# grid variable: analysed_sst

## DHW ##
# parameter name: Degree_Heating_Weeks
# sensor dataset: Coral Reef Watch
# frequency: daily
# URL: https://oceanwatch.pifsc.noaa.gov/erddap/
# dataset ID: CRW_dhw_v1_0
# grid variable: degree_heating_week

ERDDAP_Node="https://oceanwatch.pifsc.noaa.gov/erddap/" # set the node for data pulling 
out_df <- NULL # define dataframe that the for loops will populate
test_df <- NULL

# here is are nested for loops that cycle through each island and year extracts sst + dhw data, then merges it all together
for(i in c(1:length(names))){
  for(j in c(1:length(years))){
  print(paste("Extracting data for ", names[i], " ", years[j], "...", sep = ""))
  isl_2 <- isl[ which(isl$ISLAND_CD == names[i]),] # subset to island as specified by names[i]
  isl_2=st_coordinates(isl_2) # retrieves coordinates in matrix form for that island that we want data for

  poly <- as.data.frame(isl_2[,1:2]) # turn into a data frame
  names(poly) <- c("lon","lat") # add column names to coordinates

  # convert the longitudes to 0-360º to make it easier to work across the dateline
  I=which(poly$lon<0)
  poly$lon[I]=poly$lon[I]+360

  ## extraction begins ##
  xcoord <- poly$lon # set boundaries
  ycoord <- poly$lat
  if(years[j] == "1985"){ # set the period we pull data for, need to add time in for 1985 to not throw error
    tcoord <- c(paste(years[j],"-03-25T12:00:00", sep = ""), paste(years[j],"-12-31", sep = ""))
  }
  if(years[j] > "1985"){
    tcoord <- c(paste(years[j],"-01-01", sep = ""), paste(years[j],"-12-31", sep = ""))
  }

  # set data info for sst
  dataInfo <- rerddap::info('CRW_sst_v1_0', url=ERDDAP_Node) # define what dataset we want
  dataInfo$variable$variable_name # choose from these variables: "analysed_sst" "sea_ice_fraction"
  parameter=dataInfo$variable$variable_name[1] # select sst as parameter of interest
  sst <- rxtractogon(dataInfo, parameter=parameter, xcoord=xcoord, ycoord=ycoord, tcoord=tcoord) # extract sst data
  # image(sst$analysed_sst[,,1]) # this is how you can view the grid cell data

  # set data info for dhw
  dataInfo_2 <- rerddap::info('CRW_dhw_v1_0', url=ERDDAP_Node)
  dataInfo_2$variable$variable_name
  parameter_2=dataInfo_2$variable$variable_name[1]
  dhw <- rxtractogon(dataInfo_2, parameter=parameter_2, xcoord=xcoord, ycoord=ycoord, tcoord=tcoord) # extract dhw data
  # image(dhw$degree_heating_week[,,1])

  # calculate the spatial mean per timestep:
  ts_sst <- aaply(sst$analysed_sst, 3, mean, na.rm = TRUE)# calculate spatial average -- the '3' refers to the date column, which is what we want an average for
  ts_dhw <- aaply(dhw$degree_heating_week, 3, mean, na.rm = TRUE)

  # return a dataframe with the spatial mean per timestep:
  this_df <- data.frame(date = sst$time, year = years[j], isl = names[i], sst_mn = ts_sst, dhw_mn = ts_dhw, stringsAsFactors = FALSE)
  out_df <- rbind(out_df, this_df) # now merge all the dataframes so that we have a single data frame with the mean SST per timestep, per island
  print(paste("SST + DHW data added for ", names[i], years[j])) # print out status (it takes a while to download all this data!)
  }
}


#### plot time series figures per island ####
setwd("~/Morgan Winston/Projects/Satellite Course/Time Series Figures")
str(out_df)
summary(out_df)
saved_df <- out_df # save a copy
write.csv(saved_df, "DHW_SST_1985.2019_6km.buffer.haw.csv")
df_plot <- out_df 
df_plot <- read.csv("DHW_SST_1985.2019_6km.buffer.haw.csv")

# need to add a month + day column
df_plot$date<- as.Date(df_plot$date, format = "%Y-%m-%d")
df_plot$MM <- as.numeric(format(df_plot$date, "%m"))
df_plot$DD <- as.numeric(format(df_plot$date, "%d"))
df_plot$MM0=df_plot$MM-1
MonthDays=ddply(df_plot,.(MM0),summarize,MaxDays=max(DD))
df_plot$FracMonth=df_plot$MM0+df_plot$DD/MonthDays$MaxDays[match(df_plot$MM0,MonthDays$MM0)]
df_plot=df_plot[order(df_plot$date),]
df_plot$frame=1:nrow(df_plot)
df_plot$month=df_plot$FracMonth

# calculate maximum of monthly mean SST climatology (MMM) and bleaching threshold SST (1degC above MMM)
mean_month <- aggregate(list(sst_mn = df_plot$sst_mn), by = list(isl = df_plot$isl, month = df_plot$MM), mean)
MMM_isl <- aggregate(list(MMM = mean_month$sst_mn), by = list(isl = mean_month$isl), max)
df_plot <- merge(df_plot, MMM_isl)
df_plot$ble_thresh <- NA
for(i in c(1:nrow(df_plot))){
  df_plot$ble_thresh[i] <- df_plot$MMM[i] + 1
}

## NWHI
nwhi <- df_plot[ which(df_plot$isl %in% c("FFS", "LIS", "PHR", "KUR")),]
nwhi$isl <- droplevels(nwhi$isl)
levels(nwhi$isl) <- c("French Frigate Shoals", "Kure", "Lisianski", "Pearl & Hermes")
summary(nwhi)
# range of sst = 17.82-29.46 
# range of dhw = 0-18.280

# when was the peak DHW experienced in the NWHI?
aggregate(nwhi$dhw_mn, by = list(isl = nwhi$isl), max)
range(nwhi[which(nwhi$isl == "French Frigate Shoals" & nwhi$dhw_mn == 12.21167),]$date)

nwhi_plot <- ggplot(nwhi, aes(x = month, y = sst_mn)) +
  geom_line(aes(group = year, color = "1985-2019")) +
  geom_line(aes(group = year, color = "2019"), data=subset(nwhi,year=="2019"),size=0.9) +
  geom_line(aes(group = year, color = "2014"), data=subset(nwhi,year=="2014"),size=0.9) +
  geom_line(aes(y = dhw_mn/.8 + 4, group = year, color = "1985-2019")) +
  geom_line(aes(y = dhw_mn/.8 + 4, group = year, color = "2014"), data=subset(nwhi,year=="2014"),size=1.6) +
  geom_line(aes(y = dhw_mn/.8 + 4, group = year, color = "2019"), data=subset(nwhi,year=="2019"),size=1.6) +
  facet_wrap("isl") +
  geom_line(aes(y = MMM),color="red", size = 1) +
  geom_line(aes(y = ble_thresh), color = "red", linetype = "dashed", size = 1) +
  scale_y_continuous(name="Sea Surface Temperature (°C)\n",
                     sec.axis = sec_axis(~ . *.8 - 2, name = "Degree Heating Weeks (°C-weeks)\n"),
                     limits = c(4,30))+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,12),
                     breaks = seq(0, 11, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug","Sep", "Oct", "Nov", "Dec"),
                     name="Month") +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) + 
  scale_color_manual(values = c(
    "1985-2019" = "gray75",
    "2019" = "dodgerblue3",
    "2014" = "goldenrod1")) +
  labs(color = "Year")

ggsave(plot=nwhi_plot,
       filename="nwhi_sst_dw.png",
       height = 8,width=8)

## MHI
mhi <- df_plot[ which(df_plot$isl %in% c("HAW", "KAU", "LAN", "MAI", "MOL", "OAH")),]
mhi$isl <- droplevels(mhi$isl)
levels(mhi$isl) <- c("Hawaii", "Kauai", "Lanai", "Maui", "Molokai", "Oahu")
# range of sst = 21.88826-29.17182
# range of dhw = 0-8.700667

mhi.2 <- mhi[ which(mhi$year %in% c(2015,2019) & mhi$dhw_mn > 0),]
mhi.max.dhw <- aggregate(mhi.2$dhw_mn, by = list(yr = mhi.2$year, isl = mhi.2$isl), max)

# hawaii - 2.2
# kauai - 2.1
# lanai - 1.1
# maui - 1.2
# molokai - 2.1
# oahu - 8.2

mhi_plot <- ggplot(mhi, aes(x = month, y = sst_mn)) +
  geom_line(aes(group = year, color = "1985-2019")) +
  geom_line(aes(group = year, color = "2019"), data=subset(mhi,year=="2019"),size=0.9) +
  geom_line(aes(group = year, color = "2015"), data=subset(mhi,year=="2015"),size=0.9) +
  geom_line(aes(y = dhw_mn/.5 + 7, group = year, color = "1985-2019")) +
  geom_line(aes(y = dhw_mn/.5 + 7, group = year, color = "2015"), data=subset(mhi,year=="2015"),size=1.6) +
  geom_line(aes(y = dhw_mn/.5 + 7, group = year, color = "2019"), data=subset(mhi,year=="2019"),size=1.6) +
  facet_wrap("isl", nrow = 3) +
  geom_line(aes(y = MMM),color="red", size = 1) +
  geom_line(aes(y = ble_thresh), color = "red", linetype = "dashed", size = 1) +
  scale_y_continuous(name="Sea Surface Temperature (°C)\n",
                     sec.axis = sec_axis(~ . *.5 - 3.5, name = "Degree Heating Weeks (°C-weeks)\n"),
                     limits = c(7,30))+
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,12),
                     breaks = seq(0, 11, 1),
                     labels = c("Jan", "Feb", "Mar", "Apr","May", "Jun", "Jul", "Aug","Sep", "Oct", "Nov", "Dec"),
                     name="Month") +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  ) + 
  scale_color_manual(values = c(
    "1985-2019" = "gray75",
    "2019" = "dodgerblue3",
    "2015" = "sienna2")) +
  labs(color = "Year")

ggsave(plot=mhi_plot,
       filename="mhi_sst_dw.png",
       height = 12,width=8)


