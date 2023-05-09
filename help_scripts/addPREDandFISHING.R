#### adding predation and fisheries data ####


#### SEAL ####

df_juv$seal = NA

for(y in 1989:2020){ # loop through all years of seal data
  
  # load seal data
  seals = raster(paste0( "//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/SealMapsConsumption/Seal consumtion_", y,  ".tif"))
  
  # subset data to year
  dat = df_juv[year(df_juv$date) == y , c("longitude", "latitude")]
  
  # create spatial points object and fix projection
  C = SpatialPoints(
    coords = cbind(dat$longitude, dat$latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  C2 = spTransform(
    C,
    
    CRSobj =  crs(seals)
  )
  
  # add data into main frame
  df_juv$seal[year(df_juv$date) == y ] = extract(seals, C2)
  
  
  # distance to available data
  df_juv$dist_seal[year(df_juv$date) == y & !is.na(df_juv$seal)] = 0
  
  
  
  # values will be NA either due to absence of seals or due to resolution of data layer
  # extract closest value within 10k of point and save distance
  
  # subset data to year and cases where values are NA
  dat = df_juv[year(df_juv$date) == y & is.na(df_juv$seal), c("longitude", "latitude")]
  
  # create spatial points object and fix projection
  C = SpatialPoints(
    coords = cbind(dat$longitude, dat$latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  C2 = spTransform(
    C,
    CRSobj =  crs(seals)
  )
  
  # identify nearest points with not-NA values and turn into spatial points object
  C3 = nearestLand(coordinates(C2), seals, 10000)
  C4 = SpatialPoints(coords = C3[!is.na(C3[ , 1]),], proj4string = crs(seals))
  
  # calculate distance between new and old points
  df_juv$dist_seal[year(df_juv$date) == y & is.na(df_juv$seal)][!is.na(C3[ , 1])] = 
    spDists(
      spTransform(C2[!is.na(C3[ , 1]),], CRSobj = "+proj=longlat +datum=WGS84") , 
      spTransform(C4, CRSobj = "+proj=longlat +datum=WGS84"), 
      longlat = TRUE,
      diagonal = TRUE)*1000
  
  
  # extract seal values where missing
  df_juv$seal[year(df_juv$date) == y & is.na(df_juv$seal)][!is.na(C3[ , 1])] = extract(seals, C4)
  
  
  
}




# since all of coast is surveyed, NA = 0
df_juv$seal[is.na(df_juv$seal)] = 0



#### CORMORANT ####

# load the two cormorant rasters
corm2006 = raster("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/Cormorant2006/cormor06pyh/hdr.adf")
corm2012 = raster("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/Cormorant2012/inv2012pyh/hdr.adf")


# since all of coast is surveyed, NA = 0
raster::values(corm2006)[is.na(raster::values(corm2006))] = 0
raster::values(corm2012)[is.na(raster::values(corm2012))] = 0


# make sure the rasters match up
raster::origin(corm2012) = raster::origin(corm2006)
corm2006_NEW = raster::crop(corm2006, extent(corm2012), snap = "out") 
corm2012_NEW = raster::crop(corm2012, extent(corm2006_NEW))


# average the two surveys
cormAVG = (corm2006_NEW + corm2012_NEW)/2


# extract data from the average raster
df_juv$cormorant = NA
C = SpatialPoints(
  coords = cbind(df_juv$longitude, df_juv$latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
C2 = spTransform(
  C,
  
  CRSobj =  crs(cormAVG)
)

df_juv$cormorant = extract(cormAVG, C2) # NA outside raster, 0 = within raster but no corms

# where data are available distance to data = 0
df_juv$dist_corm = NA
df_juv$dist_corm[!is.na(df_juv$cormorant) & df_juv$cormorant != 0] = 0


# values will be NA either due to absence of cormorants or due to resolution of data layer
# extract closest value within 10k of point and save distance

# create spatial points object of missing values
C = SpatialPoints(
  coords = cbind(df_juv$longitude[df_juv$cormorant == 0 & !is.na(df_juv$cormorant)], df_juv$latitude[df_juv$cormorant == 0 & !is.na(df_juv$cormorant)]),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

C2 = spTransform(
  C,
  CRSobj =  crs(cormAVG)
)

# identify nearest points with not-NA values and turn into spatial points object
C3 = nearestLand_0(coordinates(C2), cormAVG, 10000)
C4 = SpatialPoints(coords = C3[!is.na(C3[ , 1]),], proj4string = crs(cormAVG))


# calculate distance between new and old points
df_juv$dist_corm[df_juv$cormorant == 0 & !is.na(df_juv$cormorant)][!is.na(C3[ , 1])] = 
  spDists(
    spTransform(C2[!is.na(C3[ , 1]),], CRSobj = "+proj=longlat +datum=WGS84") , 
    spTransform(C4, CRSobj = "+proj=longlat +datum=WGS84"), 
    longlat = TRUE,
    diagonal = TRUE)*1000


# extract cormorant values where missing
df_juv$cormorant[df_juv$cormorant == 0 & !is.na(df_juv$cormorant)][!is.na(C3[ , 1])] = extract(cormAVG, C4)




# make a plot to show extent to which cormorant densities differ between surveys
samps = which( !is.na(values(corm2006_NEW) ))[ c(TRUE,rep(FALSE, 99)) ] # sample every 100th value

corm_df = 
  data.frame(
    cell = samps,
    value2006 = values(corm2006_NEW)[samps],
    value2012 = values(corm2012_NEW)[samps],
    valueMean = values(cormAVG)[samps]
  )


ggplot(corm_df, aes(x = value2006, y = value2012, colour = valueMean)) +
  geom_point() +
  scale_colour_gradientn(colours = rev(met.brewer("Greek")), name = "Average consumption") +
  labs(x = "Cormorant fish consumption 2006", y = "Cormorant fish consumption 2012") +
  theme_bw() +
  theme_sets



### alt cormorant ###

df_juv$corm2 = NA
cormDIFF = (corm2012_NEW - corm2006_NEW)

for(y in 2001:2020){
  
  print(y)
  
  sub = df_juv[df_juv$year == y, ]
    
  C = SpatialPoints(
    coords = cbind(sub$longitude, sub$latitude),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  C2 = spTransform(
    C,
    
    CRSobj =  crs(cormAVG)
  )
  
  if(y <= 2006){
  df_juv$corm2[df_juv$year == y] = extract(corm2006_NEW, C2) # NA outside raster, 0 = within raster but no corms
  
  plot(corm2006_NEW, main = y)
  }
  
  if(y >= 2012){
    df_juv$corm2[df_juv$year == y] = extract(corm2012_NEW, C2) # NA outside raster, 0 = within raster but no corms
  
    plot(corm2012_NEW, main = y)
    }
  
  if(y %in% 2007:2011){
    
    corm_int = corm2006_NEW + cormDIFF*(y-2006)/(length(2007:2012))
    
    df_juv$corm2[df_juv$year == y] = extract(corm_int, C2) # NA outside raster, 0 = within raster but no corms
    
    plot(corm_int, main = y)
  }
  
  
}

plot(df_juv$cormorant, df_juv$corm2)

#### TOTAL PREDATORS ####

df_juv$totalTopPred = df_juv$corm2 + df_juv$seal



#### FISHERIES ####

## reacreational ##

# load output from survey/population data
abbo06 = raster("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/HobbyFishery/Abborre0608/abbo0608/hdr.adf")
abbo10 = raster("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/HobbyFishery/Abborre2010/abbo10/hdr.adf")
abbo13 = raster("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/HobbyFishery/Abborre1314/abbo1314/hdr.adf")

# calculate mean of the three rasters
mean_abbo = (abbo06 + abbo10 + abbo13)/3

# make a plot to show extent to which estimated fishing pressure differ between surveys
samps = which(!is.na(values(mean_abbo)))[ c(TRUE,rep(FALSE, 99)) ]  # sample every 500th value

abbo_df = 
  data.frame(
    cell = rep(samps, 3),
    year = rep(c(2006, 2010, 2013), each = length(samps)),
    value = c(values(abbo06)[samps],
    values(abbo10)[samps],
    values(abbo13)[samps]),
    valueMean = rep(values(mean_abbo)[samps], 3)
  )


ggplot(abbo_df, aes(x = year, y = value, colour = valueMean, group = cell)) +
  geom_line() +
  scale_colour_gradientn(colours = rev(met.brewer("Greek")), name = expression(Average ~ fishing ~ pressure ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +
  labs(x = "Year", y = expression(Recreational ~ fishing ~ pressure ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +
  theme_bw() +
  theme_sets +
  theme(legend.position = "bottom")

ggsave("suppFigures/ConsistencyFisheryRecreational.jpg", width = 18.5, height = 12, units = "cm")




## commercial ##

# load commercial perch data
abboC = st_read("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/CommercialFishery/Yrkesfiske_abborre_kg_per_ICESruta_och_år_1999_2015.shp")


# rasterize each yearly polygon (1999-2015) and calculate the mean
abboCras =
  (
    fasterize(abboC, mean_abbo, field = c("X99perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X00perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X01perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X02perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X03perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X04perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X05perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X06perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X07perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X08perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X09perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X10perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X11perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X12perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X13perKM2")) +
      fasterize(abboC, mean_abbo, field = c("X14perKM2")) + 
      fasterize(abboC, mean_abbo, field = c("X15perKM2")) 
  )/length(c(1999:2015))



# make a plot to show extent to which estimated fishing pressure differ between surveys

comm_df = 
  data.frame(
    value = c(
      abboC$X99perKM2,
      abboC$X00perKM2,
      abboC$X01perKM2,
      abboC$X02perKM2,
      abboC$X03perKM2,
      abboC$X04perKM2,
      abboC$X05perKM2,
      abboC$X06perKM2,
      abboC$X07perKM2,
      abboC$X08perKM2,
      abboC$X09perKM2,
      abboC$X10perKM2,
      abboC$X11perKM2,
      abboC$X12perKM2,
      abboC$X13perKM2,
      abboC$X14perKM2,
      abboC$X15perKM2)
    )

comm_df$year = rep(1999:2015, each = nrow(comm_df)/length(1999:2015))
comm_df$loc = rep(1:sum(comm_df$year == 1999))

means = aggregate(value ~ loc, comm_df, mean); names(means) = c("loc", "mean")
comm_df = merge(comm_df, means, by = "loc")

ggplot(comm_df, aes(x = year, y = value, colour = mean, group = loc)) +
  geom_line() +
  scale_colour_gradientn(colours = rev(met.brewer("Greek")), name = expression(Average ~ fishing ~ pressure ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +
  labs(x = "Year", y = expression(Commercial ~ fishing ~ pressure ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +
  theme_bw() +
  theme_sets +
  theme(legend.position = "bottom")

ggsave("suppFigures/ConsistencyFisheryCommercial.jpg", width = 18.5, height = 12, units = "cm")







## total ##

# add up the two types of fishing pressure
totFish = mean_abbo + abboCras

# make plot
test_spdf <- as(totFish, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)

ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=layer), alpha=0.8) + 
  scale_fill_gradientn(colours = met.brewer("VanGogh3"), name = expression(Average ~ fishing ~ pressure ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +
  labs(x = " ", y = " ") +
  coord_equal() +
  theme_bw() +
  theme_sets + 
  theme(legend.position = "bottom",
        axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )



ggsave("suppFigures/FisheryPressure.jpg", width = 18.5, height = 24, units = "cm")





# extract for each bay-year 
C = SpatialPoints(
  coords = cbind(df_juv$longitude, df_juv$latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

C2 = spTransform(
  C,
  CRSobj =  crs(totFish)
)

df_juv$fishing = extract(totFish, C2)

df_juv$dist_fish = NA
df_juv$dist_fish[!is.na(df_juv$fishing)] = 0

# fix missing values (likely on land due to different raster projections)

# find nearest land within 6 km
C = SpatialPoints(
  coords = cbind(df_juv$longitude[is.na(df_juv$fishing)], df_juv$latitude[is.na(df_juv$fishing)]),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

C2 = spTransform(
  C,
  CRSobj =  crs(totFish)
)

C3 = nearestLand(coordinates(C2), totFish, 6000)
C3 = SpatialPoints(coords = C3, proj4string = crs(cormAVG))

# calculate distance between new and old points
df_juv$dist_fish[is.na(df_juv$fishing)] = 
  spDists(
    spTransform(C2, CRSobj = "+proj=longlat +datum=WGS84") , 
    spTransform(C3, CRSobj = "+proj=longlat +datum=WGS84"), 
    longlat = TRUE,
    diagonal = TRUE)*1000



# add in values for new points
df_juv$fishing[is.na(df_juv$fishing)] = extract(totFish, C3)


