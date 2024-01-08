# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
#### make a smooth of BIAS data to be compared w juvenile data ####
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 

# only include potential spawners
BIAS = BIAS[BIAS$length >= 5.5, ]

# translate length into weight
BIAS$weight = 5*10^(-6)*(BIAS$length*10)^3.09


# sum product of abundance and weight for all length groups within each rectangle and year
bm_rec = aggregate(abundance*weight ~ rectangle + year, data = BIAS, FUN = sum)


# extract the area of each rectangle
areas = aggregate(area ~ rectangle + year, data = BIAS, FUN = mean)


# merge biomass and area info
BIAS = merge(bm_rec, areas)
names(BIAS) = c("rectangle", "year", "biomass_ton", "area")


# translate from grams to tonnes
BIAS$biomass_ton = BIAS$biomass_ton/1000000


# translate into densities (3.4299 is constant for translating from nm2 to km2)
BIAS$biomass_density_tonnes_per_km2 = BIAS$biomass_ton/(BIAS$area*3.4299)

# merge BIAS data with spatial info on recs
BIAS = merge(BIAS, recs, by = "rectangle", all.x = T, all.y = F)


#### try to make a smooth ####

# smooth paras
max_dist = 150 # https://link.springer.com/content/pdf/10.1007/s13280-021-01684-x.pdf
min_samp = 3

# add column to store calculated BIAS abundances
df_juv$BIASmean = NA

# this is the baseline in point format
locs = SpatialPoints(
  coords = cbind(baseline$long, baseline$lat),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

# paras for weighting
sigma = max_dist/1.96

plot(1:250, dnorm(1:250, 0, sigma)/dnorm(0, 0, sigma), type = "l",
      xlab = "distance (km)", ylab = "weight")


# create a progress bar
print("Creating BIAS smooth")
pb = txtProgressBar(min = 1, 
                    max = nrow(df_juv), 
                    initial = 0,
                    style = 3) 

for(i in 1:nrow(df_juv)){ # loop through all juvenile data
  
  # progress bar
  setTxtProgressBar(pb, i)
  
  # save year of data point
  y = df_juv$year[i] 
  
  
  # subset BIAS data to this year + make into spatial points object
  subB = BIAS[BIAS$year == y-1,]
  
  if(nrow(subB) > 0){
    
    locs_B = SpatialPoints(
      coords = cbind(subB$x_mid, subB$y_mid),
      proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    
    # create spatial points object of current data row
    C = SpatialPoints(
      coords = cbind(df_juv$longitude[i], df_juv$latitude[i]),
      proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    # extract closest point on baseline
    base.pt = baseline[which.min(spDistsN1(locs, C, longlat = T)),] # this assumes straight distances 
    
    # make baseline point into spatial points object
    C = SpatialPoints( 
      coords = cbind(base.pt$long, base.pt$lat),
      proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    
    # calculate distance from baseline point to all BIAS rectangles
    subB$distance = spDistsN1(locs_B, C, longlat = T)
    
    if(sum(subB$distance < max_dist) < min_samp) df_juv$BIASmean[i] = NA
    if(sum(subB$distance < max_dist) >= min_samp) df_juv$BIASmean[i] = weighted.mean(subB$biomass_density_tonnes_per_km2, w = dnorm(subB$distance, 0, sigma))
    
  }
}

# save so that it could be loaded later
write.csv(df_juv, "data/df_BIAS.csv", row.names = FALSE)

