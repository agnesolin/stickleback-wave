# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### script adds wave exposure and distance to baseline for each detonation point ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### add SWM ####

# create spatial points object from detonation coordinates
C = SpatialPoints(
  coords = cbind(df_juv$longitude, df_juv$latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
C2 = spTransform(
  C,
  CRSobj =  crs(SWM)
)

# extract values from SWM
df_juv$SWM = extract(SWM, C2)


#### missing SWM values ####

# in some cases, SWM values are missing, e.g. if point coordinates are a little bit off
# in this case, extract SWM values from closest point


# create spatial points object of missing values
C = SpatialPoints(
  coords = cbind(df_juv$longitude[is.na(df_juv$SWM)], df_juv$latitude[is.na(df_juv$SWM)]),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

C2 = spTransform(
  C,
  CRSobj =  crs(SWM)
)

# identify nearest points with not-NA values and turn into spatial points object
C3 = nearestLand(coordinates(C2), SWM, 2000)
C3 = SpatialPoints(coords = C3, proj4string = crs(SWM))

# calculate distance between new and old points
df_juv$dist_SWM[is.na(df_juv$SWM)] = 
 spDists(
  spTransform(C2, CRSobj = "+proj=longlat +datum=WGS84") , 
  spTransform(C3, CRSobj = "+proj=longlat +datum=WGS84"), 
  longlat = TRUE,
  diagonal = TRUE)*1000

df_juv$dist_SWM[is.na(df_juv$dist_SWM)] = 0

# extract SWM values where missing
df_juv$SWM[is.na(df_juv$SWM)] = extract(SWM, C3)



#### add distances ####

# create point layer
C = SpatialPoints(
  coords = cbind(df_juv$longitude, df_juv$latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

# reproject spatial points to match distance layer
C2 = spTransform(
  C,
  CRSobj =  crs(distance_baseline)
)

# extract distance values
df_juv$distance = extract(distance_baseline, C2)









#### fix closed-in bays ####

# in some cases, the distance layer (at 5m resolution) results in bays with narrow openings being closed off (they get a value of 0 when distances are extracted)
# if so, we extract the closest value and also calculate the distance to this value, and adds it on
# Gotland is not included in the distance layer and is therefore a separate issue dealt with below

# index of problem bays
index = 
  which(
    df_juv$distance <= 0 &
      !(df_juv$longitude > 18.1 & # Gotland
          df_juv$latitude > 57.05 &  
          df_juv$latitude < 58 )
  )

# create spatial points object of missing values
C = SpatialPoints(
  coords = cbind(df_juv$longitude[index], df_juv$latitude[index]),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

C2 = spTransform(
  C,
  CRSobj =  crs(distance_baseline)
)

# identify nearest points with not-NA values and turn into spatial points object
source("help_scripts/nearestLand_0.R") # modified function which find closest non-0 value (rather than non-NA)
C3 = nearestLand_0(coordinates(C2), distance_baseline, 2000)
C3 = SpatialPoints(coords = C3, proj4string = crs(distance_baseline))

# calculate distance between new and old points
df_juv$dist_distBAS = NA
df_juv$dist_distBAS[index] = 
  spDists(
    spTransform(C2, CRSobj = "+proj=longlat +datum=WGS84") , 
    spTransform(C3, CRSobj = "+proj=longlat +datum=WGS84"), 
    longlat = TRUE,
    diagonal = TRUE)*1000

df_juv$dist_distBAS[is.na(df_juv$distBAS)] = 0

# extract distance values where missing
df_juv$distance[index] = extract(distance_baseline, C3) + df_juv$dist_distBAS[index]






#### sort out values for Gotland ####
# Gotland is not included in the distance to baseline layer


C = SpatialPoints(
  coords = cbind(df_juv$longitude, df_juv$latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

# reproject spatial points to match 5m land layer
C2 = spTransform(
  C,
  CRSobj =  crs(land5)
)

# this is the baseline in the form of points (also covers Gotland)
locs = SpatialPoints(
  coords = cbind(baseline$long, baseline$lat),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
# reproject to match 5m land layer
locs2 = spTransform(locs,
                    CRSobj =  crs(land5))


Gotland_index = which(df_juv$longitude > 18.1 & df_juv$latitude > 57.05 &  df_juv$latitude < 58)


add_on = 10000 # maximum vertical/horizontal search area


# create a progress bar
print("Adding distance to baseline values to Gotland")
pb = txtProgressBar(min = 0, 
                    max = length(Gotland_index), 
                    initial = 0,
                    style = 3) 
stepi = 1


for(i in Gotland_index){
  
  # progress bar
  setTxtProgressBar(pb,stepi)
  stepi = stepi +1
  
  # extract coordinates
  x = coordinates(C2)[i,1]
  y = coordinates(C2)[i,2]
  
  # crop out subset of land
  e = extent(x - add_on, x + add_on, y - add_on, y + add_on)
  subland = crop(land5, e)
  
  # identify cell number of closest cell to location
  longs = seq(extent(subland)[1],
              extent(subland)[2],
              length.out = dim(subland)[2])
  lats = seq(extent(subland)[3],
             extent(subland)[4],
             length.out = dim(subland)[1])
  cell.no = cellFromRowCol(subland, which.min(abs(rev(lats) - y)), which.min(abs(longs - x)))
  
  # calculate waterway distances to all water points in land layer
  subland[cell.no] = -10 # set location to 2 (so it can be used as origin in gridDistance)
  distances = gridDistance(subland, origin = -10, omit = 0)
  
  # extract distance values for points corresponding to baseline
  baseline_dists = extract(distances, locs2)
  
  # store extracted values
  if (sum(!is.na(baseline_dists)) != 0)
    df_juv$distance[i] = min(baseline_dists, na.rm = T) # take smallest value
  if (sum(!is.na(baseline_dists)) == 0)
    df_juv$distance[i] = -1 # enclosed water body
  
  
  
}

# only case with enclosed bays in Gotland (lergrav)
# use distance of closest point (~100 m away)

dists = spDistsN1(cbind(df_juv$longitude[df_juv$distance != -1], 
                        df_juv$latitude[df_juv$distance != -1]), 
                  c(df_juv$longitude[df_juv$distance == -1], 
                    df_juv$latitude[df_juv$distance == -1]), longlat = T)

df_juv$dist_distBAS[df_juv$distance == -1] = min(dists)
df_juv$distance[df_juv$distance == -1] =
  df_juv$distance[df_juv$distance != -1][which.min(dists)] + min(dists)


# rest of distances to non-NA values set to 0
df_juv$dist_distBAS[is.na(df_juv$dist_distBAS)] = 0



