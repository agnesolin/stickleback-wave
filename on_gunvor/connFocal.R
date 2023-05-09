#### add focal analysis connectivity ####

# load libraries
library(raster)
library(sp)
library(sf)
library(igraph)
library(parallel)
library(rgdal)
library(nlme)
library(doParallel)

#### sort out coordinates ####

# load data points
df_sub = read.csv("data/df_sub.csv")

# load landmask
land5 = raster("data/land5crop.tif")

# only make calculation once per location
conn_cords = unique(data.frame(longitude = df_sub$longitude, latitude = df_sub$latitude))

# save id column for later matching
ids = data.frame(point_id =
  df_sub$point_id[!duplicated(data.frame(longitude = df_sub$longitude, latitude = df_sub$latitude))])
write.csv(ids, "backup/conn_foc_id.csv", row.names = FALSE)

# convert into spatial points
C = SpatialPoints(
  coords = cbind(conn_cords$longitude, conn_cords$latitude),
  proj4string = CRS("+proj=longlat +datum=WGS84")
)
C2 = spTransform(
  C,
  CRSobj =  crs(land5)
)

# at a maximum distance of 100 m, move to water if on land
source("backup/nearestWater.R")
#newCoords5 = nearestWater(points = coordinates(C2), raster = land5, max_distance = 100)
#write.csv(newCoords5, "data/newCoords5.csv", row.names = F)
newCoords = read.csv("data/newCoords5.csv")

library(terra)
detach("package:raster", unload=TRUE)



#### set up for focal connectivity ####

# max distance for connecting habitats 
max_dist = 10000 # according to model, 95% of less than this distance

# weighting model based on migration distances
dist = c(0, 2.5, 7.5, 12.5, 17.5)

p1 = c(0, cumsum(c(74, 11, 5, 3))/100)
p2 = c(0, cumsum(c(73, 15, 9, 3))/100)
p3 = c(0, cumsum(c(83, 17, 0, 0))/100)
p4 = c(0, cumsum(c(75, 25, 0, 0))/100)
p5 = c(0, cumsum(c(78, 0, 22, 0))/100)

distWeight = data.frame(
  resp = c(p1, p2, p3, p4, p5),
  dist = rep(dist,5)
)

mod = nls(resp ~ SSasymp(dist, Asym, R0, lrc), distWeight)


# making chunks for processing (this is just to track progress and to make sure that not all is lost if process crashes)
start = seq(1, 9001, 500)
end = seq(500, 9500, 500)
end[length(end)] = nrow(newCoords)


#### using 3.2 as a limit ####


# function which calculates focal connectivity
connFUN = function(x, y){
  
  land5 = rast("data/land5crop.tif")
  habitat = rast("data/raster32.tif")
  
  
  # extent 10k in each direction
  e = ext(x - max_dist, x + max_dist, y - max_dist, y + max_dist)
  subland = crop(land5, e)
  
  # get coordinates for subset
  longs = seq(ext(subland)[1],
              ext(subland)[2],
              length.out = dim(subland)[2])
  lats = seq(ext(subland)[3],
             ext(subland)[4],
             length.out = dim(subland)[1])
  
  # find cell number of corresponding to coordinates of detonation point
  cell.no = cellFromRowCol(subland, which.min(abs(rev(lats) - y)), which.min(abs(longs - x)))
  subland[cell.no] = -10 # set location to -10 (so it can be used as origin in gridDistance)
  
  # calculate water distances detonation point
  distances = gridDistance(subland, omit = 0, origin = -10)
  
  # if closed in
  if(max(values(distances), na.rm = T) < max_dist){ 
    return(NA)
    
    # if not closed in  
  } else {
    
    # set distances further away than max distance to NA
    values(distances)[values(distances) > max_dist] = NA
    
    # translate distances into weight
    weights = 1 - predict(mod, newdata = data.frame(dist =  values(distances)/1000))
    
    # multiply weights with habitat matrix
    subhab = crop(habitat, e)
    return(sum(values(subhab)*weights, na.rm = T))
    
  }
  
}




res32_df = data.frame(x = newCoords[, 1], y = newCoords[, 2], conn32 = NA)

for(i in 1:length(start)){
  
  start_time = Sys.time()
  
  
  # extract coordinates of points
  x = newCoords[start[i]:end[i] , 1]
  y = newCoords[start[i]:end[i] , 2]
  
  cl = parallel::makeCluster(5, outfile="")
  registerDoParallel(cl)
  result = foreach(x = x, y = y,
                   .combine='list',
                   .packages = c("terra"),
                   .errorhandling = "pass") %dopar% {
                     connFUN(x,y)
                   }
  parallel::stopCluster(cl)
  
  res32_df$conn32[start[i]:end[i]] = unlist(result)
  write.csv(res32_df, "backup/res32_df.csv", row.names = FALSE)  
  
  print(paste(i,   Sys.time()-start_time ))
  
}






#### using 3.5 as a limit ####


# function which calculates focal connectivity
connFUN = function(x, y){
  
  land5 = rast("data/land5crop.tif")
  habitat = rast("data/raster35.tif")
  
  
  # extent 10k in each direction
  e = ext(x - max_dist, x + max_dist, y - max_dist, y + max_dist)
  subland = crop(land5, e)
  
  # get coordinates for subset
  longs = seq(ext(subland)[1],
              ext(subland)[2],
              length.out = dim(subland)[2])
  lats = seq(ext(subland)[3],
             ext(subland)[4],
             length.out = dim(subland)[1])
  
  # find cell number of corresponding to coordinates of detonation point
  cell.no = cellFromRowCol(subland, which.min(abs(rev(lats) - y)), which.min(abs(longs - x)))
  subland[cell.no] = -10 # set location to -10 (so it can be used as origin in gridDistance)
  
  # calculate water distances detonation point
  distances = gridDistance(subland, omit = 0, origin = -10)
  
  # if closed in
  if(max(values(distances), na.rm = T) < max_dist){ 
    return(NA)
    
    # if not closed in  
  } else {
    
    # set distances further away than max distance to NA
    values(distances)[values(distances) > max_dist] = NA
    
    # translate distances into weight
    weights = 1 - predict(mod, newdata = data.frame(dist =  values(distances)/1000))
    
    # multiply weights with habitat matrix
    subhab = crop(habitat, e)
    return(sum(values(subhab)*weights, na.rm = T))
    
  }
  
}




res35_df = data.frame(x = newCoords[, 1], y = newCoords[, 2], conn35 = NA)

for(i in 1:length(start)){
  
  start_time = Sys.time()
  

  
  # extract coordinates of points
  x = newCoords[start[i]:end[i] , 1]
  y = newCoords[start[i]:end[i] , 2]
  
  cl = parallel::makeCluster(5, outfile="")
  registerDoParallel(cl)
  result = foreach(x = x, y = y,
                   .combine='list',
                   .packages = c("terra"),
                   .errorhandling = "pass") %dopar% {
                     connFUN(x,y)
                   }
  parallel::stopCluster(cl)
  
  res35_df$conn35[start[i]:end[i]] = unlist(result)
  write.csv(res35_df, "backup/res35_df.csv", row.names = FALSE)  
  
  print(paste(i,   Sys.time()-start_time ))
  
}
