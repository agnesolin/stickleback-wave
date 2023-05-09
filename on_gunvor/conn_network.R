# ~ # ~ # ~ # ~ # ~ # ~ # 
##### network setup #####
# ~ # ~ # ~ # ~ # ~ # ~ #

library(igraph)
library(sf)
library(terra)
library(parallel)
library(doParallel)

polygons32 = st_read("data/polygons32.shp")
polygons35 = st_read("data/polygons35.shp")

# centroids for the polygons, also from QGIS
centroids32 = st_read("data/centroids32.shp")
centroids35 = st_read("data/centroids35.shp")





network_calc = function(centroids, nodes_df, filename){
  
  
  
  loopFUN = function(x,y){  
    
    
    # max distance for connecting habitats 
    max_dist = 10000 # according to model, 95% of less than this distance
    
    
    # load land layer
    land5 = rast("data/land5crop.tif")
    
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
    
    # find cell number of corresponding to coordinates of focal node cebtroid
    cell.no = cellFromRowCol(subland, which.min(abs(rev(lats) - y)), which.min(abs(longs - x)))
    subland[cell.no] = -10 # set location to -10 (so it can be used as origin in gridDistance)
    
    # calculate water distances from centroid of node
    distances = gridDistance(subland, omit = 0, origin = -10)
    
    # get coordinates and ids of all centroids within the cropped 10km radius area
    sub_x = coords[,1] [ coords[,1] >= e[1] & 
                           coords[,1] <= e[2] &
                           coords[,2] >= e[3] & 
                           coords[,2] <= e[4]]
    sub_y = coords[,2] [ coords[,1] >= e[1] & 
                           coords[,1] <= e[2] &
                           coords[,2] >= e[3] & 
                           coords[,2] <= e[4]]
    sub_id = nodes_df$id [ coords[,1] >= e[1] & 
                             coords[,1] <= e[2] &
                             coords[,2] >= e[3] & 
                             coords[,2] <= e[4]]
    
    
    ### sorting out link dataframe ###
    
    # set up empty data frame to store links
    res = data.frame(
      a = nodes_df$id[which(x == coords[,1] & y == coords[,2])],
      b = sub_id, 
      dist = extract(distances, cbind(sub_x, sub_y)) #  extract distances for all closeby node centroids
    ); names(res) = c("a", "b", "dist")
    
    # remove links longer than 10000 and own node
    res = res[res$dist <= 10000 & res$a != res$b,]
    
    # remove NAs
    res = res[!is.na(res$a),]
    
    # bind together
    return(res)
    
    
  }
  
  
  #### links ####
  
  coords = st_coordinates(centroids)
  
  run_seq_start = seq(1, nrow(centroids), by = 500)
  run_seq_end = c(seq(500, nrow(centroids), by = 500), nrow(centroids))

  final_res = data.frame()
  
  for(i in 1:length(run_seq_start)){
    
    print(i / length(run_seq_start) )
    
    x = coords[run_seq_start[i]:run_seq_end[i], 1]
    y = coords[run_seq_start[i]:run_seq_end[i], 2]
    
    
    cl = parallel::makeCluster(5, outfile="")
    registerDoParallel(cl)
    result = foreach(x = x, y = y,
                     .combine='rbind',
                     .packages = c("terra"),
                     .errorhandling = "pass") %dopar% {
                       loopFUN(x,y)
                     }
    parallel::stopCluster(cl)
    
    final_res = rbind(final_res, result)
    write.csv(final_res, filename, row.names = FALSE)  
    
  }
  
  
}


#### base 32 ####

centroids = centroids32
nodes_df = as.data.frame(polygons32)

nodes_df$DN = NULL; nodes_df$geometry = NULL
names(nodes_df) = c("id", "area")


network_calc(centroids, nodes_df, "backup/res32_network.csv")



#### base 35 ####

centroids = centroids35
nodes_df = as.data.frame(polygons35)

nodes_df$DN = NULL; nodes_df$geometry = NULL
names(nodes_df) = c("id", "area")


network_calc(centroids, nodes_df, "backup/res35_network.csv")


















# if not connected to any other nodes save and check

# can instead do it checking for a 

if(sum(!is.na(res$dist)) == 0) not_work = c(not_work, i)


not_work_original = not_work


# first check which ones are just loners

not_work_max_dist = rep(NA, length(not_work))

for(i in not_work){
  
  
  # re-run distance calculations
  x = coords[i,1]
  y = coords[i,2]
  
  e = extent(x - max_dist, x + max_dist, y - max_dist, y + max_dist)
  
  subland = crop(land5, e)
  
  longs = seq(extent(subland)[1],
              extent(subland)[2],
              length.out = dim(subland)[2])
  lats = seq(extent(subland)[3],
             extent(subland)[4],
             length.out = dim(subland)[1])
  
  cell.no = cellFromRowCol(subland, which.min(abs(rev(lats) - y)), which.min(abs(longs - x)))
  subland[cell.no] = -10 # set location to -10 (so it can be used as origin in gridDistance)
  
  gc()
  distances = gridDistance(subland, omit = 0, origin = -10)
  
  print(which(i == not_work)/ length(not_work))
  
  # save distance
  not_work_max_dist[which(i == not_work)] =  max(values(distances), na.rm = T)  
  
}


# remove those that include longer distances (ie issue is not being closed in but being unconnected) (if not enclosed minimum maximum distance should be diagonal to the corner)
not_work = not_work[not_work_max_dist < sqrt(5000^2+10000^2)]


# fix remaining enclosed ones
for(i in not_work[8]){
  
  x = coords[i,1]
  y = coords[i,2]
  
  e = extent(x - 12000, x + 12000, y - 12000, y + 12000) # add slightly more than max_dist since node centroid is moved
  
  subland = crop(land5, e)
  
  # identify connected water mass
  subland[median(which(values(subland) == 1))] = -10 # set one of the water cells in the middle to focal cell
  gc()
  water = gridDistance(subland, omit = 0, origin = -10)
  
  # find closest coordinates in connected water
  new_coords = nearestLand(cbind(x,y), water, max_distance = 5000)
  
  # extract new_coordinates
  x_new = new_coords[1]
  y_new = new_coords[2]
  
  # turn new and old coordinates into spatial points
  N1 = SpatialPoints(coords = cbind(x, y),
                     proj4string = crs(habitat))
  N2 = SpatialPoints(coords = cbind(x_new, y_new),
                     proj4string = crs(habitat))
  
  # calculate and print distance between old and new coordinates
  distNewOld = spDistsN1(N1,N2)
  print(paste(i, distNewOld))
  
  # repeat distance calculataions 
  e = extent(x_new - max_dist, x_new + max_dist, y_new - max_dist, y_new + max_dist)
  
  subland = crop(land5, e)
  
  longs = seq(extent(subland)[1],
              extent(subland)[2],
              length.out = dim(subland)[2])
  lats = seq(extent(subland)[3],
             extent(subland)[4],
             length.out = dim(subland)[1])
  
  cell.no = cellFromRowCol(subland, which.min(abs(rev(lats) - y_new)), which.min(abs(longs - x_new)))
  subland[cell.no] = -10 # set location to -10 (so it can be used as origin in gridDistance)
  
  gc()
  distances = gridDistance(subland, omit = 0, origin = -10)
  
  distances = distances + distNewOld
  
  
  
  sub_x = coords[,1] [ coords[,1] >= e[1] & 
                         coords[,1] <= e[2] &
                         coords[,2] >= e[3] & 
                         coords[,2] <= e[4]]
  sub_y = coords[,2] [ coords[,1] >= e[1] & 
                         coords[,1] <= e[2] &
                         coords[,2] >= e[3] & 
                         coords[,2] <= e[4]]
  sub_id = nodes_df$id [ coords[,1] >= e[1] & 
                           coords[,1] <= e[2] &
                           coords[,2] >= e[3] & 
                           coords[,2] <= e[4]]
  
  res = data.frame(
    a = nodes_df$id[i],
    b = sub_id, 
    dist = rep(NA, length(sub_id)) # calculates distance from colony in km
  )
  
  for(j in 1:length(sub_id)){
    
    cell.no = cellFromRowCol(distances, which.min(abs(rev(lats) - sub_y[j])), which.min(abs(longs - sub_x[j])))
    res$dist[j] = distances[cell.no]
  }
  
  res = res[res$dist <= 10000 & res$dist != 0,]
  
  if(sum(!is.na(res$dist)) == 0) not_work = c(not_work, i)
  res = res[!is.na(res$a),]
  
  links = rbind(links, res)
  
}

## clean up and remove duplicates ##

links_original = links # save original just in case

links$dist = round(links$dist) # round to whole meters
links = links[links$dist <= 10000,] # remove links longer than 10k

# sorting out node names so that smallest no is always first
links$a2 = apply(cbind(links$a, links$b), 1, min) 
links$b2 = apply(cbind(links$a, links$b), 1, max)
links$a = links$a2
links$b = links$b2
links$a2 = NULL; links$b2 = NULL

# remove cases of exact duplication
links$dupl = duplicated(links[, c("a", "b")]) # same node pair
links$dupl2 = duplicated(links[, c("a", "b", "dist")]) # same node pair and distance
links = links[order(links$a, links$b),]
links = links[links$dupl2 == FALSE,] 
links$dupl2 = NULL

# look at cases where same node pair is getting different values
links$diff = c(NA, diff(links$dist))
links$diff[links$dupl == FALSE] = NA
sum(!is.na(links$diff)) 
diffs = abs(links$diff[!is.na(links$diff)])
max(diffs)
mean(diffs) 
table(diffs)

# there seem to be instances where distances differ depending on direction (a->b, b->a), use shortest distance (doesn't happen often)
links = aggregate(links$dist, list(links$a, links$b), min)
names(links) = c("a", "b", "dist")

# add info on location of nodes
links = merge(links,
              data.frame(
                a = nodes_df$id,
                a_X = nodes_df$x,
                a_Y = nodes_df$y
              ),
              by = "a")

links = merge(links,
              data.frame(
                b = nodes_df$id,
                b_X = nodes_df$x,
                b_Y = nodes_df$y
              ),
              by = "b")




#### add weights ####


dist = c(0, 2.5, 7.5, 12.5, 17.5, 20)

p1 = c(0, cumsum(c(74, 11, 5, 3))/100, 1)
p2 = c(0, cumsum(c(73, 15, 9, 3))/100, 1)
p3 = c(0, cumsum(c(83, 17, 0, 0))/100, 1)
p4 = c(0, cumsum(c(75, 25, 0, 0))/100, 1)
p5 = c(0, cumsum(c(78, 0, 22, 0))/100, 1)


link_plot = data.frame(
  resp = c(p1, p2, p3, p4, p5),
  dist = rep(dist,5)
)

mod = drm(resp ~ dist, data = link_plot, fct = MM.2()) # could potentially still be improved
weight_pred = data.frame(
  dist = seq(0,20,0.01), 
  pred = predict(mod, newdata = data.frame(dist = seq(0,20,0.01))))




