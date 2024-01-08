

# function to calculate connected area
calc_connected_area = function(nodes, centroids, links){
  
  # extracts node ids and areas
  nodes_df = data.frame(
    id = as.data.frame(nodes)$fid,
    area = as.data.frame(nodes)$area
  )
  
  # adds coordinates from centroids
  coords = st_coordinates(centroids)[, 1:2]
  nodes_df = cbind(nodes_df, coords)
  
  # add empty variable for connected area
  nodes_df$connected_area = NA
  
  test = data.frame()
  
  # loop through all nodes  
  for(i in 1:nrow(nodes_df)){
    
    # subset link info
    id = nodes_df$id[i]
    link_sub = links[links$a == id, ]
    
    # merge with area info
    link_sub = 
      merge(link_sub, 
            nodes_df[, c("id", "area")],
            by.x = "b", by.y = "id",
            all.x = T,
            all.y = F)
    
    # calculate distance-based weight based on model fit earlier
    link_sub$weight = 
      1-predict(mod, newdata = data.frame(dist = link_sub$dist/1000))
    
    test = rbind(test, link_sub)
    
    # calculate connected area
    nodes_df$connected_area[i] = sum(link_sub$weight*link_sub$area)
    
  }
  return(nodes_df) 
}





#### calculate values for different network representations ####
nodes = polygons32
centroids = centroids32
links = network_32
res_network_32 = calc_connected_area(nodes, centroids, links)


nodes = polygons35
centroids = centroids35
links = network_35
res_network_35 = calc_connected_area(nodes, centroids, links)



#### make some histograms ####


# make histogram of area

maxV = max(c(res_network_35$area, res_network_32$area), na.rm = T)/10000

p1 = ggplot(data = res_network_35, aes(area/10000)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "Patch size (ha)", y = "Count", title = "Limit = 3.5") +
  lims(x = c(0, maxV), y = c(0, 610)) +
  theme_sets 

p2 = ggplot(data = res_network_32, aes(area/10000)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "Patch size (ha)", y = "Count", title = "Limit = 3.2") +
  lims(x = c(0, maxV), y = c(0, 610)) +
  theme_sets 


ggarrange(p1, p2,  labels = c("a.", "b."), font.label = list(size = 14, family = "sans", face = "plain"))
ggsave("suppFigures/patchAreaHist.jpg", width = 18.5, height = 9, units = "cm")


mean(res_network_35$area/10000)
quantile(res_network_35$area/10000, 0.95)
range(res_network_35$area/10000)

mean(res_network_32$area/10000)
quantile(res_network_32$area/10000, 0.95)
range(res_network_32$area/10000)


(10000*10000*pi)/10000


# make histogram of link length

p1 = ggplot(data = network_35, aes(dist/1000)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  labs(x = "Link distance (km)", y = "Count", title = "Limit = 3.5") +
  lims(x = c(0, 10), y = c(0, 2500)) +
  theme_sets 

p2 = ggplot(data = network_32, aes(dist/1000)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  labs(x = "Link distance (km)", y = "Count", title = "Limit = 3.2") +
  lims(x = c(0, 10), y = c(0, 2500)) +
  theme_sets 



p3 = ggplot(data = network_35, aes(1-predict(mod, newdata = data.frame(dist = dist/1000)))) +
  geom_histogram(binwidth = 0.001) +
  theme_bw() +
  labs(x = "Link weight", y = "Count", title = "Limit = 3.5") +
  lims(x = c(0, 1), y = c(0, 3000)) +
  theme_sets 


p4 = ggplot(data = network_32, aes(1-predict(mod, newdata = data.frame(dist = dist/1000)))) +
  geom_histogram(binwidth = 0.001) +
  theme_bw() +
  labs(x = "Link weight", y = "Count", title = "Limit = 3.2") +
  lims(x = c(0, 1), y = c(0, 3000)) +
  theme_sets 

ggarrange(p1, p2, p3, p4, labels = c("a.", "b.", "c.", "d."), font.label = list(size = 14, family = "sans", face = "plain"))
ggsave("suppFigures/linkLengthHist.jpg", width = 18.5, height = 18.5, units = "cm")








#### extract values ####



## 3.2 basic ##

# find closest nodes
closest = nn2(
  res_network_32[, c("X", "Y")],
  data.frame(X = df_sub$X, Y = df_sub$Y),
  k = 1
)

# extract metrics
index_closest = closest$nn.idx
df_sub$area_net_32 = res_network_32$area[index_closest]
df_sub$connected_area_net_32 = res_network_32$connected_area[index_closest]
df_sub$node_id_32 = res_network_32$id[index_closest]



## 3.5 basic ##

# find closest nodes
closest = nn2(
  res_network_35[, c("X", "Y")],
  data.frame(X = df_sub$X, Y = df_sub$Y),
  k = 1
)

# extract metrics
index_closest = closest$nn.idx
df_sub$area_net_35 = res_network_35$area[index_closest]
df_sub$connected_area_net_35 = res_network_35$connected_area[index_closest]
df_sub$node_id_35 = res_network_35$id[index_closest]




## histogram of extracted values


maxV = max(c(df_sub$connected_area_net_32[!is.infinite(df_sub$connected_area_net_32)], df_sub$connected_area_net_35[!is.infinite(df_sub$connected_area_net_35)]), na.rm = T)/10000

p1 = ggplot(data = df_sub, aes(connected_area_net_35/10000)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "Connected area (ha)", y = "Count", title = "Limit = 3.5") +
  lims(x = c(0, maxV), y = c(0, 900)) +
  theme_sets

p2 = ggplot(data = df_sub, aes(connected_area_net_32/10000)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "Connected area (ha)", y = "Count", title = "Limit = 3.2") +
  lims(x = c(0, maxV), y = c(0, 900)) +
  theme_sets 


ggarrange(p1, p2,  labels = c("a.", "b."), font.label = list(size = 14, family = "sans", face = "plain"))
ggsave("suppFigures/netConnHist.jpg", width = 18.5, height = 9, units = "cm")




#### comparison habitat availability and connectivity ####


p1 = ggplot(data = df_sub, aes(x = conn35/10000, y = connected_area_net_35/10000)) +
  geom_point(size = 0.5) +
  theme_bw() +
  labs(x = "Connected area (ha) - focal", y = "Connected area (ha) - network", title = "Limit = 3.5") +
  theme_sets 

p2 = ggplot(data = df_sub, aes(x = conn32/10000, y = connected_area_net_32/10000)) +
  geom_point(size = 0.5) +
  theme_bw() +
  labs(x = "Connected area (ha) - focal", y = "Connected area (ha) - network", title = "Limit = 3.2") +
  theme_sets 

p3 = ggplot(data = df_sub, aes(x = conn35/10000, y = conn32/10000)) +
  geom_point(size = 0.5) +
  theme_bw() +
  labs(x = "Limit = 3.5", y = "Limit = 3.2", title = "Connected area (ha) - focal") +
  theme_sets 


p4 = ggplot(data = df_sub, aes(x = connected_area_net_35/10000, y = connected_area_net_32/10000)) +
  geom_point(size = 0.5) +
  theme_bw() +
  labs(x = "Limit = 3.5", y = "Limit = 3.2", title = "Connected area (ha) - network") +
  theme_sets 

ggarrange(p1, p2, p3, p4, labels = c("a.", "b.", "c.", "d."), font.label = list(size = 14, family = "sans", face = "plain"))
ggsave("suppFigures/connComparison.jpg", width = 18.5, height = 18.5, units = "cm")













