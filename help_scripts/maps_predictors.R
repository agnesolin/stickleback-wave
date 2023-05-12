
# ~~~~~~~~~~~~~~~~~~~~~~~~ #
#### MAPS OF PREDICTORS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~ #


## sort out coordinates ##

C = SpatialPoints(
  coords = cbind(df_sub$X, df_sub$Y),
  proj4string = crs(land25)
)
C2 = spTransform(
  C,
  
  CRSobj =  crs(geom_coast)
)
df_sub$X_map = coordinates(C2)[, 1]
df_sub$Y_map = coordinates(C2)[, 2]

xlims = range(df_sub$X_map)
ylims = range(df_sub$Y_map)




## add colour scale ##

cols = colorRampPalette(brewer.pal(9, "YlGnBu"))(70)


## SWM ##

ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$SWM),], aes(x = X_map, y = Y_map, colour = log10(SWM))) + 
  
  scale_colour_gradientn(colours = cols, name = "log10(wave exposure)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_SWM.jpg", width = 18.5, height = 24, units = "cm")






## distance to offshore ##


ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$distance),], aes(x = X_map, y = Y_map, colour = distance/1000)) + 
  
  scale_colour_gradientn(colours = cols, name = "Distance from open sea (km)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_dist.jpg", width = 18.5, height = 24, units = "cm")



## BIAS ##


ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$BIASmean),], aes(x = X_map, y = Y_map, colour = BIASmean)) + 
  
  scale_colour_gradientn(colours = cols, name =  expression(Stickleback ~ open ~ sea ~ (tonnes ~ km ^ -2  ))) + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_BIAS.jpg", width = 18.5, height = 24, units = "cm")





## seals ##


ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$seal),], aes(x = X_map, y = Y_map, colour = seal)) + 
  
  scale_colour_gradientn(colours = cols, name = expression(Seal ~ predation ~ (kg ~ km ^ -2 ~ year ^ -1 ))) + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_seals.jpg", width = 18.5, height = 24, units = "cm")





## cormorants ##

ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$corm2),], aes(x = X_map, y = Y_map, colour = corm2)) + 
  
  scale_colour_gradientn(colours = cols, name = expression(Cormorant ~ predation ~ (kg ~ km ^ -2 ~ year ^ -1 ))) + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_corm.jpg", width = 18.5, height = 24, units = "cm")





## total predators ##


ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$totalTopPred),], aes(x = X_map, y = Y_map, colour = totalTopPred)) + 
  
  scale_colour_gradientn(colours = cols, name = expression(Total ~ predation ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +  
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_pred.jpg", width = 18.5, height = 24, units = "cm")







## fishing ##

ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$fishing),], aes(x = X_map, y = Y_map, colour = fishing)) + 
  
  scale_colour_gradientn(colours = cols, name = expression(Fishing ~ (kg ~ km ^ -2 ~ year ^ -1 ))) +  
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_fish.jpg", width = 18.5, height = 24, units = "cm")



## temperature ##

ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$DD),], aes(x = X_map, y = Y_map, colour = DD)) + 
  
  scale_colour_gradientn(colours = cols, name = "Temperature (degree days)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_temp.jpg", width = 18.5, height = 24, units = "cm")



## connectivity ##


## focal 3.5 ##

ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$conn35),], aes(x = X_map, y = Y_map, colour = conn35/10000)) + 
  
  scale_colour_gradientn(colours = cols, name = "Connected area (ha)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_conn35.jpg", width = 18.5, height = 24, units = "cm")


## focal 3.2 ##




ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$conn32),], aes(x = X_map, y = Y_map, colour = conn32/10000)) + 
  
  scale_colour_gradientn(colours = cols, name = "Connected area (ha)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_conn32.jpg", width = 18.5, height = 24, units = "cm")





## network 3.5 ##

ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$connected_area_net_35),], aes(x = X_map, y = Y_map, colour = connected_area_net_35/10000)) + 
  
  scale_colour_gradientn(colours = cols, name = "Connected area (ha)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )


ggsave("suppFigures/map_net35.jpg", width = 18.5, height = 24, units = "cm")





## network 3.2 ##


ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_sub[!is.na(df_sub$connected_area_net_32),], aes(x = X_map, y = Y_map, colour = connected_area_net_32/10000)) + 
  
  scale_colour_gradientn(colours = cols, name = "Connected area (ha)") + 
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )


ggsave("suppFigures/map_net32.jpg", width = 18.5, height = 24, units = "cm")





