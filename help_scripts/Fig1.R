

p1 = ggplot(df_Eklof, aes(RPD)) +
  stat_bin(binwidth = 0.1, boundary = 0.05, fill = "#e7e5e1", colour = "#615c4e") + #, fill = cols
  labs(x = "Relative predator dominance", y = "Count") +
  theme_bw(base_size = 8.5) +
  theme_sets +
  theme(plot.margin = margin(0.1, 0.1, 0.6, 0.1, "cm"))


p1






## for line drawing later
# links_sub = subset(network_32, a %in% c(13215, 13190, 13253, 13288) & b %in% c(13215, 13190, 13253, 13288)) # 13215 top left, 13190 top right, 13253, midleft, 13288, bottom left
# 
# links_sub$weight = 1-predict(mod, newdata = data.frame(dist = links_sub$dist/1000))
# 
# links_sub$linewidth = NA
# links_sub$linewidth[links_sub$weight ==  max(links_sub$weight )] = 2
# links_sub$linewidth[links_sub$weight ==  min(links_sub$weight )] = 0.2
# lin_mod = lm(linewidth ~weight, data = links_sub)
# links_sub$linewidth = round(predict(lin_mod, newdata = data.frame(weight = links_sub$weight)), digits = 3)


# focal point 
x = 4814426.77
y = 4095513.22
max_dist = 10000

e = extent(x - 4000, x + 8000, y - 8000, y + 3000)
land_bg = crop(distance_baseline, e)
land_bg[land_bg>0] = 1

# get coordinates for subset
longs = seq(extent(land_bg)[1],
            extent(land_bg)[2],
            length.out = dim(land_bg)[2])
lats = seq(extent(land_bg)[3],
           extent(land_bg)[4],
           length.out = dim(land_bg)[1])

# find cell number of corresponding to coordinates of detonation point
cell.no = cellFromRowCol(land_bg, which.min(abs(rev(lats) - y)), which.min(abs(longs - x)))
land_bg[cell.no] = -10 # set location to -10 (so it can be used as origin in gridDistance)

# calculate water distances detonation point
distances = gridDistance(land_bg, omit = 0, origin = -10)
distances[distances > 10000] = NA
plot(distances)


polys = st_transform(polygons35, crs(distances))

m = distances %>% mask(polys) %>%
  crop(polys)

m_df = data.frame(rasterToPoints(m, spatial = TRUE))

m_df$conn = 1-predict(mod, newdata = data.frame(dist = m_df$layer/1000))


p2 = ggplot() +
  
  geom_raster(data = m_df , aes(x = x, y = y, fill = conn)) +
  scale_fill_gradientn(colours=rev(sequential_hcl(10, palette = "Heat 2")), name = "Connectivity strength", limits = c(0,1), guide = guide_colourbar(title.position="top", title.hjust = 0.5, ticks.colour = "black",  frame.colour = "black")) +
  

  new_scale_fill() +
  
  geom_sf(data = geom_coast,  colour = "#615c4e", fill = "#e7e5e1", size = 0) + 
  coord_sf(crs = crs(distance_baseline), xlim =  c(e[1], e[2]), ylim = c(e[3], e[4])) +

  labs(x = "", y = "") +
  
  geom_point(aes(x = x, y = y), pch = 1, size = 8) +
  
  theme_bw(base_size = 8) +
  theme_sets +
  theme(plot.margin = margin(0.5, 0.5, 0, 0.2, "cm"),
        legend.position = "bottom",
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.key.width = unit(0.5, "cm")) 
p2


## temporal trend ##
df_Eklof = df_Eklof[df_Eklof$Included == 1,]

orig_wave_model = glm(
  RPD ~
    log10(SWM) +
    DD_N +
    Distance*Year,
  data = df_Eklof,
  family = binomial)

dist_mini = crop(distance_baseline, extent(4800000, 4866000, 4074640, 4107299))
pred_df = as.data.frame(dist_mini, xy=TRUE)
coordinates(pred_df)=~x+y
names(pred_df) = "Distance"
pred_df$Distance[pred_df$Distance == 0] = NA

coordPoints = SpatialPoints(coords = coordinates(pred_df), proj4string = crs(dist_mini))
coordPoints2 = spTransform(coordPoints,   "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

pred_df$DD_N = coordinates(coordPoints2)[,2]

SWM_mini = crop(SWM5, extent(690000, 770000, 6550000, 6640000))

SWM_mini = projectRaster(SWM_mini,
                         crs = crs(dist_mini))
SWM_mini = crop(SWM_mini, extent(dist_mini))

pred_df$SWM = extract(SWM_mini, coordPoints)

empty_ras = dist_mini
values(empty_ras) = NA


pred_df_1985 = as.data.frame(pred_df)
pred_df_1985$Year = 1985
preds_1985 = predict(orig_wave_model, newdata = pred_df_1985, type = "response")
preds_1985_ras = empty_ras
values(preds_1985_ras) = preds_1985
equal_1985 = preds_1985_ras > 0.478 &  preds_1985_ras < 0.522  
equal_1985[equal_1985 != 1 & !is.na(equal_1985)] = NA

pred_df_2000 = pred_df_1985
pred_df_2000$Year = 2000
preds_2000 = predict(orig_wave_model, newdata = pred_df_2000, type = "response")
preds_2000_ras = empty_ras
values(preds_2000_ras) = preds_2000
equal_2000 = preds_2000_ras > 0.479 & preds_2000_ras < 0.521 
equal_2000[equal_2000 != 1 & !is.na(equal_2000)] = NA

pred_df_2015 = pred_df_1985
pred_df_2015$Year = 2015
preds_2015 = predict(orig_wave_model, newdata = pred_df_2015, type = "response")
preds_2015_ras = empty_ras
values(preds_2015_ras) = preds_2015
equal_2015 = preds_2015_ras > 0.48 & preds_2015_ras < 0.52
equal_2015[equal_2015 != 1 & !is.na(equal_2015)] = NA


df_85 = data.frame(rasterToPoints( equal_1985, spatial = TRUE))
df_00 = data.frame(rasterToPoints( equal_2000, spatial = TRUE))
df_15 = data.frame(rasterToPoints( equal_2015, spatial = TRUE))


p3 = ggplot() +
  geom_raster(data = df_85 , aes(x = x, y = y, fill = as.factor(layer))) +
  scale_fill_manual(values = c("#559c9e")) +
  new_scale_fill() +
  geom_raster(data = df_00 , aes(x = x, y = y, fill = as.factor(layer))) +
  scale_fill_manual(values = c("#235d72")) +
  new_scale_fill() +
  geom_raster(data = df_15 , aes(x = x, y = y, fill = as.factor(layer))) +
  scale_fill_manual(values = c("#123f5a")) +
  
  geom_sf(data = geom_coast,  colour = "#615c4e", fill = "#e7e5e1", size = 0) + 
  coord_sf(crs = crs(dist_mini), xlim =  c(4810000, 4864000), ylim = c(4076500, 4107299)) +
  
  annotate("text", label = "mainland", x = 4812000, y = 4107000, size = 3.5, colour = "#615c4e") +
  annotate("text", label = "open sea", x = 4860000, y = 4107000, size = 3.5, colour = "#615c4e") +

  annotate("text", label = "1985", x = 4855000, y = 4078000, size = 3.5, colour = "#559c9e") +
  annotate("text", label = "2000", x = 4840000, y = 4078000, size = 3.5, colour = "#235d72") +
  annotate("text", label = "2015", x = 4827000, y = 4078000, size = 3.5, colour = "#123f5a") +
  
  geom_rect(aes(xmin = e[1], xmax =  e[2], ymin = e[3], ymax = e[4]), fill = NA, color="black", alpha=0.5, linetype = "dotted") +
  
  labs(x = "", y = "") +
  
  theme_bw(base_size = 8) +
  theme_sets +
  theme(legend.position = "none") 






row2 = ggarrange(p3, p2, labels = c("b", "c"), font.label = list(face = "bold", size = 14), widths = c(1, 0.5))
ggarrange(p1, row2, ncol = 1, nrow = 2, labels = c("a", ""), font.label = list(face = "bold", size = 14))
ggsave("figures/Fig1.pdf", width = 18, height = 16, unit = "cm")


