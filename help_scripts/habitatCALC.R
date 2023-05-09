#### create underlying map for spawning habitat ####


# fix resolutions
SWM5 = resample(x = SWM, y = land5, method = "bilinear")
writeRaster(SWM5, "data/SWM5.TIF")
tmpFiles(remove = T)
rm(SWM5)

# then use QGIS to fill in empty values that seem to result from slightly different landmask being used
# raster > xx > FillNA (maximum distance = 100*5m)

SWM5 = raster("data/SWM5_NAfill.TIF")


depth5 = resample(x = depth, y = land5, method = "bilinear")
writeRaster(depth5, "//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/depth5.TIF")
tmpFiles(remove = T)
rm(depth5)
depth5 = raster("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/depth5.TIF")

#land5 = raster("data/Landmask5m.TIF") #landmask 5m res

habitat35 = overlay(land5, SWM5, depth5, fun = function(x, y, z) x == 1 & log10(y) <= 3.5 & z >= -3, progress = "text")
raster::writeRaster(habitat35, "data/habitat35.TIF", progress = "text")

habitat32 = overlay(land5, SWM5, depth5, fun = function(x, y, z) x == 1 & log10(y) <= 3.2 & z >= -3, progress = "text")
raster::writeRaster(habitat32, "data/habitat32.TIF", progress = "text")
