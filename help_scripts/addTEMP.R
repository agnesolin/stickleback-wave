# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~ ADDING TEMPERATURE DATA ~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


yrs = 1982:2020 # 1982 first year of temp data, 2020 last year of detonation data


df = data.frame()

# start of sampling period
end_date = 205


for(y in yrs){
  print(y)
  
  df_y = df_juv[year(as.Date(df_juv$date)) == y, ]
  temp_data = as.data.frame(fread(paste0("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/temperature/temp_output/temp_", y ,"_processed.csv")))
  temp_data_locs = temp_data[!duplicated(temp_data[c("long", "lat")]), c("long", "lat")]
  
  locs = SpatialPoints(coords = cbind(temp_data_locs$long, temp_data_locs$lat), proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  for(i in 1:nrow(df_y)){
    C = SpatialPoints(coords = cbind(df_y$longitude[i], df_y$latitude[i]), proj4string=CRS("+proj=longlat +datum=WGS84"))
    temp_data_locs$distance = spDistsN1(locs, C, longlat = T) # calculates distance sample in km
    
    temps = temp_data$temp[temp_data$lat == temp_data_locs$lat[temp_data_locs$distance == min(temp_data_locs$distance)] &
                             temp_data$long == temp_data_locs$long[temp_data_locs$distance == min(temp_data_locs$distance)] &
                             temp_data$julian_day <= end_date]
    

    df_y$DD[i] = sum(temps[temps > 10]- 10)
    
    
    df_y$temp_comp[i] = temp_data$temp[temp_data$lat == temp_data_locs$lat[temp_data_locs$distance == min(temp_data_locs$distance)] &
                                        temp_data$long == temp_data_locs$long[temp_data_locs$distance == min(temp_data_locs$distance)] &
                                        temp_data$julian_day == yday(df_y$date[i]) ]
  }
  
  
  df = rbind(df, df_y)
  
}

# save df
write.csv(df, "data/df_juv_temp.csv", row.names = F)










