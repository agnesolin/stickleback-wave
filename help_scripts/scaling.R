#### scaling script ####

df_mod$lat_sc = scale(df_mod$latitude)
df_mod$year_sc = scale(df_mod$year)
df_mod$distance_sc = scale(df_mod$distance)
df_mod$BIAS_sc = scale(df_mod$BIASmean)
df_mod$swm_sc = scale(log10(df_mod$SWM))
df_mod$fishing_sc = scale(df_mod$fishing)
df_mod$pred_sc = scale(df_mod$totalTopPred)
df_mod$corm_sc = scale(df_mod$corm2)
df_mod$seal_sc = scale(df_mod$seal)
df_mod$temp_sc = scale(df_mod$DD)
df_mod$conn32_sc = scale(df_mod$conn32)
df_mod$conn35_sc = scale(df_mod$conn35)
df_mod$net32_sc = scale(df_mod$connected_area_net_32)
df_mod$net35_sc = scale(df_mod$connected_area_net_35)
df_mod$loc32_sc = scale(df_mod$area_net_32)
df_mod$loc35_sc = scale(df_mod$area_net_35)

