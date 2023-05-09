# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # 
# ~~~ processing of detonation data ~~~ #
# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # 



## pre-R cleaning included ##
# Removing biotestsjön (not natural conditions)
# removed some duplicates from 2019
# Changed some F/FS and associated abundances according to comments. Deleted some occasions where there were unclarities regarding number of detonations a data point was based on.
# Deleted sunken data when visibility was so poor that it was deemed to impact collection.
# Also removed few others based on notes regarding uncertainties.


#### import data ####

df = read.csv("data/YNGELDATA_20210202.csv", sep = ";", fileEncoding="UTF-8-BOM")


#### cleaning up file ####

# remove G?vlebukten missing info
df = df[df$bay != "GB",] 

# assume that inventories in Kalmar late 80s and 90s take place end of September (based on Karås reports)
df$date = as.Date(df$date)
df$date[is.na(df$date)] = as.Date(paste0(df$year[is.na(df$date)], "-09-30"))

# set stickleback to NA when FS-delvis as not clear whether S counted or not (abborre and gädda both counted floating and sunken)
df[df$FloatSink == "FS-delvis", which(names(df) == "noll_stsp_tot"):which(names(df) == "stsp_totalt")] = NA



#### calculating FS conversion factors ####

for( i in which(names(df) == "noll_abbo_tot"):which(names(df) == "noll_gadd_sjunk")){
  df[,i] = as.numeric(df[,i])} # make sure all fish abundances are numeric

df_FS = df[df$FloatSink == "FS",]

# abborre
ratio_abbo = df_FS$noll_abbo_flyt/
  df_FS$noll_abbo_tot
sum(!is.nan(ratio_abbo) & !is.infinite(ratio_abbo))
FS_abbo = round(mean(ratio_abbo[!is.nan(ratio_abbo) & !is.infinite(ratio_abbo)], na.rm = T), digits = 2) # 0.42 (0.39 in file) 

# gadda
ratio_gadd = df_FS$noll_gadd_flyt/
  df_FS$noll_gadd_tot
sum(!is.nan(ratio_gadd) & !is.infinite(ratio_gadd))
FS_gadd = round(mean(ratio_gadd[!is.nan(ratio_gadd) & !is.infinite(ratio_gadd)], na.rm = T), digits = 2)  # 0.64 (0.64 in file) 

# spigg
ratio_stsp = (df_FS$noll_stsp_flyt + df_FS$juvad_stsp_flyt)/
  (df_FS$noll_stsp_tot + df_FS$juvad_stsp_tot)
sum(!is.nan(ratio_stsp) & !is.infinite(ratio_stsp))
FS_stsp = round(mean(ratio_stsp[!is.nan(ratio_stsp) & !is.infinite(ratio_stsp)], na.rm = T), digits = 2) #0.17 (0.19 in file)



#### apply correction factors ####

# abborre
df$noll_abbo_corrected = df$noll_abbo_tot
df$noll_abbo_corrected[df$FloatSink == "F"] = df$noll_abbo_tot[df$FloatSink == "F"]/FS_abbo
df$noll_abbo_corrected = df$noll_abbo_corrected*df$Conversion.factor.Pentex25.10g

# gadda
df$noll_gadd_corrected = df$noll_gadd_tot
df$noll_gadd_corrected[df$FloatSink == "F"] = df$noll_gadd_tot[df$FloatSink == "F"]/FS_gadd
df$noll_gadd_corrected = df$noll_gadd_corrected*df$Conversion.factor.Pentex25.10g

# spigg
df$tot_stsp_corrected = df$noll_stsp_tot + df$juvad_stsp_tot
df$tot_stsp_corrected[df$FloatSink == "F"] = (df$noll_stsp_tot + df$juvad_stsp_tot)[df$FloatSink == "F"]/FS_stsp
df$tot_stsp_corrected = df$tot_stsp_corrected*df$Conversion.factor.Pentex25.10g


#### look at seasonality ####

## range to use ##
start_date = 205
end_date = 275

## spread of samples ##
range(yday(as.Date(df$date)), na.rm = T)
par(mfrow = c(1,2))
plot(table(yday(as.Date(df$date)))); abline(v = c(start_date, end_date), lty = 3)
plot(table(month(as.Date(df$date)))); abline(v = c(start_date, end_date), lty = 3)


## seasonality of presence ##
par(mfrow = c(3,2))
plot(aggregate(df$noll_abbo_corrected != 0, list(yday(as.Date(df$date))), mean),  xlab = "datum", ylab = "abborryngel", cex.lab = 1.6, pch = 16); abline(v = c(start_date, end_date), lty = 3)
plot(aggregate(df$noll_abbo_corrected, list(yday(as.Date(df$date))), mean),  xlab = "datum", ylab = "abborryngel", cex.lab = 1.6, pch = 16); abline(v = c(start_date, end_date), lty = 3)
plot(aggregate(df$noll_gadd_corrected != 0, list(yday(as.Date(df$date))), mean),  xlab = "datum", ylab = "g?ddyngel", cex.lab = 1.6, pch = 16); abline(v = c(start_date, end_date), lty = 3)
plot(aggregate(df$noll_gadd_corrected, list(yday(as.Date(df$date))), mean),  xlab = "datum", ylab = "g?ddyngel", cex.lab = 1.6, pch = 16); abline(v = c(start_date, end_date), lty = 3)
plot(aggregate(df$tot_stsp_corrected != 0, list(yday(as.Date(df$date))), mean),  xlab = "datum", ylab = "spiggyngel", cex.lab = 1.6, pch = 16); abline(v = c(start_date, end_date), lty = 3)
plot(aggregate(df$tot_stsp_corrected, list(yday(as.Date(df$date))), mean),  xlab = "datum", ylab = "spiggyngel", cex.lab = 1.6, pch = 16, ylim = c(0,250)); abline(v = c(start_date, end_date), lty = 3)
par(mfrow = c(1,1))


## variability in sampling timing ##
scatter.smooth(year(as.Date(df$date)), yday(as.Date(df$date))); abline(h = c(start_date, end_date), lty = 3)
scatter.smooth(year(as.Date(df$date[year(as.Date(df$date)) >= 2000])), 
               yday(as.Date(df$date[year(as.Date(df$date)) >= 2000]))); abline(h = c(start_date, end_date), lty = 3)
plot(aggregate(yday(as.Date(df$date)),
               list(year(as.Date(df$date))),
               mean)); abline(h = c(start_date, end_date), lty = 3)
scatter.smooth(df$DD_N, yday(as.Date(df$date))); abline(h = c(start_date, end_date), lty = 3)

sum(yday(as.Date(df$date)) >= start_date & yday(as.Date(df$date)) <= end_date)/nrow(df)

df = df[ yday(as.Date(df$date)) >= start_date & yday(as.Date(df$date)) <= end_date,]


## clean up columns and save df ##
df = df[, c("date", "bay", "DD_N", "DD_E", 
            "noll_abbo_corrected", "noll_gadd_corrected", "tot_stsp_corrected")]
names(df) = c("date", "bay", "latitude", "longitude", "perch", "pike", "stickleback")
write.csv(df, "data/yngel_final.csv", row.names = F) 



