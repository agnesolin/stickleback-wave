# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~ PROCESS TEMPERATURE DATA ~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# load libraries
library(data.table)
library(ncdf4)

# data are from https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=SST_BAL_SST_L4_REP_OBSERVATIONS_010_016
# co-ordinate limits used are (53,66) and (13,25)

year = 2020

# import data
temp_nc = nc_open(paste0("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/temperature/temp_input/temp_", year, ".nc"))

temp = ncvar_get(temp_nc,
                 varid = "analysed_sst")

# make into kelvin
temp = temp - 273.15

# make into easier data frame format
df_temp = data.frame(
  julian_day = rep(1:dim(temp)[3], each = dim(temp)[1]*dim(temp)[2]), 
  lat = rep(rep(temp_nc$dim$lat$vals, each = dim(temp)[1]), dim(temp)[3]),
  long = temp_nc$dim$lon$vals
)

# add temps to data frame
df_temp$temp = as.numeric(temp)
df_temp = df_temp[!is.na(df_temp$temp) , ]

# save new data frame
fwrite(df_temp, paste0("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/temperature/temp_output/temp_", year, "_processed.csv"))

