#### data filtering ####

## predators ##

## seals - cleaning the data based on distance from available seal info
table(round(df_sub$dist_seal/1000))
plot(df_sub$dist_seal, df_sub$seal)

df_sub$seal[is.na(df_sub$dist_seal)] = 0 # if further away than 10 km very likely to be zero
df_sub$seal[df_sub$dist_seal > 4000 & !is.na(df_sub$dist_seal) & df_sub$seal < 10] = 0 # if further away than 4km from data and less than 10 -> practically 0
df_sub$seal[df_sub$dist_seal > 2000 & !is.na(df_sub$dist_seal) & df_sub$seal > 4000] = NA # some odd outliers

## cormorants - cleaning the data based on distance from available cormorant info
table(round(df_sub$dist_corm/1000))
plot(df_sub$dist_corm, df_sub$cormorant)

df_sub$cormorant[is.na(df_sub$dist_corm)] = 0 # if further away than 10 km very likely to be zero
df_sub$cormorant[df_sub$dist_corm > 4000 & !is.na(df_sub$dist_corm) & df_sub$cormorant < 150] = 0 # if further away than 4km from data and less than 150 -> practically 0
df_sub$cormorant[df_sub$dist_corm > 2000 & !is.na(df_sub$dist_corm) & df_sub$cormorant > 1250] = NA # some odd outliers


## new total predator variable
df_sub$totalTopPred = df_sub$cormorant + df_sub$seal


## fishing ##
table(round(df_sub$dist_fish/1000))
plot(df_sub$dist_fish, df_sub$fishing)
# when taking data from cell more than 3 km away fishing pressure is low, so probs ok to keep as is


## SWM ##
table(round(df_sub$dist_SWM/1000))
plot(df_sub$dist_SWM, log10(df_sub$SWM))

df_sub$SWM[df_sub$dist_SWM > 100] = NA # SWM layer very high resolution, if more than a very short distance away, likely some form of inland bay/wetland-ish -> set to NA


## distance to offshore ##
table(round(df_sub$dist_distBAS/1000))
plot(df_sub$dist_distBAS, df_sub$distance)

# assuming a straight distance to closest point is fine (all short and most 0) -> keep all

