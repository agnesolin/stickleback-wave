### check distribution of response variables ###

hist(df_sub$RPD, breaks = seq(0, 1, 0.1))
ggplot(df_sub, aes(RPD)) + geom_histogram()

hist(df_sub$stickleback_int, seq(floor(min(df_sub$stickleback_int, na.rm = T)), ceiling(max(df_sub$stickleback_int, na.rm = T)), 1))
hist(df_sub$stickleback_int, seq(floor(min(df_sub$stickleback_int, na.rm = T)), ceiling(max(df_sub$stickleback_int, na.rm = T)), 1), xlim = c(0, 100)) # very many zeroes

hist(df_sub$fishPred, seq(floor(min(df_sub$fishPred, na.rm = T)), ceiling(max(df_sub$fishPred, na.rm = T)), 1))
hist(df_sub$fishPred, seq(floor(min(df_sub$fishPred, na.rm = T)), ceiling(max(df_sub$fishPred, na.rm = T)), 1), xlim = c(0, 100))



#### spatial autocorrelation ####

par(mfrow = c(1,1))

### response variables ###

## RPD ##
vario = variogram( RPD~1 , data=df_sub[!is.na(df_sub$RPD),], locations= ~X+Y)
plot(vario)
vario = variogram( RPD~1 , data=df_sub[!is.na(df_sub$RPD),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# very clear spatial structure, scale up to 3000 m or longer (120000)


## predators ##
vario = variogram( fishPred~1 , data=df_sub[!is.na(df_sub$fishPred),], locations= ~X+Y)
plot(vario)
vario = variogram( fishPred~1 , data=df_sub[!is.na(df_sub$fishPred),], locations= ~X+Y, cutoff = 3000)
plot(vario)
vario = variogram( fishPred~1 , data=df_sub[!is.na(df_sub$fishPred),], locations= ~X+Y, cutoff = 1000)
plot(vario)
# not at all the same clear structure, but possibly an increase up to ca 400 m



## stickleback ##
vario = variogram(stickleback_int~1 , data=df_sub[!is.na(df_sub$stickleback_int),], locations= ~X+Y)
plot(vario)
vario = variogram(stickleback_int~1 , data=df_sub[!is.na(df_sub$stickleback_int),], locations= ~X+Y, cutoff = 5000)
plot(vario)
vario = variogram(stickleback_int~1 , data=df_sub[!is.na(df_sub$stickleback_int),], locations= ~X+Y, cutoff = 1000)
plot(vario)
# also not a very clear structure but possibly increase up to ca 600 m


### predictors ###

## BIAS stickleback ##
vario = variogram(BIASmean~1 , data = df_sub[!is.na(df_sub$BIASmean),], locations= ~X+Y)
plot(vario)
vario = variogram(BIASmean~1 , data=df_sub[!is.na(df_sub$BIASmean),], locations= ~X+Y, cutoff = 5000)
plot(vario)
# increase over space but not very clear, ca 5000 m


## distance from offshore ##
vario = variogram(distance~1 , data = df_sub[!is.na(df_sub$distance),], locations= ~X+Y)
plot(vario)
vario = variogram(distance~1 , data=df_sub[!is.na(df_sub$distance),], locations= ~X+Y, cutoff = 100000)
plot(vario)
# very clear up to ca 100000 m


## SWM ##
vario = variogram(SWM~1 , data = df_sub[!is.na(df_sub$SWM),], locations= ~X+Y)
plot(vario)
vario = variogram(SWM~1 , data=df_sub[!is.na(df_sub$SWM),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# not very clear but still a bit of an increase up to 3000 m


## predation ##
vario = variogram(totalTopPred~1 , data = df_sub[!is.na(df_sub$totalTopPred),], locations= ~X+Y)
plot(vario)
vario = variogram(totalTopPred~1 , data=df_sub[!is.na(df_sub$totalTopPred),], locations= ~X+Y, cutoff = 20000)
plot(vario)
# clear increase over up to ca 100000 m (makes sense with how the variable is created)


## fishing ##
vario = variogram(fishing~1 , data = df_sub[!is.na(df_sub$fishing),], locations= ~X+Y)
plot(vario)
vario = variogram(fishing~1 , data=df_sub[!is.na(df_sub$fishing),], locations= ~X+Y, cutoff = 6000)
plot(vario)
# no clear pattern over larger scales, bu clear at short scales Up to ca 6000 m


## conn 35 focal ##
vario = variogram(conn35~1 , data = df_sub[!is.na(df_sub$conn35),], locations= ~X+Y)
plot(vario)
vario = variogram(conn35~1 , data=df_sub[!is.na(df_sub$conn35),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# clear increase up  to ca 100000 m

## conn 32 focal ##
vario = variogram(conn32~1 , data = df_sub[!is.na(df_sub$conn32),], locations= ~X+Y)
plot(vario)
vario = variogram(conn32~1 , data=df_sub[!is.na(df_sub$conn32),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# clear increase up  to ca 150000 m

## conn 35 net ##
vario = variogram(connected_area_net_35~1 , data = df_sub[!is.na(df_sub$connected_area_net_35),], locations= ~X+Y)
plot(vario)
vario = variogram(connected_area_net_35~1 , data=df_sub[!is.na(df_sub$connected_area_net_35),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# clear increase up  to ca 150000 m


## conn 32 net ##
vario = variogram(connected_area_net_32~1 , data = df_sub[!is.na(df_sub$connected_area_net_32),], locations= ~X+Y)
plot(vario)
vario = variogram(connected_area_net_32~1 , data=df_sub[!is.na(df_sub$connected_area_net_32),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# clear increase up  to ca 150000 m


## conn 35 net local ##
vario = variogram(area_net_35~1 , data = df_sub[!is.na(df_sub$area_net_35),], locations= ~X+Y)
plot(vario)
vario = variogram(area_net_35~1 , data=df_sub[!is.na(df_sub$area_net_35),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# no clear pattern over larger scales (potentially increase), but increases up to 3000 m

## conn 32 net local ##
vario = variogram(area_net_32~1 , data = df_sub[!is.na(df_sub$area_net_32),], locations= ~X+Y)
plot(vario)
vario = variogram(area_net_32~1 , data=df_sub[!is.na(df_sub$area_net_32),], locations= ~X+Y, cutoff = 3000)
plot(vario)
# clear increase up  to ca 150000 m



#### look at distances between samples ####

par(mfrow = c(4,5))


## how many samples per latitudinal bin ##

for(y in sort(unique(df_sub$year))){
  
  sub = df_sub[df_sub$year == y,]
  hist(sub$X, xlab = "Latitude", main = y, breaks = seq(459000, 749000, 1000), ylim = c(0,100))
  
}
# pretty good spread but lack of low-latitude samples in early years


par(mfrow = c(1,1))

for(y in sort(unique(df_sub$year))){
  
  
  sub = df_sub[df_sub$year == y,]
  
  dists = spDists(as.matrix(cbind(sub$X, sub$Y)), longlat = FALSE, diagonal = FALSE)
  sub$point_id[colSums(dists == 0) > 1]
  
  diag(dists) = NA
  
  
  print(sum(dists < 1000 & !is.na(dists))/( sum(!is.na(dists))/1 ))
  
  plot(table(round(dists/1000)), main = y, xlab = "km", ylab = "Freq")
  
  plot(table(round(dists[dists < 1000])), main = y, xlab = "m", ylab = "Freq")
  
 
}



#### collinearity between explanatory variables ####

pairs(cbind(dist = df_sub$distance, 
            
            BIAS = df_sub$BIASmean,
            
            SWM = log10(df_sub$SWM),
            
            top_pred = df_sub$totalTopPred,
            
            fishing = df_sub$fishing,
            
            connectivity = df_sub$conn35 
))

pairs(cbind(dist = df_sub$distance, 
            
            BIAS = df_sub$BIASmean,
            
            SWM = log10(df_sub$SWM),
            
            top_pred = df_sub$totalTopPred,
            
            fishing = df_sub$fishing,
            
            connectivity = df_sub$conn32 
))


pairs(cbind(dist = df_sub$distance, 
            
            BIAS = df_sub$BIASmean,
            
            SWM = log10(df_sub$SWM),
            
            top_pred = df_sub$totalTopPred,
            
            fishing = df_sub$fishing,
            
            connectivity = df_sub$connected_area_net_35
))


pairs(cbind(dist = df_sub$distance, 
            
            BIAS = df_sub$BIASmean,
            
            SWM = log10(df_sub$SWM),
            
            top_pred = df_sub$totalTopPred,
            
            fishing = df_sub$fishing,
            
            connectivity = df_sub$connected_area_net_32
))


pairs(cbind(dist = df_sub$distance, 
            
            BIAS = df_sub$BIASmean,
            
            SWM = log10(df_sub$SWM),
            
            top_pred = df_sub$totalTopPred,
            
            fishing = df_sub$fishing,
            
            loc_area = df_sub$area_net_35
))


pairs(cbind(dist = df_sub$distance, 
            
            BIAS = df_sub$BIASmean,
            
            SWM = log10(df_sub$SWM),
            
            top_pred = df_sub$totalTopPred,
            
            fishing = df_sub$fishing,
            
            loc_area = df_sub$area_net_32
))



# cor tables # 
df_cor = data.frame(conn35 = df_sub$conn35,
                              conn32 = df_sub$conn32,
                              connected_area_net_35 = df_sub$connected_area_net_35,
                              connected_area_net_32 = df_sub$connected_area_net_32,
                              year = df_sub$year, 
                              BIAS = df_sub$BIASmean, 
                              distance = df_sub$distance, 
                              latitude = df_sub$latitude, 
                              SWM = df_sub$SWM, 
                              predation = df_sub$totalTopPred, 
                              fishing = df_sub$fishing, 
                              temperature = df_sub$DD)


df_cor = df_cor[rowSums(is.na(df_cor)) == 0 ,]
nrow(df_cor)

# conn 35
cor.var = data.frame(connectivity = df_cor$conn35,
                     year = df_cor$year, 
                     BIAS = df_cor$BIAS, 
                     distance = df_cor$distance, 
                     latitude = df_cor$latitude, 
                     SWM = df_cor$SWM, 
                     predation = df_cor$predation, 
                     fishing = df_cor$fishing, 
                     temperature = df_cor$temperature)

cors = as.data.frame(cor(cor.var))
rcorr(as.matrix(cor.var),type="pearson")


tab_df(cors,show.rownames = TRUE,
       file="suppTables/cors_conn35.doc")


# conn 32
cor.var = data.frame(connectivity = df_cor$conn32,
                     year = df_cor$year, 
                     BIAS = df_cor$BIAS, 
                     distance = df_cor$distance, 
                     latitude = df_cor$latitude, 
                     SWM = df_cor$SWM, 
                     predation = df_cor$predation, 
                     fishing = df_cor$fishing, 
                     temperature = df_cor$temperature)


cors = as.data.frame(cor(cor.var))
rcorr(as.matrix(cor.var),type="pearson")


tab_df(cors,show.rownames = TRUE,
       file="suppTables/cors_conn32.doc")



# net 35
cor.var = data.frame(connectivity = df_cor$connected_area_net_35,
                     year = df_cor$year, 
                     BIAS = df_cor$BIAS, 
                     distance = df_cor$distance, 
                     latitude = df_cor$latitude, 
                     SWM = df_cor$SWM, 
                     predation = df_cor$predation, 
                     fishing = df_cor$fishing, 
                     temperature = df_cor$temperature)


cors = as.data.frame(cor(cor.var))
rcorr(as.matrix(cor.var),type="pearson")


tab_df(cors,show.rownames = TRUE,
       file="suppTables/cors_net35.doc")




# net 32
cor.var = data.frame(connectivity = df_cor$connected_area_net_32,
                     year = df_cor$year, 
                     BIAS = df_cor$BIAS, 
                     distance = df_cor$distance, 
                     latitude = df_cor$latitude, 
                     SWM = df_cor$SWM, 
                     predation = df_cor$predation, 
                     fishing = df_cor$fishing, 
                     temperature = df_cor$temperature)



cors = as.data.frame(cor(cor.var))
rcorr(as.matrix(cor.var),type="pearson")


tab_df(cors,show.rownames = TRUE,
       file="suppTables/cors_net32.doc")







