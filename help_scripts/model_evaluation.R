## model validation ##

# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html



# simulate residuals
simres = simulateResiduals(mod)


# look at residuals/ appropriateness of model structure
plot(simres)
testUniformity(simres)

hist(simres$scaledResiduals)
testDispersion(simres)
testDispersion(simres, alternative = "greater") # overdispersion
testDispersion(simres, alternative = "less") # underdispersion


# zero/ones inflation
testZeroInflation(simres)

countZeroes <- function(x) sum(x == 0)  # testing for number of 0s
testGeneric(simres, summary = countZeroes, alternative = "less") # 0-deficit
testGeneric(simres, summary = countZeroes, alternative = "greater") # 0-inflation

countOnes <- function(x) sum(x == 1)  # testing for number of 1s
testGeneric(simres, summary = countOnes, alternative = "less") # 1-deficit
testGeneric(simres, summary = countOnes, alternative = "greater") # 1-inflation


# predictor vs residuals
testQuantiles(simres, predictor = df_mod$year)
testQuantiles(simres, predictor = df_mod$distance)
testQuantiles(simres, predictor = df_mod$latitude)
testQuantiles(simres, predictor = df_mod$BIASmean)
testQuantiles(simres, predictor = log10(df_mod$SWM))
testQuantiles(simres, predictor = df_mod$totalTopPred)
testQuantiles(simres, predictor = df_mod$fishing)
testQuantiles(simres, predictor = df_mod$conn35)
testQuantiles(simres, predictor = df_mod$conn32)
testQuantiles(simres, predictor = df_mod$connected_area_net_35)
testQuantiles(simres, predictor = df_mod$connected_area_net_32)
testQuantiles(simres, predictor = df_mod$local_area_net_35)
testQuantiles(simres, predictor = df_mod$local_area_net_32)
testQuantiles(simres, predictor = df_mod$DD)

# testing for temporal autocorrelation
res = recalculateResiduals(simres, group = df_mod$year)
testTemporalAutocorrelation(res, time = unique(df_mod$year))


# testing for spatial autocorrelation
groupLocations = aggregate(df_mod[, c("X", "Y")], list(factor(paste0(df_mod$X, df_mod$Y))), mean)
res = recalculateResiduals(simres, group = factor(paste0(df_mod$X, df_mod$Y)))
testSpatialAutocorrelation(res, groupLocations$X, groupLocations$Y)


# multicollinearity
mctest(mod)
imcdiag(mod, corr = TRUE)

if(qqmod){
  
  
  if(resp == "RPD") jpeg("suppFigures/qq_RPD.jpg", width = 18.5, height = 18.5, units = "cm", res = 300)
  if(resp == "stick") jpeg("suppFigures/qq_stick.jpg", width = 18.5, height = 18.5, units = "cm", res = 300)
  if(resp == "pp")  jpeg("suppFigures/qq_pp.jpg", width = 18.5, height = 18.5, units = "cm", res = 300)
  
  plotQQunif(simres, 
             testUniformity = F, testOutliers = F, testDispersion = F, 
             main = "QQ-plot base model",
             family = "sans")
  
  dev.off()
  
}





