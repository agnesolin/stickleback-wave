
#### model evaluation for sdmTMB models ####

# generate residuals
if(rf) pred_fixed = mod$family$linkinv(predict(mod)$est_non_rf)
if(rf == FALSE) pred_fixed = mod$family$linkinv(predict(mod)$est)
dharm = simulate(mod, nsim = 500)
dharm = DHARMa::createDHARMa(
  simulatedResponse = dharm,
  observedResponse = observedResp,
  fittedPredictedResponse = pred_fixed
)


# look at residuals/ appropriateness of model structure
plot(dharm)
testUniformity(dharm)

hist(dharm$scaledResiduals)
testDispersion(dharm)
testDispersion(dharm, alternative = "greater") # overdispersion
testDispersion(dharm, alternative = "less") # underdispersion


# zero/ones inflation
testZeroInflation(dharm)

countZeroes <- function(x) sum(x == 0)  # testing for number of 0s
testGeneric(dharm, summary = countZeroes, alternative = "less") # 0-deficit
testGeneric(dharm, summary = countZeroes, alternative = "greater") # 0-inflation

countOnes <- function(x) sum(x == 1)  # testing for number of 1s
testGeneric(dharm, summary = countOnes, alternative = "less") # 1-deficit
testGeneric(dharm, summary = countOnes, alternative = "greater") # 1-inflation


# predictor vs residuals
testQuantiles(dharm, predictor = df_mod$distance)
testQuantiles(dharm, predictor = df_mod$BIASmean)
testQuantiles(dharm, predictor = log10(df_mod$SWM))


# testing for temporal autocorrelation
res = recalculateResiduals(dharm, group = df_mod$year)
testTemporalAutocorrelation(res, time = unique(df_mod$year))

# testing for spatial autocorrelation
groupLocations = aggregate(df_mod[, c("X_new", "Y_new")], list(factor(paste0(df_mod$X_new, df_mod$Y_new))), mean)
res = recalculateResiduals(dharm, group = factor(paste0(df_mod$X_new, df_mod$Y_new)))
DHARMa::testSpatialAutocorrelation(res, groupLocations$X_new, groupLocations$Y_new)


# multicollinearity
col_mod = glmmTMB(RPD ~ # doesn't seem to work well with sdmTMB - refit
            
            BIAS_sc +
            distance_sc +
            BIAS_sc*distance_sc  +
            swm_sc +
            (1 | yearF),
          
          data = df_mod,
          family = binomial)

mctest(col_mod)
imcdiag(col_mod, corr = TRUE)


