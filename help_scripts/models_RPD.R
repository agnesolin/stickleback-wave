
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### RPD: wave, spatio-temporal patterns, drivers of residuals ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#### create model data frame ####


# subset so that no values are NA
df_mod = df_sub[!is.na(df_sub$totalTopPred) & !is.na(df_sub$SWM) & !is.na(df_sub$RPD) & !is.na(df_sub$BIASmean) & !is.na(df_sub$conn35),]
colSums(is.na(df_mod))
nrow(df_mod)

# save as figure source data
write.csv(df_mod, 
          "source-data-figures/source_raw-data-RPD_fig2-3.csv",
          row.names = FALSE)

# scale values
source("help_scripts/scaling.R")

#### create mesh for fitting spatial models ####


# load polygon defining extent of water
water_poly = readOGR("//storage.slu.se/home_2$/agin0002/Documents/SLU-Agnes/general_GIS/water_polygons/water_poly.shp")

# check if data points are in water
wat = extract(land25, cbind(df_mod$X, df_mod$Y))

needs_fixing = which(!is.na(wat))

# move points on land to closest point in water
source("help_scripts/nearestLand_NA.R")
in_water = nearestLand_NA(points = cbind(df_mod$X[needs_fixing], df_mod$Y[needs_fixing]), raster = land25, max_distance = 150)

df_mod$X_new = df_mod$X
df_mod$X_new[needs_fixing] = in_water[,1]
df_mod$Y_new = df_mod$Y
df_mod$Y_new[needs_fixing] = in_water[,2]


# parameter values for mesh based on https://rpubs.com/jafet089/886687 and https://haakonbakkagit.github.io/btopic104.html#what-is-needed-for-a-good-mesh
bound.outer = diff(range(df_mod$X))/3
max.edge = bound.outer/5

# define mesh
inla_mesh = inla.mesh.2d(
  
  # defines boundaries and make sure no spatial dependencies across land
  boundary = water_poly,
  
  # location of points
  loc = cbind(df_mod$X_new, df_mod$Y_new),
  
  # max edge = longest allowed triangle length, first value is inner domain, second is outer extension
  max.edge = c(1,2)*max.edge, 
  
  # this controls the shortest allowed distance between points (clusters points close together)
  cutoff = 100, # reasonably small and comparable to resolution of boundary layer
  
  # controls how much domain should be extended to avoid boundary effects
  offset = c(max.edge, bound.outer)
)

mesh = make_mesh(df_mod, xy_cols = c("X_new", "Y_new"), mesh = inla_mesh)

mesh = add_barrier_mesh(
  mesh,
  st_as_sf(water_poly),
  plot = TRUE
)


#### fit models ####

# create factor year
df_mod$yearF = as.factor(df_mod$year)

# set up cluster structure
clust = sample(1:5, size = nrow(df_mod), replace = T)

#### BASE MODEL ####

# fit model
base_model = sdmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "off",
  
  family = binomial(),
  weights = rep(1, nrow(df_mod)))

# estimates and AIC
tidy(base_model, conf.int = T)
AIC(base_model) # 4008.1

# model evaluation
# mod = base_model
# observedResp = df_mod$RPD
# rf = FALSE
# source("help_scripts/model_evaluation_sdmTMB.R")

#cross-validation
m_cv = sdmTMB_cv(
  RPD ~
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "off",
  
  fold_ids = clust,
  k_folds = length(unique(clust)),
  
  family = binomial())

m_cv$elpd # -0.6235506



#### BASIC RANDOM YEAR EFFECT ####

# fit model
temp_model = sdmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    (1 | yearF),
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "off",
  
  family = binomial(),
  weights = rep(1, nrow(df_mod)))

# estimates and AIC
tidy(temp_model, conf.int = T)
AIC(temp_model) # 3785.353

# model evaluation
# mod = temp_model
# observedResp = df_mod$RPD
# rf = FALSE
# source("help_scripts/model_evaluation_sdmTMB.R")


# cross-validation
m2_cv = sdmTMB_cv(
  RPD ~
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc+
    (1 | yearF),
  
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "off",
  
  fold_ids = clust,
  k_folds = length(unique(clust)),
  
  family = binomial())

m2_cv$elpd # -0.5688644




#### CONSTANTS SPATIAL FIELDS ####


# fit (or load) model
# spat_model = sdmTMB(
#   RPD ~
# 
#     BIAS_sc +
#     distance_sc +
#     BIAS_sc*distance_sc  +
#     swm_sc,
# 
#   data = df_mod,
# 
#   mesh = mesh,
#   spatial = "on",
# 
#   family = binomial(),
#   weights = rep(1, nrow(df_mod)))
# 
# save(spat_model , file = "models/spat_model.Rdata")

load("models/spat_model.Rdata")


# estimates and AIC
tidy(spat_model, conf.int = T) 
AIC(spat_model) # 2812.367

# model evaluation
# mod = spat_model
# observedResp = df_mod$RPD
# rf = TRUE
# source("help_scripts/model_evaluation_sdmTMB.R")



# cross-validation
# m3_cv = sdmTMB_cv(
#   RPD ~
#     
#     BIAS_sc +
#     distance_sc +
#     BIAS_sc*distance_sc  +
#     swm_sc,
#   
#   data = df_mod,
#   
#   mesh = mesh,
#   spatial = "on",
#   
#   fold_ids = clust,
#   k_folds = length(unique(clust)),
#   
#   family = binomial())
# 
# 
# m3_cv$elpd # -0.3726271




#### SPATIAL FIELDS WITH CROSS-AREA RANDOM YEAR EFFECT ####


# fit (or load) model
# spat_temp_model_1 = sdmTMB(
#   RPD ~
# 
#     BIAS_sc +
#     distance_sc +
#     BIAS_sc*distance_sc  +
#     swm_sc +
#     (1 | yearF),
# 
#   data = df_mod,
# 
#   mesh = mesh,
#   spatial = "on",
# 
# 
#   family = binomial(),
#   weights = rep(1, nrow(df_mod)))
# 
# 
# save(spat_temp_model_1, file = "models/spat_temp_model_1.Rdata")
load("models/spat_temp_model_1.Rdata")


# estimates and AIC
tidy(spat_temp_model_1, conf.int = T)
AIC(spat_temp_model_1) # 2598.503




# model evaluation
# mod = spat_temp_model_1
# observedResp = df_mod$RPD
# rf = TRUE
# source("help_scripts/model_evaluation_sdmTMB.R")


#cross-validation
# plan(multisession, workers = 5)
# 
# m4_cv = sdmTMB_cv(
#   RPD ~
#     
#     BIAS_sc +
#     distance_sc +
#     BIAS_sc*distance_sc  +
#     swm_sc +
#     (1 | yearF),
#   
#   data = df_mod,
#   
#   mesh = mesh,
#   spatial = "on",
#   
#   
#   fold_ids = clust,
#   k_folds = length(unique(clust)),
#   
#   
#   family = binomial())
# 
# plan(sequential)
# 
# m4_cv$elpd # -0.3376445



#### SPATIO-TEMPORAL FIELDS ####


# fit (or load) model
# spat_temp_model_2 = sdmTMB(
#   RPD ~
#     
#     BIAS_sc +
#     distance_sc +
#     BIAS_sc*distance_sc  +
#     swm_sc,
#   
#   data = df_mod,
#   
#   time = "year",
#   spatiotemporal = "iid",
#   
#   mesh = mesh,
#   spatial = "off",
#   
#   
#   family = binomial(),
#   weights = rep(1, nrow(df_mod)))
# 
# 
# save(spat_temp_model_2 , file = "models/spat_temp_model_2.Rdata")

load("models/spat_temp_model_2.Rdata")

# estimates and AIC
tidy(spat_temp_model_2, conf.int = T)
AIC(spat_temp_model_2) # 2581.013

# model evaluation
# mod = spat_temp_model_2
# observedResp = df_mod$RPD
# rf = TRUE
# source("help_scripts/model_evaluation_sdmTMB.R")

# cross-validation
# plan(multisession, workers = 5)
# 
# m5_cv = sdmTMB_cv(
#   RPD ~
#     
#     BIAS_sc +
#     distance_sc +
#     BIAS_sc*distance_sc  +
#     swm_sc,
#   
#   data = df_mod,
#   
#   time = "year",
#   spatiotemporal = "iid",
#   
#   mesh = mesh,
#   spatial = "off",
#   
#   parallel = TRUE,
#   
#   fold_ids = clust,
#   k_folds = length(unique(clust)),
#   
#   family = binomial())
# 
# plan(sequential)
# 
# m5_cv$elpd




#### saving AICs and elpds ####

# RPD_spat_temp_model_comp = data.frame(
#   RPD_AIC = c(AIC(base_model) - AIC(base_model, temp_model, spat_model, spat_temp_model_1, spat_temp_model_2)),
#   RPD_elpds = c(m_cv$elpd, m2_cv$elpd, m3_cv$elpd, m4_cv$elpd, NA) 
#   
# )




# #### visualisation  ####

cols = c("darkslategrey", "darkolivegreen4", "darkgoldenrod")


range(df_mod$distance)
dist_breaks = seq(0, 35000, 10000)
dist_breaks_scaled = (dist_breaks-attributes(df_mod$distance_sc)$`scaled:center`)/
  attributes(df_mod$distance_sc)$`scaled:scale`
dist_breaks = dist_breaks/1000


quant_breaks = quantile(df_mod$BIAS_sc, c(0, 0.33, 0.66, 1))
quant_breaks = quant_breaks + c(-1, 0, 0, 1)
df_mod$BIAS_cut = as.factor(cut(df_mod$BIAS_sc, breaks = quant_breaks))



visreg(spat_temp_model_2,
       xvar = "distance_sc",
       by = "BIAS_sc",
       scale = "response", overlay = T, 
       gg = T, nn = 1000, line.par = list(size = 1, lty = rep(c(1,2,6), each = 1000)), partial = F, rug = F) +
  
  scale_colour_manual(values = cols, labels = c("low", "mid", "high"), name = "Offshore stickleback",   guide = guide_legend(
    override.aes = list(
      linetype = c("solid", "dashed", "twodash")
    ))) +
  scale_fill_manual(values = alpha(cols, 0.5), labels = c("low", "mid", "high"), name = "Offshore stickleback") +
  
  new_scale_colour() +
  
  geom_point(data = df_mod, aes(x = distance_sc, y = RPD, colour = BIAS_cut), alpha = 0.4, size = 0.8) +
  scale_colour_manual(values = cols, labels = c("low", "mid", "high"), name = "Offshore stickleback") +
  
  
  scale_x_continuous(breaks = dist_breaks_scaled, labels = dist_breaks) +
  
  labs(x = "Distance from the open sea (km)", y = "Probability of predator dominance") +
  
  
  coord_cartesian(clip = 'off') +
  ylim(c(0,1)) +
  
  
  theme_bw(base_size = 9)+
  theme_sets +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.7, 'cm') 
  ) 



ggsave("suppFigures/stickleback_wave.png", width = 9, height = 10, units = "cm")



#### ADD IN POTENTIAL DRIVERS ####



mod_ref = glmmTMB( # ref model with only stickleback wave-variables
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    (1 | year),
  
  data = df_mod,
  
  family = binomial)


# mod = mod_ref
# qqmod = FALSE
# source("help_scripts/model_evaluation.R")


mod_ref2 = glmmTMB( # second reference model with bo stickleback wave-drivers
  RPD ~ 
    
    
    conn35_sc +
    pred_sc +
    fishing_sc +
    conn35_sc*pred_sc +
    conn35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    (1 | year),
  
  data = df_mod,
  
  family = binomial)


# mod = mod_ref2
# qqmod = FALSE
# source("help_scripts/model_evaluation.R")


#### conn 35  ####


mod_RPD_conn35_full = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    conn35_sc +
    pred_sc +
    fishing_sc +
    conn35_sc*pred_sc +
    conn35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)


mod = mod_RPD_conn35_full
qqmod = TRUE
resp = "RPD"
source("help_scripts/model_evaluation.R")

# dredge_RPD_conn35_full = dredge(mod_RPD_conn35_full, trace = 2)
# tab_df(dredge_RPD_conn35_full,
#        file="suppTables/dredge_conn35_full.doc")


mod_RPD_conn35_noFish = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    conn35_sc +
    pred_sc +
    conn35_sc*pred_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  data = df_mod,
  
  family = binomial)

tidy(mod_RPD_conn35_noFish, conf.int = T)
tidy(mod_RPD_conn35_full, conf.int = T)[c(1:6,8:10,12:13),]




mod_RPD_conn35_noConn = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    pred_sc +
    fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)


#### conn 32  ####


mod_RPD_conn32_full = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    
    conn32_sc +
    pred_sc +
    fishing_sc +
    conn32_sc*pred_sc +
    conn32_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  data = df_mod,
  
  na.action = "na.fail",
  
  family = binomial)

# mod = mod_RPD_conn32_full
# qqmod = FALSE
# source("help_scripts/model_evaluation.R")

# dredge_RPD_conn32_full = dredge(mod_RPD_conn32_full, trace = 2)
# tab_df(dredge_RPD_conn32_full,
#        file="suppTables/dredge_conn32_full.doc")


mod_RPD_conn32_noFish = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    
    conn32_sc +
    pred_sc +
    conn32_sc*pred_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  data = df_mod,
  
  family = binomial)

tidy(mod_RPD_conn32_noFish, conf.int = T)
tidy(mod_RPD_conn32_full, conf.int = T)[c(1:6,8:10,12:13),]


#### net 35  ####


mod_RPD_net35_full = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    net35_sc +
    pred_sc +
    fishing_sc +
    net35_sc*pred_sc +
    net35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    
    (1 | year),
  
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)

# mod = mod_RPD_net35_full
# qqmod = FALSE
# source("help_scripts/model_evaluation.R")

# dredge_RPD_net35_full = dredge(mod_RPD_net35_full, trace = 2)
# tab_df(dredge_RPD_net35_full,
#        file="suppTables/dredge_net35_full.doc")


mod_RPD_net35_noFish = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    net35_sc +
    pred_sc +
    net35_sc*pred_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    
    (1 | year),
  
  data = df_mod,
  
  family = binomial)

tidy(mod_RPD_net35_noFish, conf.int = T)
tidy(mod_RPD_net35_full, conf.int = T)[c(1:6,8:10,12:13),]




#### net 32  ####


mod_RPD_net32_full = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    net32_sc +
    pred_sc +
    fishing_sc +
    net32_sc*pred_sc +
    net32_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    
    (1 | year),
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)

# mod = mod_RPD_net32_full
# qqmod = FALSE
#source("help_scripts/model_evaluation.R")


# dredge_RPD_net32_full = dredge(mod_RPD_net32_full, trace = 2)
# tab_df(dredge_RPD_net32_full,
#        file="suppTables/dredge_net32_full.doc")


mod_RPD_net32_noFish = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    net32_sc +
    pred_sc +
    net32_sc*pred_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    
    (1 | year),
  
  data = df_mod,
  
  family = binomial)

tidy(mod_RPD_net32_noFish, conf.int = T)
tidy(mod_RPD_net32_full, conf.int = T)[c(1:6,8:10,12:13),]





# SAVE ESTIMATES AND COEFFICIENTS #

# AIC and R2 table #

mod_sel_tab = AIC(mod_ref, mod_ref2, mod_RPD_conn35_full, mod_RPD_conn32_full, mod_RPD_net35_full, mod_RPD_net32_full)
mod_sel_tab$R2 = c(
  r.squaredGLMM(mod_ref)[2,1],
  r.squaredGLMM(mod_ref2)[2,1],
  r.squaredGLMM(mod_RPD_conn35_full)[2,1],
  r.squaredGLMM(mod_RPD_conn32_full)[2,1],
  r.squaredGLMM(mod_RPD_net35_full)[2,1],
  r.squaredGLMM(mod_RPD_net32_full)[2,1])
mod_sel_tab$dAIC = mod_sel_tab$AIC-mod_sel_tab$AIC[1]
mod_sel_tab = mod_sel_tab[, c(3:4)]

tab_df(mod_sel_tab,
       file="suppTables/mod_sel_RPD_local_drivers.doc")

r.squaredGLMM(mod_RPD_conn35_noConn)[2,1]

# sort out table with all coefficients #
row_order = c(3,2, 9, 4, 5, 6, 7, 10, 11, 8, 12)
coefs_tab = data.frame(varab = tidy(mod_RPD_conn35_full, conf.int = T)[row_order, 4])


vals = round(tidy(mod_RPD_conn35_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$conn35 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

vals = round(tidy(mod_RPD_conn32_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$conn32 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

vals = round(tidy(mod_RPD_net35_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$net35 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

vals = round(tidy(mod_RPD_net32_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$net32 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

tab_df(coefs_tab,
       file="suppTables/coefs_RPD_local_drivers")




#### visualisation ####

cols = c("darkslategrey", "darkolivegreen4", "darkgoldenrod")



## pred vs conn ##

range(df_mod$conn35)
conn_breaks = seq(0,  300000, 50000)
conn_breaks_scaled = (conn_breaks - attributes(df_mod$conn35_sc)$`scaled:center`)/
  attributes(df_mod$conn35_sc)$`scaled:scale`
conn_breaks = conn_breaks/10000


pred_conn_1 = visreg(mod_RPD_conn35_full, xvar = "conn35_sc", by= "pred_sc",
                     #breaks = c(min(df_mod$pred_sc), median(df_mod$pred_sc), max(df_mod$pred_sc)),
                     scale = "response", overlay = T, 
                     gg = T, nn = 1000, line.par = list(size = 1, lty = rep(c(1,2,6), each = 1000)), partial = F, rug = 2) +
  
  scale_colour_manual(values = cols, labels = c("low", "mid", "high"), name = "Predation pressure",   guide = guide_legend(
    override.aes = list(
      linetype = c("solid", "dashed", "twodash")
    ))) +
  scale_fill_manual(values = alpha(cols, 0.5), labels = c("low", "mid", "high"), name = "Predation pressure") +
  
  labs(x = "Connected area (ha)", y = "Probability of predator dominance") +
  
  scale_x_continuous(breaks = conn_breaks_scaled, labels = conn_breaks) +
  
  theme_bw(base_size = 7)+
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.20, 0.20, "cm"))


# save as figure source data
plot_data = pred_conn_1$data
plot_data$connectivity = plot_data$x*attributes(df_mod$conn35_sc)$`scaled:scale` + attributes(df_mod$conn35_sc)$`scaled:center` # back to original scale
plot_data$predation_pressure = as.numeric(plot_data$pred_sc)*attributes(df_mod$pred_sc)$`scaled:scale` + attributes(df_mod$pred_sc)$`scaled:center` # back to original scale

write.csv(plot_data, 
          "source-data-figures/source_predictions-RPD_fig2.csv",
          row.names = FALSE)


## temp vs dist ##

range(df_mod$DD)
temp_breaks = seq(100,  500, 100)
temp_breaks_scaled = (temp_breaks - attributes(df_mod$temp_sc)$`scaled:center`)/
  attributes(df_mod$temp_sc)$`scaled:scale`


temp_dist_1 = visreg(mod_RPD_conn35_full, xvar = "temp_sc", by = "distance_sc", 
                     #breaks = c(min(df_mod$temp_sc), median(df_mod$temp_sc), max(df_mod$temp_sc)),
                     scale = "response", overlay = T, 
                     gg = T, nn = 1000, line.par = list(size = 1, lty = rep(c(1,2,6), each = 1000)), partial = F, rug = 2) +
  
  
  scale_colour_manual(values = cols, 
                      labels = c("short", "mid", "long"), 
                      name = "Distance from open sea",   guide = guide_legend(
                        override.aes = list(
                          linetype = c("solid", "dashed", "twodash")
                        ))) +
  scale_fill_manual(values = alpha(cols, 0.5), 
                    labels = c("short", "mid", "long"), 
                    name = "Distance from open sea") +
  
  
  
  labs(x = "Temperature (degree days)", y = "Probability of predator dominance") +
  scale_x_continuous(breaks = temp_breaks_scaled, labels = temp_breaks) +
  
  theme_bw(base_size = 7) +
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.35, 0.20, "cm"))

# save as figure source data
plot_data = temp_dist_1$data
plot_data$DD = plot_data$x*attributes(df_mod$temp_sc)$`scaled:scale` + attributes(df_mod$temp_sc)$`scaled:center` # back to original scale
plot_data$distance = as.numeric(plot_data$distance_sc)*attributes(df_mod$distance_sc)$`scaled:scale` + attributes(df_mod$distance_sc)$`scaled:center` # back to original scale

write.csv(plot_data, 
          "source-data-figures/source_predictions-RPD_fig3.csv",
          row.names = FALSE)

## fish vs conn ##

fish_conn_1 = visreg(mod_RPD_conn35_full, xvar = "conn35_sc", by = "fishing_sc",
                     scale = "response", overlay = T, 
                     gg = T, nn = 1000, line.par = list(size = 1, lty = rep(c(1,2,6), each = 1000)), partial = F, rug = 2) +
  
  scale_colour_manual(values = cols, labels = c("low", "mid", "high"), name = "Fishing pressure",   guide = guide_legend(
    override.aes = list(
      linetype = c("solid", "dashed", "twodash")
    ))) +
  scale_fill_manual(values = alpha(cols, 0.5), labels = c("low", "mid", "high"), name = "Fishing pressure") +
  
  labs(x = "Connected area (ha)", y = "Probability of predator dominance") +
  
  scale_x_continuous(breaks = conn_breaks_scaled, labels = conn_breaks) +
  
  theme_bw(base_size = 9)+
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.20, 0.20, "cm"))


