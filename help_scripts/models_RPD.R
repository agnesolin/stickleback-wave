
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### RPD: wave, spatio-temporal patterns, drivers of residuals ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#### create model data frame ####


# subset so that no values are NA
df_mod = df_sub[!is.na(df_sub$totalTopPred) & !is.na(df_sub$SWM) & !is.na(df_sub$RPD) & !is.na(df_sub$BIASmean) & !is.na(df_sub$conn35),]
colSums(is.na(df_mod))
nrow(df_mod)

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
# source("help_scripts/model_evaluation_sdmTMB.R")


# cross-validation
m3_cv = sdmTMB_cv(
  RPD ~
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "on",
  
  fold_ids = clust,
  k_folds = length(unique(clust)),
  
  family = binomial())


m3_cv$elpd # -0.3726271




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
# source("help_scripts/model_evaluation_sdmTMB.R")


#cross-validation
plan(multisession, workers = 5)

m4_cv = sdmTMB_cv(
  RPD ~
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    (1 | yearF),
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "on",
  
  
  fold_ids = clust,
  k_folds = length(unique(clust)),
  
  
  family = binomial())

plan(sequential)

m4_cv$elpd # -0.3376445



#### SPATIO-TEMPORAL FIELDS ####


# fit (or load) model
spat_temp_model_2 = sdmTMB(
  RPD ~
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,
  
  data = df_mod,
  
  time = "year",
  spatiotemporal = "iid",
  
  mesh = mesh,
  spatial = "off",
  
  
  family = binomial(),
  weights = rep(1, nrow(df_mod)))


save(spat_temp_model_2 , file = "models/spat_temp_model_2.Rdata")

load("models/spat_temp_model_2.Rdata")

# estimates and AIC
tidy(spat_temp_model_2, conf.int = T)
AIC(spat_temp_model_2) # 2581.013

# model evaluation
# mod = spat_temp_model_2 
# source("help_scripts/model_evaluation_sdmTMB.R")

# cross-validation
plan(multisession, workers = 5)

m5_cv = sdmTMB_cv(
  RPD ~
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,
  
  data = df_mod,
  
  time = "year",
  spatiotemporal = "iid",
  
  mesh = mesh,
  spatial = "off",
  
  parallel = TRUE,
  
  fold_ids = clust,
  k_folds = length(unique(clust)),
  
  family = binomial())

plan(sequential)

m5_cv$elpd




#### saving AICs and elpds ####

RPD_spat_temp_model_comp = data.frame(
  RPD_AIC = c(AIC(base_model) - AIC(base_model, temp_model, spat_model, spat_temp_model_1, spat_temp_model_2)),
  RPD_elpds = c(m_cv$elpd, m2_cv$elpd, m3_cv$elpd, m4_cv$elpd,NA) # m5_cv$elpd)
  
)




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
    #plot.margin = margin(0.5,0.3,0.3,0.5, "cm"),
    legend.key.width = unit(0.7, 'cm') 
  ) 

# annotation_raster(imgS, ymin = -0.12, ymax = -0.02, xmin = -2.15, xmax = -1.55) +
# annotation_raster(imgP, ymin = 1.02, ymax = 1.12, xmin = -2.25, xmax = -1.5)



ggsave("suppFigures/stickleback_wave.png", width = 180, height = 9, units = "cm")

# 
# #### random effects ####
# 
# # make a grid based on where data are collected and predict based on this
# 
# 
# # coastline
# coast = st_read("data/Europe_coastline_poly.shp")
# geom_coast = st_geometry(coast)
# 
# # create a grid for prediction 
# r = raster(xmn = min(df_mod$X_new),
#            xmx = max(df_mod$X_new), 
#            ymn = min(df_mod$Y_new), 
#            ymx = max(df_mod$Y_new), 
#            res = 5000)
# 
# 
# crs(r) = crs(land5)
# 
# pts = SpatialPointsDataFrame(coords = cbind(df_mod$X_new, df_mod$Y_new), data = data.frame(year = df_mod$year), proj4string = crs(geom_coast))
# 
# 
# countr = rasterize(pts, r, "year", fun=function(x,...){length(unique(x))})
# countr =  as.data.frame(rasterToPoints(countr))
# 
# countr = countr[countr$layer > 1 ,] # only include cells with data from at least 2 years
# 
# # predict
# nd = data.frame(X_new = countr$x, Y_new = countr$y, 
#                 BIAS_sc = mean(df_mod$BIAS_sc), distance_sc = mean(df_mod$distance_sc), swm_sc = mean(df_mod$swm_sc))
# 
# p = predict(spat_model, newdata = nd)
# 
# # sort out projection
# C = SpatialPoints(
#   coords = cbind(p$X_new, p$Y_new),
#   proj4string = crs(land25)
# )
# C2 = spTransform(
#   C,
#   
#   CRSobj =  crs(geom_coast)
# )
# p$X_map = coordinates(C2)[, 1]
# p$Y_map = coordinates(C2)[, 2]
# 
# xlims = range(coordinates(C2)[, 1])
# ylims = range(coordinates(C2)[, 2])
# 
# 
# 
# 
# ggplot(data = geom_coast) + # coastline
#   geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
#   coord_sf(xlim = xlims, ylim = ylims) +
#   
#   geom_tile(data = p, 
#             aes(x = X_map, y = Y_map, fill = omega_s), width = 5000, height = 5000) +
#   scale_fill_continuous_divergingx(palette = "Fall", mid = 0, rev = TRUE, name = "Residual variation") + 
#   #  annotate("text", xlims[1] + 0.2*(xlims[2]-xlims[1]), ylims[1] + 0.08*(ylims[2]-ylims[1]), label = "Blekinge", family = "sans", size = 13) +
#   
#   #  annotate("text", xlims[1] + 0.3*(xlims[2]-xlims[1]), ylims[1] + 0.2*(ylims[2]-ylims[1]), label = "Kalmarsund", family = "sans", size = 13) +
#   
#   #  annotate("text", xlims[1] + 0.38*(xlims[2]-xlims[1]), ylims[1] + 0.57*(ylims[2]-ylims[1]), label = "St Anna", family = "sans", size = 13) +
#   
#   #  annotate("text", xlims[1] + 0.55*(xlims[2]-xlims[1]), ylims[1] + 0.97*(ylims[2]-ylims[1]), label = "Forsmark", family = "sans", size = 13) +
#   
#   theme_bw() +
#   theme_sets +
#   theme(
#     legend.position = "bottom",
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_line(),
#     panel.grid.minor = element_line(),
#   )
# 
# 
# ggsave("figures/resid_spat_Johan.jpg", width = 18.5, height = 22, units = "cm")
# 
# 
# 
# #### yearly maps ####
# 
# # predict on a yearly basis
# p = predict(spat_temp_model_2)
# 
# # sort out projection
# C = SpatialPoints(
#   coords = cbind(p$X_new, p$Y_new),
#   proj4string = crs(land25)
# )
# C2 = spTransform(
#   C,
#   
#   CRSobj =  crs(geom_coast)
# )
# 
# p$X_map = coordinates(C2)[, 1]
# p$Y_map = coordinates(C2)[, 2]
# 
# xlims = range(coordinates(C2)[, 1])
# ylims = range(coordinates(C2)[, 2])
# 
# 
# library(colorspace)
# 
# ggplot(data = geom_coast) + # coastline
#   geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
#   coord_sf(xlim = xlims, ylim = ylims) +
#   
#   geom_point(data = p, 
#              aes(x = X_map, y = Y_map, colour = epsilon_st)) +
#   scale_colour_continuous_divergingx(palette = "Fall", mid = 0, rev = TRUE, name = "Residual variation") + 
#   
#   facet_wrap(~year) +
#   
#   theme_bw() +
#   theme_sets +
#   theme(
#     legend.position = "bottom",
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     panel.grid.major = element_line(),
#     panel.grid.minor = element_line(),
#   )
# 
# 
# ggsave("suppFigures/resid_spat_year.jpg", width = 18.5, height = 22, units = "cm")
# 



#### ADD IN POTENTIAL DRIVERS ####



mod_ref = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    (1 | year),
  
  data = df_mod,
  
  family = binomial)

mod = mod_ref
#source("help_scripts/model_evaluation.R")


mod_ref2 = glmmTMB(
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

mod = mod_ref2
#source("help_scripts/model_evaluation.R")

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
#source("help_scripts/model_evaluation.R")

# dredge_RPD_conn35_full = dredge(mod_RPD_conn35_full, trace = 2)
# tab_df(dredge_RPD_conn35_full,
#        file="result_tables/dredge_conn35_full.doc")


mod_RPD_conn35_full_noFish = glmmTMB(
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

tidy(mod_RPD_conn35_full_noFish, conf.int = T)
tidy(mod_RPD_conn35_full, conf.int = T)[c(1:6,8:10,12:13),]



mod_RPD_conn35_PRED = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    conn35_sc +
    seal_sc +
    corm_sc +
    fishing_sc +
    conn35_sc*corm_sc +
    conn35_sc*seal_sc +
    conn35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)

tidy(mod_RPD_conn35_PRED, conf.int = T)


# dredge_RPD_conn35_PRED = dredge(mod_RPD_conn35_PRED, trace = 2)
# tab_df(dredge_RPD_conn35_PRED,
#        file="result_tables/dredge_conn35_PRED.doc")



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

mod = mod_RPD_conn32_full
#source("help_scripts/model_evaluation.R")

# dredge_RPD_conn32_full = dredge(mod_RPD_conn32_full, trace = 2)
# tab_df(dredge_RPD_conn32_full,
#        file="result_tables/dredge_conn32_full.doc")


mod_RPD_conn32_full_noFish = glmmTMB(
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

tidy(mod_RPD_conn32_full_noFish, conf.int = T)
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

mod = mod_RPD_net35_full
#source("help_scripts/model_evaluation.R")

# dredge_RPD_net35_full = dredge(mod_RPD_net35_full, trace = 2)
# tab_df(dredge_RPD_net35_full,
#        file="result_tables/dredge_net35_full.doc")


mod_RPD_net35_full_noFish = glmmTMB(
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

tidy(mod_RPD_net35_full_noFish, conf.int = T)
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

mod = mod_RPD_net32_full
#source("help_scripts/model_evaluation.R")


# dredge_RPD_net32_full = dredge(mod_RPD_net32_full, trace = 2)
# tab_df(dredge_RPD_net32_full,
#        file="result_tables/dredge_net32_full.doc")


mod_RPD_net32_full_noFish = glmmTMB(
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

tidy(mod_RPD_net32_full_noFish, conf.int = T)
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
       file="result_tables/mod_sel_RPD_local_drivers.doc")



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
       file="result_tables/coefs_RPD_local_drivers")




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
  
  theme_bw(base_size = 9)+
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.20, 0.20, "cm"))




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
  
  theme_bw(base_size = 9) +
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.35, 0.20, "cm"))



## fish vs conn ##

fish_conn_1 = visreg(mod_RPD_conn35_full, xvar = "conn35_sc", by= "fishing_sc",
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









#### residual map ####

df_mod$resids = residuals(mod_RPD_conn35_full, type = "response")

C = SpatialPoints(
  coords = cbind(df_mod$X, df_mod$Y),
  proj4string = crs(land25)
)
C2 = spTransform(
  C,
  
  CRSobj =  crs(geom_coast)
)
df_mod$X_map = coordinates(C2)[, 1]
df_mod$Y_map = coordinates(C2)[, 2]

xlims = range(df_mod$X_map)
ylims = range(df_mod$Y_map)


ggplot(data = geom_coast) + # coastline
  geom_sf(fill = "lightgrey" , colour = "lightgrey" , size = 0.2) +
  coord_sf(xlim = xlims, ylim = ylims) +
  
  geom_point(data = df_mod, aes(x = X_map, y = Y_map, colour = resids)) + 
  
  scale_colour_continuous_divergingx(palette = "Fall", mid = 0, rev = TRUE, name = "Residuals") +
  
  theme_bw() +
  theme_sets +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_line(),
    panel.grid.minor = element_line(),
  )

ggsave("suppFigures/map_resid_RPD.jpg", width = 18.5, height = 24, units = "cm")



#### seals vs cormorants ####

mod_RPD_conn35_corm = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    conn35_sc +
    corm_sc +
    fishing_sc +
    conn35_sc*corm_sc +
    conn35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    
    (1 | year),
  
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)

summary(mod_RPD_conn35_corm)
summary(mod_RPD_conn35_full)
AIC(mod_RPD_conn35_corm)
AIC(mod_RPD_conn35_full)
r.squaredGLMM(mod_RPD_conn35_corm)
r.squaredGLMM(mod_RPD_conn35_full)


mod_RPD_conn35_seal = glmmTMB(
  RPD ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    conn35_sc +
    seal_sc +
    fishing_sc +
    conn35_sc*seal_sc +
    conn35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    
    (1 | year),
  
  
  na.action = "na.fail",
  
  data = df_mod,
  
  family = binomial)

summary(mod_RPD_conn35_seal)
summary(mod_RPD_conn35_full)
AIC(mod_RPD_conn35_seal)
AIC(mod_RPD_conn35_full)
r.squaredGLMM(mod_RPD_conn35_seal)
r.squaredGLMM(mod_RPD_conn35_full)


row_order = c(3,2, 9, 4, 5, 6, 7, 10, 11, 8, 12)
coefs_tab = data.frame(varab = tidy(mod_RPD_conn35_full, conf.int = T)[row_order, 4])
vals = round(tidy(mod_RPD_conn35_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$both =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )
vals = round(tidy(mod_RPD_conn35_corm, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$corm =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )
vals = round(tidy(mod_RPD_conn35_seal, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$seal =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

coefs_tab = rbind(coefs_tab, c("AIC", 
                               round(AIC(mod_RPD_conn35_full)),
                               round(AIC(mod_RPD_conn35_corm)),
                               round(AIC(mod_RPD_conn35_seal))
))

coefs_tab = rbind(coefs_tab, c("R2", 
                               round(r.squaredGLMM(mod_RPD_conn35_full)[2,1], digits = 2),
                               round(r.squaredGLMM(mod_RPD_conn35_corm)[2,1], digits = 2),
                               round(r.squaredGLMM(mod_RPD_conn35_seal)[2,1], digits = 2)
))



tab_df(coefs_tab,
       file="result_tables/sealsVScorms.doc")


