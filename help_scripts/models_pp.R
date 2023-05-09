
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### pike & perch: wave, spatio-temporal patterns, drivers of residuals ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#### create model data frame ####


# subset so that no values are NA
df_mod = df_sub[!is.na(df_sub$totalTopPred) & !is.na(df_sub$SWM) & !is.na(df_sub$fishPred) & !is.na(df_sub$BIASmean) & !is.na(df_sub$conn35),]
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
base_model_pp = sdmTMB(
  fishPred ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "off",
  
  family = nbinom2())

# estimates and AIC
tidy(base_model_pp, conf.int = T)
AIC(base_model_pp) 


# model evaluation
# mod = base_model_pp
# source("help_scripts/model_evaluation_sdmTMB.R")

# cross-validation
m_cv = sdmTMB_cv(
  fishPred ~

    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,

  data = df_mod,

  mesh = mesh,
  spatial = "off",

  fold_ids = clust,
  k_folds = length(unique(clust)),

  family = nbinom2())

m_cv$elpd # -0.6668057



#### BASIC RANDOM YEAR EFFECT ####

# fit model
temp_model_pp = sdmTMB(
  fishPred ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    (1 | yearF),
  
  data = df_mod,
  
  mesh = mesh,
  spatial = "off",
  
  family = nbinom2())

# estimates and AIC
tidy(temp_model_pp, conf.int = T)
AIC(temp_model_pp) 


# model evaluation
# mod = temp_model_pp
# source("help_scripts/model_evaluation_sdmTMB.R")


# cross-validation
m2_cv = sdmTMB_cv(
  fishPred ~

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

  family = nbinom2())

m2_cv$elpd # -0.6547103




#### CONSTANTS SPATIAL FIELDS ####


# fit (or load) model
spat_model_pp = sdmTMB(
  fishPred ~

    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,

  data = df_mod,

  mesh = mesh,
  spatial = "on",

  family = nbinom2())

save(spat_model_pp , file = "models/spat_model_pp.Rdata")

load("models/spat_model_pp.Rdata")


# estimates and AIC
tidy(spat_model_pp, conf.int = T) 
AIC(spat_model_pp) 



# model evaluation
# mod = spat_model_pp
# source("help_scripts/model_evaluation_sdmTMB.R")


# cross-validation
m3_cv = sdmTMB_cv(
  fishPred ~

    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,

  data = df_mod,

  mesh = mesh,
  spatial = "on",

  fold_ids = clust,
  k_folds = length(unique(clust)),

  family = nbinom2())


m3_cv$elpd # -0.6203497



#### SPATIAL FIELDS WITH CROSS-AREA RANDOM YEAR EFFECT ####


# fit (or load) model
spat_temp_model_pp_1 = sdmTMB(
  fishPred ~

    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    (1 | yearF),

  data = df_mod,

  mesh = mesh,
  spatial = "on",


  family = nbinom2())


save(spat_temp_model_pp_1, file = "models/spat_temp_model_pp_1.Rdata")
load("models/spat_temp_model_pp_1.Rdata")


# estimates and AIC
tidy(spat_temp_model_pp_1, conf.int = T)
AIC(spat_temp_model_pp_1) 


# model evaluation
# mod = spat_temp_model_pp_1 
# source("help_scripts/model_evaluation_sdmTMB.R")



# cross-validation

plan(multisession, workers = 10)

m4_cv = sdmTMB_cv(
  fishPred ~

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

  family = nbinom2())

plan(sequential)

m4_cv$elpd





#### SPATIO-TEMPORAL FIELDS ####


# fit (or load) model
spat_temp_model_pp_2 = sdmTMB(
  fishPred ~

    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,

  data = df_mod,

  time = "year",
  spatiotemporal = "iid",

  mesh = mesh,
  spatial = "off",


  family = nbinom2())


save(spat_temp_model_pp_2 , file = "models/spat_temp_model_pp_2.Rdata")

load("models/spat_temp_model_pp_2.Rdata")

# estimates and AIC
tidy(spat_temp_model_pp_2, conf.int = T)
AIC(spat_temp_model_pp_2) # 22792.02




# model evaluation
# mod = spat_temp_model_pp_2 
# source("help_scripts/model_evaluation_sdmTMB.R")

# cross-validation
plan(multisession, workers = 10)

m5_cv = sdmTMB_cv(
  fishPred ~

    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc,

  data = df_mod,

  time = "year",
  spatiotemporal = "iid",

  mesh = mesh,
  spatial = "off",

  fold_ids = clust,
  k_folds = length(unique(clust)),

  family = nbinom2())

plan(sequential)

m5_cv$elpd





#### saving AICs and elpds ####

pp_spat_temp_model_comp = data.frame(
  pp_AIC = c(AIC(base_model) - AIC(base_model, temp_model, spat_model, spat_temp_model_1, spat_temp_model_2)),
  pp_elpds = c(m_cv$elpd, m2_cv$elpd, m43_cv$elpd, m4_cv$elpd, m5_cv$elpd)

)





#### ADD IN POTENTIAL DRIVERS ####

mod_ref = glmmTMB(
  fishPred ~ 
    
    BIAS_sc +
    distance_sc +
    BIAS_sc*distance_sc  +
    swm_sc +
    
    (1 | year),
  
  data = df_mod,
  
  family = nbinom2)

mod = mod_ref
#source("help_scripts/model_evaluation.R")


mod_ref2 = glmmTMB(
  fishPred ~ 
    
    conn35_sc +
    pred_sc +
    fishing_sc +
    conn35_sc*pred_sc +
    conn35_sc*fishing_sc +
    
    temp_sc +
    temp_sc*distance_sc +
    
    (1 | year),
  
  data = df_mod,
  
  family = nbinom2)

mod = mod_ref2
#source("help_scripts/model_evaluation.R")





#### conn 35  ####


mod_fishPred_conn35_full = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

mod = mod_fishPred_conn35_full
#source("help_scripts/model_evaluation.R")

# dredge_fishPred_conn35_full = dredge(mod_fishPred_conn35_full, trace = 2)
# tab_df(dredge_fishPred_conn35_full,
#        file="result_tables/dredge_fishPred_conn35_full.doc")

mod_fishPred_conn35_full_noFish = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

tidy(mod_fishPred_conn35_full_noFish, conf.int = T)
tidy(mod_fishPred_conn35_full, conf.int = T)[c(1:6,8:10,12:13),]



#### conn 32  ####


mod_fishPred_conn32_full = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

mod = mod_fishPred_conn32_full
#source("help_scripts/model_evaluation.R")

# dredge_fishPred_conn32_full = dredge(mod_fishPred_conn32_full, trace = 2)
# tab_df(dredge_fishPred_conn32_full,
#        file="result_tables/dredge_fishPred_conn32_full.doc")



mod_fishPred_conn32_full_noFish = glmmTMB(
  fishPred ~ 
    
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
  
  
  family = nbinom2)

tidy(mod_fishPred_conn32_full_noFish, conf.int = T)
tidy(mod_fishPred_conn32_full, conf.int = T)[c(1:6,8:10,12:13),]



#### net 35  ####


mod_fishPred_net35_full = glmmTMB(
  fishPred ~ 
    
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
  
  data = df_mod,
  
  na.action = "na.fail",
  
  family = nbinom2)

mod = mod_fishPred_net35_full
#source("help_scripts/model_evaluation.R")

# dredge_fishPred_net35_full = dredge(mod_fishPred_net35_full, trace = 2)
# tab_df(dredge_fishPred_net35_full,
#        file="result_tables/dredge_fishPred_net35_full.doc")


mod_fishPred_net35_full_noFish = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

tidy(mod_fishPred_net35_full_noFish, conf.int = T)
tidy(mod_fishPred_net35_full, conf.int = T)[c(1:6,8:10,12:13),]




#### net 32  ####


mod_fishPred_net32_full = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

mod = mod_fishPred_net32_full
#source("help_scripts/model_evaluation.R")

# dredge_fishPred_net32_full = dredge(mod_fishPred_net32_full, trace = 2)
# tab_df(dredge_fishPred_net32_full,
#        file="result_tables/dredge_fishPred_net32_full.doc")



mod_fishPred_net32_full_noFish = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

tidy(mod_fishPred_net32_full_noFish, conf.int = T)
tidy(mod_fishPred_net32_full, conf.int = T)[c(1:6,8:10,12:13),]



# SAVE ESTIMATES AND COEFFICIENTS #

# AIC and R2 table #

mod_sel_tab = AIC(mod_ref, mod_ref2, mod_fishPred_conn35_full, mod_fishPred_conn32_full, mod_fishPred_net35_full, mod_fishPred_net32_full)
mod_sel_tab$R2 = c(
  r.squaredGLMM(mod_ref)[2,1],
  r.squaredGLMM(mod_ref2)[2,1],
  r.squaredGLMM(mod_fishPred_conn35_full)[2,1],
  r.squaredGLMM(mod_fishPred_conn32_full)[2,1],
  r.squaredGLMM(mod_fishPred_net35_full)[2,1],
  r.squaredGLMM(mod_fishPred_net32_full)[2,1])
mod_sel_tab$dAIC = mod_sel_tab$AIC-mod_sel_tab$AIC[1]
mod_sel_tab = mod_sel_tab[, c(3:4)]

tab_df(mod_sel_tab,
       file="result_tables/mod_sel_fishPred_local_drivers.doc")



# sort out table with all coefficients #
row_order = c(3,2, 9, 4, 5, 6, 7, 10, 11, 8, 12)
coefs_tab = data.frame(varab = tidy(mod_fishPred_conn35_full, conf.int = T)[row_order, 4])


vals = round(tidy(mod_fishPred_conn35_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$conn35 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

vals = round(tidy(mod_fishPred_conn32_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$conn32 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

vals = round(tidy(mod_fishPred_net35_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$net35 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

vals = round(tidy(mod_fishPred_net32_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$net32 =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

tab_df(coefs_tab,
       file="result_tables/coefs_fishPred_local_drivers")





#### visualisation ####


cols = c("darkslategrey", "darkolivegreen4", "darkgoldenrod")

## pred vs conn ##

range(df_mod$conn35)
conn_breaks = seq(0,  300000, 50000)
conn_breaks_scaled = (conn_breaks - attributes(df_mod$conn35_sc)$`scaled:center`)/
  attributes(df_mod$conn35_sc)$`scaled:scale`
conn_breaks = conn_breaks/10000


pred_conn_2 = visreg(mod_fishPred_conn35_full, xvar = "conn35_sc", by = "pred_sc", 
                     #breaks = c(min(df_mod$pred_sc), median(df_mod$pred_sc), max(df_mod$pred_sc)),
                     scale = "response", overlay = T, 
                     gg = T, nn = 1000, line.par = list(size = 1, lty = rep(c(1,2,6), each = 1000)), partial = F, rug = 2) +
  
  scale_colour_manual(values = cols, labels = c("low", "mid", "high"), name = "Predation pressure",   guide = guide_legend(
    override.aes = list(
      linetype = c("solid", "dashed", "twodash")
    ))) +
  scale_fill_manual(values = alpha(cols, 0.5), labels = c("low", "mid", "high"), name = "Predation pressure") +
  
  labs(x = "Connected area (ha)", y = "Predator density") +
  
  scale_x_continuous(breaks = conn_breaks_scaled, labels = conn_breaks) +
  
  ylim(0,56.5) +
  
  theme_bw(base_size = 9)+
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.20, 0.20, "cm"))



## temp vs dist ##

range(df_mod$DD)
temp_breaks = seq(100,  500, 100)
temp_breaks_scaled = (temp_breaks - attributes(df_mod$temp_sc)$`scaled:center`)/
  attributes(df_mod$temp_sc)$`scaled:scale`


temp_dist_2 = visreg(mod_fishPred_conn35_full, 
                     xvar = "temp_sc", by = "distance_sc", 
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
  
  labs(x = "Temperature (degree days)", y = "Predator density") +

  scale_x_continuous(breaks = temp_breaks_scaled, labels = temp_breaks) +
  
  ylim(0,60) +
  
  theme_bw(base_size = 9) +
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.20, 0.20, "cm"))


## fishing vs conn ##


fish_conn_2 = visreg(mod_fishPred_conn35_full, xvar = "conn35_sc", by = "fishing_sc", 
       scale = "response", overlay = T, 
       gg = T, nn = 1000, line.par = list(size = 1, lty = rep(c(1,2,6), each = 1000)), partial = F, rug = 2) +
  
  scale_colour_manual(values = cols, labels = c("low", "mid", "high"), name = "Fishing pressure",   guide = guide_legend(
    override.aes = list(
      linetype = c("solid", "dashed", "twodash")
    ))) +
  scale_fill_manual(values = alpha(cols, 0.5), labels = c("low", "mid", "high"), name = "Fishing pressure") +
  
  labs(x = "Connected area (ha)", y = "Predator density") +
  
  scale_x_continuous(breaks = conn_breaks_scaled, labels = conn_breaks) +
  
  ylim(0,17) +
  
  theme_bw(base_size = 9)+
  theme_sets +
  theme(
    plot.margin = margin(0.40, 0.20, 0.20, 0.20, "cm"))





#### seals vs cormorants ####

mod_pp_conn35_corm = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

summary(mod_pp_conn35_corm)
summary(mod_fishPred_conn35_full)
AIC(mod_pp_conn35_corm)
AIC(mod_fishPred_conn35_full)
r.squaredGLMM(mod_pp_conn35_corm)
r.squaredGLMM(mod_fishPred_conn35_full)


mod_pp_conn35_seal = glmmTMB(
  fishPred ~ 
    
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
  
  family = nbinom2)

summary(mod_pp_conn35_seal)
summary(mod_fishPred_conn35_full)
AIC(mod_pp_conn35_seal)
AIC(mod_fishPred_conn35_full)
r.squaredGLMM(mod_pp_conn35_seal)
r.squaredGLMM(mod_fishPred_conn35_full)


row_order = c(3,2, 9, 4, 5, 6, 7, 10, 11, 8, 12)
coefs_tab = data.frame(varab = tidy(mod_pp_conn35_full, conf.int = T)[row_order, 4])
vals = round(tidy(mod_pp_conn35_full, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$both =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )
vals = round(tidy(mod_pp_conn35_corm, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$corm =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )
vals = round(tidy(mod_pp_conn35_seal, conf.int = T)[row_order, c(5, 9, 10)], digits = 2)
coefs_tab$seal =  paste0(vals$estimate, " (", vals$conf.low, ";", vals$conf.high,")" )

coefs_tab = rbind(coefs_tab, c("AIC", 
                               round(AIC(mod_pp_conn35_full)),
                               round(AIC(mod_pp_conn35_corm)),
                               round(AIC(mod_pp_conn35_seal))
))

coefs_tab = rbind(coefs_tab, c("R2", 
                               round(r.squaredGLMM(mod_pp_conn35_full)[2,1], digits = 2),
                               round(r.squaredGLMM(mod_pp_conn35_corm)[2,1], digits = 2),
                               round(r.squaredGLMM(mod_pp_conn35_seal)[2,1], digits = 2)
))



tab_df(coefs_tab,
       file="result_tables/sealsVScormsPP.doc")


