library(tidyverse)
library(sf)
library(here)
library(tidymodels)

source("code/config.R")
#load(DATA_PATH)

# Import and join data -----------------------------------------------------------------

# Downloaded from UMESC 9/23/2025 with the option to download entire fish dataset
#https://umesc.usgs.gov/data_library/fisheries/LTRM_FISH_DATA_ENTIRE.zip
raw_fish <- read_csv(here("data", "raw", "ltrm_fish_data_entire.csv"))

#Import pool boundary data
pool_bounds <- st_read(here("data", "raw", "spatial", "UMRR_Pools", "UMRR_Pools.shp")) %>%
  st_transform(32615) %>%
  filter(POOL %in% c("Pool 04", "Pool 08", "Pool 13"))

# Aquatic Area boundary data
aqa_dat <- st_read(here("data", "raw","spatial", "aqa", "aqa_2010_lvl3_011918.shp")) %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -9999))) %>%
  select(OBJECTID:pool)

# Modify dates; add stratum names
fish <- raw_fish %>%
  mutate(sdate = as.Date(sdate, format = "%m/%d/%Y"),
         fdate = as.Date(fdate, format = "%m/%d/%Y"),
         year = format(sdate, "%Y"),
         stratum_full = as.factor(case_when(str_detect(stratum, "BWC") ~ "Backwater",
                                            str_detect(stratum, "MCB") ~ "Main Chann. Bord.",
                                            str_detect(stratum, "IMP") ~ "Impoundment",
                                            stratum == "SCB" ~ "Side Chann. Bord.",
                                            stratum == "SCB" ~ "Side Chann. Bord.",
                                            stratum == "TWZ" ~ "Tailwater",
                                            stratum == "UXO" ~ "Unexpl. Ord.",
                                            TRUE ~ "Other")))

# Subset data
fish_sub <- fish %>%
  filter(gear == "D",                   # Daylight electrofishing
         sitetype %in% c(0, 1),         # Stratified random sampling sites, not fixed sites
         pool %in% c("04", "08", "13"), # Upper pools
         year > 2009)                   # Post 2009

#Survey data; no fish
surveys <- fish_sub %>% 
  select(year, site:othrstrc, -barcode) %>%
  distinct() 

# Make spatial dataframe from survey data
surveys_sf <- surveys %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = st_crs(26915))

# Intersect survey and aquatic area data
surveys_aqa <- st_intersection(surveys_sf, aqa_dat) 

# Retain just the surveys found within CFL habitat; other habitat is side channel, main channel etc.
backwater_surveys <- surveys_aqa %>%
  filter(AQUA_DESC == "Contiguous Floodplain Lake") %>%
  mutate(utm_e = st_coordinates(.)[,1],
         utm_n = st_coordinates(.)[,2])

# Filter Aquatic Areas to just those with upper pool fish data
# 221 AqA with fish data
aqa_w_fish <- filter(aqa_dat, uniq_id %in% backwater_surveys$uniq_id)

# Finalize fish and covariate data --------------------------------------------

# Wide version of fish data
fish_wide <- backwater_surveys %>%
  st_drop_geometry() %>%
  select(year, sdate, utm_e, utm_n) %>% # These variables identify unique variables
  left_join(fish_sub) %>% 
  group_by(sdate, utm_e, utm_n, fishcode) %>%
  summarize(catch = sum(catch)) %>%
  pivot_wider(names_from = fishcode, values_from = catch, values_fill = 0) %>%
  select(-`NA`)

# Fish data with predictors
fish_surveys <- backwater_surveys %>% 
  left_join(fish_wide) %>%
  st_drop_geometry() %>%
  ungroup()

# Species names
species <- names(fish_surveys)[71:154]

# Presence-absence data wide
fish_binary <- fish_surveys %>% 
  mutate(across(any_of(species), ~ case_when(. > 0 ~ 1, TRUE ~ 0)))

# Long version of presence absence data
fish_long <- fish_binary %>%
  pivot_longer(cols = c(BKCP:IODR), names_to = "species", values_to = "occ")

# Prevalence of each taxa
# Retain species > 5% prevalent
occ_prop <- fish_long %>% 
  group_by(species) %>%
  summarize(prop = sum(occ)/nrow(fish_binary)) %>%
  filter(prop > .05,
         !str_detect(species, "U-")) %>% # Unidentified species
  arrange(desc(prop))

# Species list
sp_list <- occ_prop$species

# Filter to just species with > .05 occupancy
fish_common <- fish_long %>%
  filter(species %in% occ_prop$species) 

# Data for "Aquatic Area" based Cross Validation -----------------------------------------------

#  Spatial version of fish data
fish_sites <- fish_common %>%
  st_as_sf(coords = c("utm_e", "utm_n"), crs = st_crs(aqa_w_fish)) %>%
  filter(species == sp_list[1])

# Add number of sample points to aquatica area data
aqa_w_fish$n_pts <- lengths(st_intersects(aqa_w_fish, fish_sites))

# Get centroids of each AqA
centroids <- st_coordinates(st_centroid(aqa_w_fish))

# Identify connected components; 
# Need to identify because of the large spatial discontinuity between lakes
nb <- spdep::knearneigh(centroids, k = 4) %>% spdep::knn2nb() # Find nearest neighbors of each cluster
g <- igraph::graph_from_adj_list(nb) # Graph the neighbors
aqa_w_fish$component <- igraph::components(g)$membership # Find number of distinct components among the locations

# Get total number of surveys
total_pts <- sum(aqa_w_fish$n_pts)

# Get stats on components
region_allocation <- aqa_w_fish %>%
  dplyr::group_by(component) %>%
  dplyr::summarise(
    m_lakes = length(aa_num), # Lakes per component
    comp_pts = sum(n_pts), # Pts per component
    n_regions = round(comp_pts / total_pts * 25) # Approximate clusters per component
  )

ggplot(region_allocation)+
  geom_sf(aes(color = as.factor(component)))

# Add empty cluster column; start with cluster ID == 0
aqa_w_fish$cluster <- NA_integer_
cluster_id <- 0

set.seed(1)

for (cmp in unique(aqa_w_fish$component)) {
  
  sub <- aqa_w_fish[aqa_w_fish$component == cmp, ]
  
  total_pts_sub <- sum(sub$n_pts)
  
  # number of clusters for this component
  k_sub <- round(total_pts_sub / total_pts * 21)
  k_sub <- max(1, k_sub)
  
  # When there are only 1 component to choose from (e.g., upper pool 4)
  if (k_sub == 1) {
    aqa_w_fish$cluster[aqa_w_fish$component == cmp] <- cluster_id + 1
    cluster_id <- cluster_id + 1
    next
  }
  
  # polygon centroids
  coords_sub <- st_coordinates(st_centroid(sub))
  n <- nrow(coords_sub)
  
  # replicate centroids by point counts
  coords_rep <- coords_sub[rep(seq_len(n), sub$n_pts), , drop = FALSE]
  poly_id <- rep(seq_len(n), sub$n_pts)
  
  # k-means on replicated data
  km <- kmeans(coords_rep, centers = k_sub)
  
  # majority vote back to polygons
  clusters <- tapply(
    km$cluster,
    poly_id,
    function(x) as.integer(names(which.max(table(x))))
  )
  
  aqa_w_fish$cluster[aqa_w_fish$component == cmp] <- clusters + cluster_id
  cluster_id <- max(aqa_w_fish$cluster, na.rm = TRUE)
}

# What lakes are in which cluster?
# This data will be added to final data file before saving
backwater_cluster <- st_drop_geometry(aqa_w_fish) %>%
  select(aa_num, patch_cluster = cluster)

# Assess Clusters

cluster_counts <- aqa_w_fish %>%
  st_drop_geometry() %>%
  mutate(x = centroids[,1], y = centroids[,2]) %>%
  group_by(pool, component, cluster) %>%
  summarize(x = mean(x), y = mean(y),
            tot_pts = sum(n_pts),
            lakes_cluster = length(aa_num)) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(aqa_w_fish))

hist(cluster_counts$tot_pts)  

pool_plot <- "P04"

aqa_w_fish %>%
  filter(pool == pool_plot) %>%
  ggplot()+
  geom_sf(aes(color = as.factor(cluster)))+
  geom_sf_text(data = filter(cluster_counts, pool == "P04"), aes(label = tot_pts))+
  guides(color = "none")+
  labs(title = "Pool 4") +
  theme_void()

# Data for "Pool" based Cross Validation -----------------------------------------------

# Divide data into 6 spatial clusters for cross-validation
fish_6fold <- fish_common %>%
  mutate(cluster = case_when(pool == "04" & utm_n > 4930000 ~ "04u",
                             pool == "08" & utm_n > 4847000 & utm_e < 641000 ~ "08u",
                             pool == "13" & utm_n > 4666000 ~ "13u",
                             TRUE ~ pool),
         cluster_num = as.numeric(as.factor(cluster)))

# Plot to assess pool clusters
fish_6fold %>%
  filter(species == sp_list[1]) %>%
  ggplot()+
  geom_point(aes(utm_e, utm_n, fill = cluster), color = "black", pch=21)+
  facet_wrap(~paste("Pool", pool), scales = "free")+
  scale_y_continuous(n.breaks = 20)+
  scale_x_continuous(n.breaks = 10) + 
  labs(fill = "Cluster")+
  theme_void()+
  theme(strip.text = element_text(size = 12))

fish_6fold %>%
  filter(species == sp_list[1]) %>%
  pull(cluster) %>%
  table()

# Fish data with all spatial cv scenarios
fish_all_cv <- fish_6fold %>% left_join(backwater_cluster)

write_csv(fish_all_cv, "data/cleaned/fish_all_cv.csv")

#sp_list <- unique(fish_all_cv$species)
#if (!is.null(SP_SUBSET)) sp_list <- SP_SUBSET

# ── Pool-level spatial folds (per species, stored as list) ─────────────────
spatial_folds_list <- lapply(sp_list, function(sp) {
  fish_train <- fish_all_cv %>%
    filter(species == sp) %>%
    dplyr::select(all_of(VARS_REST)) %>%
    mutate(across(c(occ, period, cluster, vegd, substrt, snag), as.factor)) %>%
    mutate(occ = factor(occ, levels = c("1","0")))
  
  list(
    data   = fish_train,
    folds  = group_vfold_cv(fish_train, group = cluster, v = N_FOLDS)
  )
})
names(spatial_folds_list) <- sp_list

# ── Patch-cluster repeated folds ───────────────────────────────────────────
cluster_folds_list <- lapply(sp_list, function(sp) {
  fish_train <- fish_all_cv %>%
    filter(species == sp) %>%
    dplyr::select(-cluster, cluster = patch_cluster) %>%
    dplyr::select(all_of(VARS_REST)) %>%
    mutate(across(c(occ, period, cluster, vegd, substrt, snag), as.factor)) %>%
    mutate(occ = factor(occ, levels = c("1","0")))
  
  list(
    data  = fish_train,
    folds = group_vfold_cv(fish_train, group = cluster,
                           v = N_CLUSTERS / 2, repeats = N_REPEATS)
  )
})
names(cluster_folds_list) <- sp_list

# ── Random folds ───────────────────────────────────────────────────────────
random_folds_list <- lapply(sp_list, function(sp) {
  fish_train <- spatial_folds_list[[sp]]$data
  list(
    data  = fish_train,
    folds = vfold_cv(fish_train, v = N_FOLDS, strata = occ)
  )
})
names(random_folds_list) <- sp_list

save(sp_list,
     spatial_folds_list,
     cluster_folds_list,
     random_folds_list,
     file = paste0(OUTPUT_PATH, "models/01_folds.RData"))

cat("01_data_prep.R complete:", length(sp_list), "species,",
    "3 CV schemes\n")
