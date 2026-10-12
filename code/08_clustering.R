library(tidyverse)
library(here)
library(purrr)
library(ggrepel)
library(ggdendro)
library(cluster) # For silhouette function
#library(dendextend)
library(fpc) # For clusterboot
library(patchwork)

source("code/config.R")
load(paste0(OUTPUT_PATH, "results/04_importance_pdp_pool_based.RData"))

#── Variable Importance Data ──────────────────────────────────
sp_guilds <- read_csv(here("data", "sp_guilds.csv")) %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(Common_Name = str_to_title(Common_Name))

# Variable importance
imp_rank <- importance_df_all %>%
  left_join(sp_guilds) %>%
  group_by(species) %>%
  mutate(rank = rank(desc(importance)),
         var_name = case_when(variable == "current" ~ "Current\nvelocity",
                              variable == "temp" ~ "Temperature",
                              variable == "depth" ~ "Depth",
                              variable == "period" ~ "Period",
                              variable == "do" ~ "Dissolved\noxygen",
                              variable == "substrt" ~ "Substrate",
                              variable == "snag" ~ "Snag",
                              variable == "vegd" ~ "Vegetation\ndensity"))

imp_rank %>% group_by(var_name) %>% summarize(mn_rank = mean(rank))

outliers <- imp_rank %>%
  group_by(variable, rank) %>%
  mutate(Rank_Count = n(),
         rank_Pct = Rank_Count / 33) %>% # 33 total species
  ungroup() %>%
  filter(rank_Pct < 0.10,
         !(species %in% c("GDSN", "BKCP", "CNCF") & variable == "do"),
         !(species %in% c("PNSD", "WDSN", "NTPK") & variable == "substrt"),
         !(species %in% c("CNCF", "JYDR", "YLBH") & variable == "current"),
         !(species %in% c("BLGL", "YWPH", "RKBS") & variable == "temp"),
         !(species %in% c("LGPH", "GNSF") & variable == "vegd")
         ) %>%
  dplyr::select(variable, var_name, Common_Name, rank, rank_Pct)

ggplot(imp_rank, aes(rank, reorder(var_name, desc(rank))))+
  ggridges::geom_density_ridges(aes(height = after_stat(density)), scale = .95, stat = "density", trim = TRUE, )+
  geom_point(data = outliers)+
  geom_text_repel(data = outliers, aes(label = Common_Name),
                  size = 2.75,
                  force = 3,
                  force_pull = 3,
                  nudge_y = 0.2,       # Starts labels 0.5 units above the point
                  box.padding = unit(0.4, "lines"))+ # More padding specifically on Y) +
  labs(x = "Variable importance rank", y = "")+
  scale_x_continuous(breaks = 1:8)+
  theme_classic(base_size = 12)

ggsave("output/figures/fig3_var_imp.jpeg", height = 5, width = 6)

#── Clustering helper functions ──────────────────────────────────

extract_categorical_features <- function(pd_df, var_name) {
  levels_df <- pd_df %>%
    arrange(.data[[var_name]]) %>%
    mutate(level = as.character(.data[[var_name]]))
  
  y <- levels_df$yhat
  n <- length(y)
  
  if (n < 2) return(tibble(
    cor_linear = NA_real_, 
    peak_location = NA_real_,
    #trough_location = NA_real_, 
    #response_range = NA_real_
  ))
  
  tibble(
    cor_linear      = cor(seq_along(y), y, method = "spearman"),
    peak_location   = (which.max(y) - 1) / (n - 1),
    trough_location = (which.min(y) - 1) / (n - 1),
    #response_range  = max(y) - min(y)
  )
}

extract_continuous_features <- function(pd_df, var_name) {
  x <- pd_df[[var_name]]
  y <- pd_df$yhat
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  
  if (length(unique(x)) < 5) return(tibble(
    cor_linear = NA_real_, peak_location = NA_real_,
    trough_location = NA_real_, response_range = NA_real_
  ))
  
  ord <- order(x); x <- x[ord]; y <- y[ord]
  x_norm <- (x - min(x)) / (max(x) - min(x))
  
  tibble(
    cor_linear      = cor(x, y, method = "spearman"),
    peak_location   = x_norm[which.max(y)],
    max_x = max(x),
    min_x = min(x),
    trough_location = x_norm[which.min(y)],
  )
}

extract_pdp_shape <- function(pd_df, group_keys) {
  var_name <- group_keys$variable
  if (var_name %in% CATEGORICAL_VARS) {
    extract_categorical_features(pd_df, var_name)
  } else {
    extract_continuous_features(pd_df, var_name)
  }
}

# ── Extract unweighted shape features ──────────────────────────────────

pdp_shape_features <- pdp_all %>%
  group_by(species, variable) %>%
  group_modify(extract_pdp_shape) %>%
  ungroup()

# ── Scaled importance (relative within species) ───────────────────────────────

importance_scaled_df <- importance_df_all %>%
  group_by(species) %>%
  mutate(importance_scaled = importance / sum(importance)) %>%
  ungroup() %>%
  dplyr::select(species, variable, importance_scaled)

shape_with_importance <- pdp_shape_features %>%
  left_join(importance_scaled_df, by = c("species", "variable"))

# ════════════════════════════════════════════════════════════════════════════
# Custom pairwise distance with species-specific importance weighting
# ════════════════════════════════════════════════════════════════════════════

# Normalize each shape feature to a comparable 0-1 scale across species,
# computed per-variable (Gower-style range normalization)
shape_normalized <- shape_with_importance %>%
  group_by(variable) %>%
  mutate(
    cor_linear_n = {
      rng <- range(cor_linear, na.rm = TRUE)
      if (diff(rng) < 1e-9) 0.5 else (cor_linear - rng[1]) / diff(rng) # Subtract minimum value; divide by range
    },
    peak_location_n   = peak_location,    # already 0-1
    trough_location_n = trough_location,  # already 0-1
  ) %>%
  ungroup() %>%
  mutate(across(c(cor_linear_n, 
                  peak_location_n, 
                  ),
                ~ replace_na(.x, 0.5)))  # neutral fill for missing values

species_list <- sort(unique(shape_normalized$species))
n_sp          <- length(species_list)

dist_mat_custom <- matrix(0, n_sp, n_sp,
                          dimnames = list(species_list, species_list))

get_sp_data <- function(sp) shape_normalized %>% filter(species == sp)
sp_data_cache <- setNames(lapply(species_list, get_sp_data), species_list)

for (i in seq_len(n_sp - 1)) {
  for (j in seq(i + 1, n_sp)) {
    
    sp_i <- species_list[i]
    sp_j <- species_list[j]
    
    merged <- inner_join(
      sp_data_cache[[sp_i]], sp_data_cache[[sp_j]],
      by = "variable", suffix = c("_i", "_j")
    )
    
    merged <- merged %>%
      mutate(
        d_cor    = abs(cor_linear_n_i      - cor_linear_n_j),
        d_peak   = abs(peak_location_n_i   - peak_location_n_j),
        d_trough = abs(trough_location_n_i - trough_location_n_j),
        d_var    = case_when(variable %in% c("current", "depth", "do", "temp", "vegd") ~ (d_cor + d_peak #)/ 2, 
                                                                                          + d_trough)/ 3,
                             TRUE ~ (d_peak + d_trough)/2),
        w_var    = (importance_scaled_i + importance_scaled_j) / 2
      )
    
    # Weighted average distance across predictors
    d_ij <- sum(merged$d_var * merged$w_var) / sum(merged$w_var)
    
    dist_mat_custom[i, j] <- d_ij
    dist_mat_custom[j, i] <- d_ij
  }
}

dist_custom    <- as.dist(dist_mat_custom)
hc_ward_custom <- hclust(dist_custom, method = "ward.D2")

# Change names to match names in subsequent code
hc_ward_gower <- hc_ward_custom
gower_dist <- dist_custom

sil_df <- map_dfr(2:8, function(k) {
  cuts   <- dendextend::cutree(hc_ward_gower, k = k)
  sil    <- cluster::silhouette(cuts, gower_dist)
  tibble(
    k         = k,
    avg_sil   = mean(sil[, "sil_width"]),
    median_sil = median(sil[, "sil_width"]), 
    min_sil   = min(sil[, "sil_width"]),
    n_negative = sum(sil[, "sil_width"] < 0)
  )
})

# ── Convert dist object to a full matrix ────────────────────────────────

dist_for_boot <- dist_custom
dist_matrix   <- as.matrix(dist_for_boot)

# Confirm matrix symmetry
stopifnot(isSymmetric(dist_matrix))

# ── Sweep across k using the distance-aware clustering interface ──────────────
boot_results <- map_dfr(2:8, function(k) {
  boot_k <- clusterboot(
    dist_matrix,
    B             = 500,
    distances     = TRUE,           # tells clusterboot data is a dissimilarity matrix
    clustermethod = disthclustCBI,  # distance-matrix-aware hierarchical interface
    method        = "ward.D2",
    k             = k,
    dissolution   = 0.5,
    recover       = 0.75
  )
  tibble(
    k         = k,
    mean_jacc = mean(boot_k$bootmean),
    min_jacc  = min(boot_k$bootmean)
  )
})

print(boot_results)

sil_boot <- left_join(sil_df, boot_results, by = "k")

print(sil_boot)

# ── Detailed bootstrap at chosen k ─────────────────────────────────────────────
boot_stab <- clusterboot(
  dist_matrix,
  B             = 1000,
  distances     = TRUE,
  clustermethod = disthclustCBI,
  method        = "ward.D2",
  k             = 2,
  dissolution   = 0.5,
  recover       = 0.75
)

print(boot_stab$bootmean)

# Plot bootstrap results and silhouette criteria together
sil_boot %>%
  pivot_longer(c(avg_sil, mean_jacc, n_negative, min_sil), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         avg_sil   = "Mean silhouette width",
                         mean_jacc = "Mean Jaccard (bootstrap)",
                         min_sil = "Minimum silhouette width",
                         n_negative = "Num. of negative sil. widths"
  )) %>%
  ggplot(aes(x = k, y = value, colour = metric)) +
  geom_line() +
  geom_point(size = 3) +
  facet_wrap(~metric, nrow = 4, scales = "free_y")+
  scale_x_continuous(breaks = 2:15) +
  labs(x = "Number of clusters (k)", y = "Evaluation metric",
       colour = NULL,
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

ggsave("output/figures/figS7_cluster_assess.jpeg", width = 8, height = 6.5)

# Select 2 based on peak in Jaccard
cuts   <- cutree(hc_ward_gower, k = 2)

# Calculate silhouette scores for these cuts; join with species data
sil_sp    <- silhouette(cuts, gower_dist) %>%
  as_tibble() %>%
  mutate(species = rownames(dist_mat_custom))

#Silhouette stats for results
range(sil_sp$sil_width)
filter(sil_sp, cluster == 1) %>% pull(sil_width) %>% length()
filter(sil_sp, cluster == 1) %>% pull(sil_width) %>% median()
filter(sil_sp, cluster == 2) %>% pull(sil_width) %>% length()
filter(sil_sp, cluster == 2) %>% pull(sil_width) %>% median()

sp_guilds <- sp_guilds %>%
  left_join(sil_sp) %>%
  mutate(name_clust = paste0(cluster, ". ", Common_Name))

# Dendrogram plot
ddata <- ggdendro::dendro_data(hc_ward_gower, type = "rectangle")

ddata_labels <- left_join(ddata$labels, sp_guilds, by = c("label" = "species")) %>%
  mutate(Habitat_Guild = str_replace(Habitat_Guild, "-", "-\n")) %>%
  mutate(Habitat_Guild = str_to_title(Habitat_Guild)) %>%   # do the casing BEFORE factoring
  mutate(Habitat_Guild = factor(Habitat_Guild, 
                                levels = c("Limnophilic", 
                                           "Limno-\nRheophilic",  
                                           "Pelagic Rheo-\nLimnophilic",                         
                                           "Rheo-\nLimnophilic", 
                                           "Rheophilic")))

n_leaves <- nrow(ddata$labels)

p1 <- 
  ggplot(segment(ddata)) +
  geom_tile(data = ddata_labels, 
            aes(x, y - .23, fill = Habitat_Guild), 
            width = 0.9, height = .46, alpha = .55) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = ddata_labels, aes(x, y - .23, label = name_clust), size = 3.5) +
  coord_flip(clip = "off", xlim = c(1, n_leaves), ylim = c(min(ddata_labels$y) - 0.35, NA))+
  labs(fill = "Habitat\nGuild") +
  scale_fill_viridis_d(option = "D") +
  theme_void(base_size = 14) +
  theme(plot.margin = margin(0, 0, 0, 10),
        legend.key.spacing.y = unit(0.2, "cm"))
  
p2 <- sil_sp %>% 
  mutate(species = factor(species, levels = ddata$labels$label)) %>%
  ggplot()+
  geom_bar(aes(x = sil_width, y = species), stat='identity') +
  scale_x_reverse(breaks = c(0, 0.3))+
  theme_classic(base_size = 11) +
  labs(y = NULL, x = "Sil. Width")+
  scale_y_discrete(position = "right", 
                   expand = c(0, 1.6)) +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) 

p2+ p1 + plot_layout(width = c(.2, 1))

ggsave(filename = "output/figures/fig4_dendro_hab_guilds.jpeg",  bg = "transparent", width = 8, height = 8, units = "in")

# Compare to other clustering versions
hc_complete <- hclust(gower_dist, method = "complete")
hc_average  <- hclust(gower_dist, method = "average")

# Compare dendrograms -- high correlation = stable clusters
cor_bakers <- dendextend::cor_bakers_gamma(
  as.dendrogram(hc_ward_gower),
  as.dendrogram(hc_average)
)
cat("Ward vs average dendrogram correlation:", round(cor_bakers, 3), "\n")

#Plot PDP data
pdp_all_long <- pdp_all %>%
  left_join(sil_sp, by = "species") %>%
  filter(variable %in% CONTINUOUS_VARS) %>%
  dplyr::select(-c(snag, substrt, period, vegd, neighbor, sil_width)) %>%
  pivot_longer(-c(variable, yhat, species, cluster), names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  left_join(importance_df_all) %>%
  mutate(var_name = recode(variable,
                           current   = "Current\nvelocity (m/s)",
                           depth   = "Depth (m)",
                           do   = "Dissolved\noxygen (mg/l)",
                           temp   = "Temperature (\u00B0C)"),
         cluster = paste("Cluster", cluster))

ggplot(pdp_all_long, aes(x = value, y = yhat, group = species, colour = factor(cluster))) +
  geom_line(alpha = 0.5, aes(lwd = importance)) +
  stat_summary(aes(group = cluster), fun = mean,
               geom = "line", linewidth = 1.5, colour = "black") +
  facet_grid(cluster ~ var_name, scales = "free", switch = "x") +
  labs(y = "Predicted probability\nof occurrence)") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
      strip.background = element_blank(),         
      strip.placement = "outside",
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank())

ggsave(filename = "output/figures/figS8_pdp_all.jpeg",  bg = "transparent", width = 11.5, height = 8, units = "in")

# Tile plot of signed importance
fixed_height <- 0.5

shape_with_importance %>%
  left_join(sil_sp) %>%
  left_join(sp_guilds) %>%
  filter(variable %in% c("current", "depth", "do", "temp")) %>%
  mutate(species = factor(species, levels = ddata$labels$label),
         Common_Name = fct_reorder(factor(Common_Name), as.numeric(species)),
                                   var_name = recode(variable,
                                                     current   = "Current\nvelocity",
                                                     depth   = "Depth",
                                                     do   = "Dissolved\noxygen",
                                                     temp   = "Temperature")) %>%
  ggplot(aes(x = var_name, 
             y = Common_Name,  # order by guild
             fill = cor_linear)) +
  geom_tile(aes(width = importance_scaled*2,      # Dynamically scales horizontally based on your column
            height = fixed_height), color = "black")+
  facet_grid(paste("Cluster",cluster) ~ ., scales = "free_y", space = "free_y") +
  scale_fill_distiller(palette = "BrBG", direction = 1) +
  labs(fill = "Linear\ncorrelation", x = "Variable", y = "Species")+
  theme_bw(base_size = 14)#+

ggsave(filename = "output/figures/fig5_cont_var_corr.jpeg",  bg = "transparent", width = 8, height = 5.5, units = "in")

# Categorical variables
# Single var
pdp_peak_cat <- pdp_all %>%
  filter(variable %in% CATEGORICAL_VARS) %>%
  dplyr::select(-c(current, depth, do, temp)) %>%
  pivot_longer(-c(variable, yhat, species), names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  left_join(importance_df_all) %>%
  group_by(species, variable) %>%
  slice_max(yhat) %>%
  dplyr::select(species, variable, value)

cat_peak <- shape_with_importance %>%
  left_join(sil_sp) %>%
  left_join(sp_guilds) %>%
  mutate(species = factor(species, levels = ddata$labels$label),
         Common_Name = fct_reorder(factor(Common_Name), as.numeric(species)),
         var_name = recode(variable,
                           period   = "Period",
                           snag   = "Snag",
                           substrt   = "Substrate",
                           vegd   = "Vegetation\nDensity")) %>% 
  left_join(pdp_peak_cat) %>%
  mutate(cat_level = factor(case_when(variable == "snag" & value == "0" ~ "Absent", 
                               variable == "snag" & value == "1" ~ "Present",
                               variable == "substrt" & value == "1" ~ "Very fine",
                               variable == "substrt" & value == "2" ~ "Fine",
                               variable == "substrt" & value == "3" ~ "Sandy",
                               variable == "substrt" & value == "4" ~ "Hard",
                               variable == "vegd" & value == "1" ~ "Sparse",
                               variable == "vegd" & value == "2" ~ "Dense", 
                               TRUE ~ value),
                            levels = c("1", "2", "3", "Absent", "Present",
                                       "Very fine", "Fine", "Sandy", "Hard",
                                       "Sparse", "Dense")))


col_pal <-  data.frame(var_name = levels(cat_peak$cat_level),
                       colors = c("#a6cee3", "#1f78b4", "#b2df8a",
                         "#1f78b4", "#b2df8a",
                         "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                         "#1f78b4", "#b2df8a"))

cat_plot_fun <- function(cat_var){
  cat_peak %>%
    filter(variable %in% c(cat_var)) %>%
    ggplot(aes(x = var_name, 
               y = Common_Name,  # order by guild
               fill = cat_level)) +
    geom_tile(aes(width = importance_scaled*2,      # Dynamically scales horizontally based on your column
                  height = fixed_height), color = "black")+
    facet_grid(paste("Cluster", cluster) ~ ., scales = "free_y", space = "free_y") +
    # scale_fill_manual(values = c("#1b9e77",
    #                               "#d95f02",
    #                               "#7570b3",
    #                               "gray40"))+
    scale_fill_manual(values = deframe(col_pal))+
    theme_bw(base_size = 14)+
    theme(legend.position = "bottom", 
          legend.direction = "vertical",
          legend.box.spacing = unit(0, "pt"),
          legend.margin = margin(t = 0, b = 0, unit = "pt"),
          #plot.margin = margin(b = 0, unit = "pt"),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 14))
}

x1 <- cat_plot_fun("period") + theme(strip.text = element_blank())
x2 <- cat_plot_fun("snag") + theme(strip.text = element_blank(),
                                   axis.title.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
x3 <- cat_plot_fun("substrt") + theme(strip.text = element_blank(),
                                      axis.title.y = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank())
x4 <- cat_plot_fun("vegd") + theme(axis.title.y = element_blank(),
                                   axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())


x1 + x2 + x3 + x4 +
  plot_layout(ncol = 4)

ggsave(filename = "output/figures/fig6_cat_var_corr.jpeg",  bg = "transparent", width = 8, height = 6, units = "in")


# OLD ##
# Other covariates
pdp_other_long <- pdp_all %>%
  left_join(sil_sp, by = "species") %>%
  filter(!variable %in% CONTINUOUS_VARS) %>%
  dplyr::select(-c(current, do, depth, temp, neighbor, sil_width)) %>%
  pivot_longer(-c(variable, yhat, species, cluster), names_to = "var2", values_to = "value") %>%
  filter(!is.na(value)) %>%
  left_join(importance_df_all)

pdp_other_long %>%
  ggplot(aes(x = value, y = yhat, colour = factor(cluster))) +
  geom_boxplot(outliers = FALSE) +
  geom_point(alpha = 0.2) +
  #geom_line(aes(group = species))+
  stat_summary(aes(group = cluster), fun = mean,
               geom = "line", linewidth = 1.5, colour = "black") +
  facet_grid(cluster ~ variable, scales = "free") +
  labs(x = "Predictor Value",
       y = "Predicted P(occurrence)",
       title = "PDP profiles by habitat guild") +
  theme_minimal() #+

# Variable Importance of other variables
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

rank_importance <- importance_df_all %>%
  group_by(species) %>%
  mutate(rank_imp = rank(-importance),
         importance_round = round(importance, 3)) %>%
  left_join(sil_sp) %>%
  group_by(cluster, variable) %>%
  mutate(outlier = ifelse(is_outlier(rank_imp), species, NA))


rank_importance %>%
  ggplot(aes(x = cluster, y = rank_imp, colour = factor(cluster))) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height = 0, width = 0.2,alpha = 0.2) +
  #geom_line(aes(group = species))+
  #stat_summary(aes(group = cluster), fun = mean,
  #            geom = "line", linewidth = 1.5, colour = "black") +
  facet_wrap(~ variable) +
  geom_text_repel(aes(label = outlier), na.rm = TRUE, hjust = -0.1) +
  scale_y_reverse()+
  labs(x = "Predictor Value",
       y = "Predicted P(occurrence)",
       title = "PDP profiles by habitat guild") +
  theme_minimal() #+

rank_importance %>%
  ggplot(aes(x = reorder(variable, rank_imp), y = rank_imp, colour = factor(cluster))) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height = 0, width = 0.2,alpha = 0.2) +
  #geom_line(aes(group = species))+
  #stat_summary(aes(group = cluster), fun = mean,
  #            geom = "line", linewidth = 1.5, colour = "black") +
  facet_wrap(~ cluster) +
  #geom_text_repel(aes(label = outlier), na.rm = TRUE, hjust = -0.1) +
  scale_y_reverse()+
  labs(x = "Predictor Value",
       y = "Predicted P(occurrence)",
       title = "PDP profiles by habitat guild") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45))

ggplot(rank_importance)+
  geom_point(aes(reorder(variable, rank_importance), rank_imp)) +
  facet_wrap(~species)+
  scale_y_reverse()

