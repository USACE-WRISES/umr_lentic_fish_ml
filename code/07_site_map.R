library(tidyverse)  # For wrangling
library(sf)         # For everything spatial
library(basemaps)   # For basemap
library(maps)       # For US boundary data
library(tidyterra)  # To convert basemap to raster 
library(terra)      # To convert basemap to raster
library(ggspatial) # North arrows; scale bars
library(patchwork) # Adding insets to main map
library(here)
#library(ggsflabel)  # To get SF labels to repel; needs to be installed from github

# Import data ####

## States ####
states <- c("wisconsin", "minnesota", "iowa", "illinois", "michigan", 
            "missouri", "north dakota", "south dakota", "nebraska", 
            "kansas", "indiana", "oklahoma", "arkansas", "kentucky")

# Import state boundaries from the maps package; filtering to relevant states not strictly necessary
us <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE)) %>%
  filter(ID %in% states) %>%
  st_transform(crs = 3857)

us_simp <- st_simplify(us, dTolerance = 1e3)

# Import pool boundaries
pool_bounds <- st_read(here("data", "raw", "spatial", "UMRR_Pools", "UMRR_Pools.shp")) %>%
  st_transform(crs = 3857) %>%
  st_simplify(dTolerance = 500)

# Import survey points
surveys <- read_csv("data/cleaned/fish_all_cv.csv")

fish_sf <- filter(surveys, species == "LMBS") %>%
sf::st_as_sf(coords = c("utm_e", "utm_n"), crs = 26915) %>%
  mutate(Pool = paste("Pool", pool)) %>%
  st_transform(crs = 3857)

state_names <- data.frame(name = c("Wisconsin", "Minnesota", "Illinois", "Iowa"),
                       lat = c(43.4,43.9,42,43),
                       lon = c(-90.4,-92.1,-89.7,-92)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3857)

cities <- data.frame(
  name = c("Minneapolis", "Lake City", "La Crosse", "Bellevue", "Clinton"),
  lat  = c(44.9778, 44.4478-.03, 43.8014, 42.2586, 41.8280),
  lon = c(-93.2650, -92.2668, -91.2396 + .05, -90.4221 - .05, -90.2496))

cities_sf <- cities %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(3857)

cities_names <- cities %>%
  mutate(lat_shift = c(-.08, -.07, 0, 0, 0),
         lon_shift = c(-.24,-.03,.29,-.26,-.23),
         lat = lat + lat_shift,
         lon = lon + lon_shift) %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3857)
  
# Import river

# Create bounding box for plot area; convoluted code because you can't st_transform bbox objects
b_box_main <- st_bbox(c(xmin = -94, xmax = -89, ymax = 45, ymin = 41.5), crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  st_transform(crs = 3857) %>%
  st_bbox()

# Download basemap from basemap package; use get_maptypes() to look at other options
set_defaults(map_service = "esri", map_type = "world_terrain_base")
#set_defaults(map_service = "esri", map_type = "usa_topo_maps")

bm2 <- basemap_terra(ext = b_box_main) %>%
  as("SpatRaster")

pal_df <- data.frame(reach = unique(fish_sf$pool),
                     color = c("#1b9e77","#d95f02", "#7570b3"))

cols <- setNames(pal_df$color, pal_df$pool)

p_main <- ggplot() + 
  geom_spatraster_rgb(data = bm2) + # Map the baseline data using function from tidyterra package
  geom_sf(data= us_simp, fill = alpha("white", 0)) +
  geom_sf(data = pool_bounds, fill = "#CCEEF9")+
  geom_sf(data = fish_sf, aes(fill = Pool), size = 1.5, pch = 21, alpha = 0.8)+
  geom_sf(data = cities_sf) +
  geom_sf_text(data = state_names, aes(label = name), size = 6, fontface = "italic")+ #nudge_x = 2.3e4)+
  geom_sf_text(data = cities_names, aes(label = name), size = 3.5)+
  scale_fill_manual(values = cols)+
  coord_sf(xlim = c(b_box_main[1], b_box_main[3]), ylim = c(b_box_main[2], b_box_main[4]), 
           expand = FALSE)+
  labs(fill = "Fish Surveys\nLocations")+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
  theme_void(base_size = 14)+
  theme(plot.background = element_rect(fill = "transparent", color = NA), # Needed to make background transparent when saving
        panel.border = element_rect(color = "gray30", fill = NA, linewidth = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
  ) +
  annotation_scale(location = "tr", height = unit(.2, "cm"), text_cex = 1.1) +
  annotation_north_arrow(location = "br", which_north = "true", width = unit(1, "cm"), height = unit(1, "cm"))
p_main

# Inset plot ####
# Create bounding box for inset plot area; convoluted code because you can't st_transform bbox objects
b_box_ins <- st_bbox(c(xmin = -100, xmax = -85, ymax = 50, ymin = 35), crs = st_crs(4326)) %>%
  st_as_sfc() %>%
  st_transform(crs = 3857) %>%
  st_bbox()

lake_sup <- st_read("data/raw/spatial/hydro_p_LakeSuperior/hydro_p_LakeSuperior.shp")%>%
  st_transform(crs = 3857) %>%
  st_union() %>% # Need to dissolve two polygons (two rows) into one shape
  st_sf()

# Make inset plot ####
#Add inset below 
rect <- st_as_sfc(
  st_bbox(
    c(
      xmin = -99.3,
      ymin = 35.3,
      xmax = -90.1,
      ymax = 36.8
    ),
    crs = 4326)) %>%
  st_transform(3857)

p_ins <- 
  ggplot() + 
  geom_sf(data= us, fill = alpha("white", 0)) +
  geom_sf(data= lake_sup, fill = alpha("white", 0)) +
  geom_sf(data = st_as_sfc(b_box_main), fill = "gray", alpha = 0.5)+
  geom_sf(data = rect, fill = "white", color = NA) +
  coord_sf(xlim = c(b_box_ins[1], b_box_ins[3]), ylim = c(b_box_ins[2], b_box_ins[4]), 
           expand = FALSE)+
  theme_void(base_size = 14)+
  theme(plot.background = element_rect(fill = "white", #"transparent", # Needed to make background transparent when saving
                                       color = "black", linewidth = 1), 
        legend.title = element_blank()) + #Remove legend title
  annotation_scale(location = "bl", height = unit(.1, "cm"), text_cex = .8)

#l,b,r,t
p_main + inset_element(p_ins, -0.4, 0.02, .65, .35, align_to = 'full')

ggsave(filename = "output/figures/fig1_umr_map.jpeg",  bg = "transparent", width = 9, height = 6, units = "in")

# Supplemental plots
set.seed(1)
fish_sf <- fish_sf %>%
  mutate(Pool = paste("Pool", pool),
         patch_cluster = factor(
           patch_cluster,
           levels = sample(unique(patch_cluster))
         ))

# Color palettes
## Pool palette
pal_pool <- data.frame(pool = unique(fish_sf$cluster),
                       color = hue_pal()(length(unique(fish_sf$cluster))))

color_pool <- setNames(pal_pool$color, pal_pool$pool)

## Patch palette
pal_patch <- data.frame(
  patch = levels(fish_sf$patch_cluster),
  color = scales::hue_pal()(length(levels(fish_sf$patch_cluster)))
)

color_patch <- setNames(pal_patch$color, pal_patch$patch)

fish_plot <- function(clust_sub, clust_type, plot_title, color_pal) {
  fish_sf %>%
    dplyr::filter(cluster %in% clust_sub) %>%
    ggplot() +
    geom_sf(aes(color = as.factor(.data[[clust_type]]))) +
    guides(color = "none") +
    labs(title = plot_title)+
    scale_color_manual(values = color_pal) +
    theme_void()+
    theme(plot.title = element_text(hjust = 0.5), panel.background = element_rect(fill = NA, color = "black", linewidth = 1))
}

## Pool CV
p4u <- fish_plot("04u", "cluster", "Upper Pool 4", color_pool)
p4l <- fish_plot("04", "cluster", "Lower Pool 4", color_pool)
p8 <- fish_plot(c("08", "08u"), "cluster", "Pool 8", color_pool)
p13 <- fish_plot(c("13", "13u"), "cluster", "Pool 13", color_pool)

layout <- "
AACCDD
AACCDD
BBCCDD
BBCCDD
"

p4u + p4l +p8+p13 +
  plot_layout(design = layout)

ggsave(filename = "output/figures/figS1_map_pool_cv.jpeg",  bg = "transparent", width = 6, height = 4, units = "in")

## Patch CV
p4u <- fish_plot("04u", "patch_cluster", "Upper Pool 4", color_patch)
p4l <- fish_plot("04", "patch_cluster", "Lower Pool 4", color_patch)
p8 <- fish_plot(c("08", "08u"), "patch_cluster", "Pool 8", color_patch)
p13 <- fish_plot(c("13", "13u"), "patch_cluster", "Pool 13", color_patch)

p4u + p4l +p8+p13 +
  plot_layout(design = layout)

ggsave(filename = "output/figures/figS2_map_patch_cv.jpeg",  bg = "transparent", width = 6, height = 4, units = "in")
