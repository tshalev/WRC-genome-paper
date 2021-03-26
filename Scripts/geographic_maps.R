library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggrepel)
library(ggmap)
library(maps)
library(ggsci)
library(ggspatial)
library(elevatr)
library(tidyverse)
source("~/UBC/GSAT/PhD/WRC/r_scripts/publication_theme.r")

setwd("~/UBC/GSAT/PhD/WRC/GS/wrc/snps/S_lines/")

land <- ne_download(scale = "large", returnclass = "sf", category = "physical", type = "land")
rivers <- ne_download(scale = "large", returnclass = "sf", category = "physical", type = "rivers_lake_centerlines")
lakes <- ne_download(scale = "large", returnclass = "sf", category = "physical", type = "lakes")
ocean <- ne_download(scale = "large", returnclass = "sf", category = "physical", type = "ocean")
canada <- ne_states(country = "canada", returnclass = "sf")
usa <- ne_states(country = "united states of america", returnclass = "sf")
naturalearth <- ne_download(scale = "large", category = "raster", type = "NE1_HR_LC_SR_W", continent = "north america")

sites <- read.table("parent_geog_locations_corrected.txt", header = T)
self_sites <- read.table("S_lines_parent_geog_locations.txt", header = T)
parent_popmap <- read.table("filtering_for_pop_gen/pop_gen_v3_snps_43929_snps/parents_strata.txt")

parent_popmap <- parent_popmap %>% 
  arrange(V1) %>% 
  rename(Parent = V1)

sites_popmap <- merge(sites, parent_popmap, by = "Parent")

col <- c(pal_nejm()(8), "#6A3D9A", "#A6CEE3")

ggplot() +
  theme_Publication() +
  geom_sf(data = land, fill ="antiquewhite") +
  geom_sf(data = canada, fill ="antiquewhite") +
  geom_sf(data = usa, fill ="antiquewhite") +
  # geom_sf(data = lakes, fill = "lightblue") +
  # geom_sf(data = rivers, colour = "lightblue") +
  # geom_sf(data = ocean, fill = "aliceblue") +
  # geom_text_repel(data = sites, aes(x = Longitude, y = Latitude), size = 2.5, label = sites$Parent) +
  geom_point(data = sites_popmap, aes(x = Longitude, y = Latitude, fill = V2), shape = 23, size = 2) +
  scale_fill_manual(name="Population",
                      breaks=c("Coastal_BC", "Haida_Gwaii", "Interior_BC", "US", "Vancouver_Island"),
                      labels=c("Coastal BC", "Haida Gwaii", "Interior BC", "Coastal NW US", "Vancouver Island"),
                      values = col) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-137, -117), ylim = c(39, 55), expand = FALSE) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed",
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  theme(legend.position = c(0.188, 0.3), legend.background = element_rect(colour = "black"))

# ggsave("all_parents_geog_locations_labels_colours.svg", width = 7.5, height = 10)

ggplot() +
  theme_Publication() +
  geom_sf(data = land, fill ="antiquewhite") +
  geom_sf(data = canada, fill ="antiquewhite") +
  geom_sf(data = usa, fill ="antiquewhite") +
  # geom_sf(data = lakes, fill = "lightblue") +
  # geom_sf(data = rivers, colour = "lightblue") +
  # geom_sf(data = ocean, fill = "aliceblue") +
  # geom_text_repel(data = sites, aes(x = Longitude, y = Latitude), size = 2.5, label = sites$Parent) +
  geom_point(data = sites_popmap, aes(x = Longitude, y = Latitude, fill = V3), shape = 23, size = 2) +
  scale_fill_manual(name="Subpopulation",
                    breaks=c("Central_Coast", "Discovery_Islands", "Haida_Gwaii", "Interior_BC", "Lower_Mainland", "Sunshine_Coast", "US",
                             "Vancouver_Island_Centre", "Vancouver_Island_North", "Vancouver_Island_South"),
                    labels=c("Central Coast", "Discovery Islands", "Haida Gwaii", "Interior BC", "Lower Mainland", "Sunshine Coast", "Coastal NW US",
                             "Vancouver Island (Centre)", "Vancouver Island (North)", "Vancouver Island (South)"),
                    values = col) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-137, -117), ylim = c(39, 55), expand = FALSE) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed",
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue")) +
  theme(legend.position = c(0.25, 0.35), legend.background = element_rect(colour = "black"))


# ggsave("all_parents_geog_locations_labels_subpop_colours.svg", width = 7.5, height = 10)

ggplot(data = world) +
  geom_sf(fill = "antiquewhite1") +
  # geom_text_repel(data = sites, aes(x = Longitude, y = Latitude), size = 2.5, label = sites$Parent) +
  geom_point(data = self_sites, aes(x = Long, y = Lat), shape = 23, size = 2, fill = "darkgreen") +
  scale_fill_manual(values = col_VI) +
  coord_sf(xlim = c(-129, -123), ylim = c(48.25, 52.25), expand = FALSE) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.5), panel.background = element_rect(fill = "aliceblue"))

# ggsave("self_parents_geog_locations.svg", width = 8, height = 9)


