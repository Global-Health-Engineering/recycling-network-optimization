
# Load libraries
library(httr)
library(sf)
library(leaflet)
library(dplyr)
library(ows4R)
library(tidyverse)

wfs_url <- "https://www.ogd.stadt-zuerich.ch/wfs/geoportal/Sammelstelle"
wfs_client <- WFSClient$new(wfs_url, serviceVersion = "1.0.0")
wfs_client$getFeatureTypes(pretty = TRUE)

# Request zusammenstellen
layer = "poi_sammelstelle_view"

url <- parse_url(wfs_url)
url$query <- list(service = "WFS",
                  version = "1.0.0", 
                  typename = layer,
                  request = "GetFeature")
request <- build_url(url)

sammelstelle_geo <- sf::st_read(request)
sammelstelle_geo <- st_set_crs(sammelstelle_geo, 2056)


plot <- ggplot(data = sammelstelle_geo) +
  geom_sf() +
  theme_minimal() +
  labs(title = "RCP locations in Zurich") +
  theme(legend.position = "bottom")

print(plot)

sammelstelle_geo <- st_transform(sammelstelle_geo, 4326)

# Calculate the mean coordinates for the initial map center
map_center <- st_coordinates(st_centroid(st_union(sammelstelle_geo))) / 10  # Adjust if necessary

# Extract coordinates
coords <- st_coordinates(sammelstelle_geo)


# Create a leaflet map
leaflet_map <- leaflet(sammelstelle_geo) %>%
  addTiles() %>%
  setView(lng = mean(coords[,1]), 
          lat = mean(coords[,2]), 
          zoom = 13) %>%
  addCircleMarkers(
    lng = ~coords[,1],
    lat = ~coords[,2],
    radius = 5,
    color = "blue",
    stroke = FALSE,
    fillOpacity = 0.7,
    popup = ~name  # Replace 'name' with the appropriate column for popups
  )
print(leaflet_map)

