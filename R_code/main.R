library(sf)
library(httr)
library(ggplot2)
library(dplyr)

# Correct URL for GeoJSON output
url <- "https://www.ogd.stadt-zuerich.ch/wfs/geoportal/Sammelstelle?service=WFS&version=2.0.0&request=GetFeature&outputFormat=GeoJSON&typeName=poi_sammelstelle_att"

response <- GET(url)

if (status_code(response) == 200) {
  # Read the GeoJSON content
  geojson_data <- content(response, "text")
  sammelstellen <- st_read(geojson_data, quiet = TRUE)
  
  # Check and set CRS if necessary
  if (is.na(st_crs(sammelstellen))) {
    st_crs(sammelstellen) <- 2056  # EPSG:2056 is CH1903+ / LV95, the new Swiss CRS
  }
  
  # Transform to WGS84 for better compatibility with ggplot2
  sammelstellen_wgs84 <- st_transform(sammelstellen, 4326)
  
  # Create the plot
  plot <- ggplot(data = sammelstellen_wgs84) +
    geom_sf() +
    theme_minimal() +
    labs(title = "Recycling Collection Points in Zurich")
  
  # Display the plot
  print(plot)
  
  print("Data successfully downloaded and map created.")
} else {
  print("Error downloading the data.")
}




