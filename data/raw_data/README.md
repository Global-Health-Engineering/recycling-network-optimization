# Raw Data

This directory contains all the raw, unprocessed data needed for the RCP optimization analysis. The actual data files are not tracked in the git repository and will be downloaded automatically by the Snakemake workflow.

## Structure

- `geodata_stadt_Zuerich/` - Geographic data from the City of Zurich
  - `3d_buildings/` - 3D building models for Zurich
    - `data/` - Contains geometric data files
  - `building_stats/` - Building statistics information
    - `data/` - Contains ssz.gwr_stzh_wohnungen.shp
  - `municipal_boundaries/` - Administrative boundaries
    - `data/` - Contains Gemeindegrenzen_-OGD.gpkg
  - `population/` - Population distribution data
    - `data/` - Contains BEVOELKERUNG_HA_F.shp
  - `recycling_sammelstellen/` - Existing recycling collection points
    - `data/` - Contains stzh.poi_sammelstelle_view.shp
  - `trees/` - Tree location data
    - `data/` - Contains tree data files
  - `vbz/` - Public transport (VBZ) data
    - `data/` - Contains VBZ data files
- `osm_data/` - OpenStreetMap extracts
  - Contains parking_lots_zurich.gpkg
- `valhalla/` - Routing data for the Valhalla routing engine (if used)

## Data Sources

All datasets are open source. For detailed metadata, please refer to the original sources.

- **Building registry**: Building statistics for Zurich  
  https://www.stadt-zuerich.ch/geodaten/download/Gebaeude__und_Wohnungsregister_der_Stadt_Zuerich__GWZ__gemaess_GWR_Datenmodell

- **Population density raster**: Population distribution in hectare cells  
  https://www.stadt-zuerich.ch/geodaten/download/63

- **Existing RCPs**: Current recycling collection points  
  https://www.stadt-zuerich.ch/geodaten/download/Sammelstelle

- **Trees**: Tree location data  
  https://www.stadt-zuerich.ch/geodaten/download/Baumkataster

- **3D building model**: 3D buildings of Zurich  
  https://www.stadt-zuerich.ch/geodaten/download/Bauten___Kombinierte_Darstellung_heute

- **Public transport infrastructure**: VBZ data  
  https://www.stadt-zuerich.ch/geodaten/download/VBZ_Infrastruktur_OGD

- **Car parks**: Parking locations from OpenStreetMap  
  https://www.openstreetmap.org/ (downloaded using QuickOSM QGIS plugin: https://docs.3liz.org/QuickOSM/)

- **Elevation model of Zurich**: Terrain elevation data  
  https://www.swisstopo.admin.ch/en/height-model-swissalti3d#swissALTI3D--Download

- **Elevation model for routing**: NASA elevation data  
  https://lpdaac.usgs.gov/products/nasadem_hgtv001

- **Slope data**: The slope_zurich file needs to be generated in GIS software  
  > Note: This file needs to be manually generated using elevation models, which are readily available for Switzerland. Use QGIS or similar GIS software to derive slope from the elevation models listed above.

- **Road network**: OpenStreetMap road data for Switzerland  
  https://download.geofabrik.de/europe/switzerland.html
