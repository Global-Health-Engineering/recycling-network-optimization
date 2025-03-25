#!/usr/bin/env python
import numpy as np
import geopandas as gpd
import fiona
import pandas as pd
import rasterio
from rasterio.mask import mask
from shapely.geometry import MultiPolygon
from rasterstats import zonal_stats
import folium
from shapely.wkt import loads, dumps
from shapely.ops import transform
import logging
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stdout)]
)

def suitability_analysis(
    buffer_dist_residential: float,
    building_footprints: gpd.GeoDataFrame,
    buffer_dist_vbz: float,
    area_threshold: float,
    buffer_trees: float,
    max_slope: float,
    parking_lots: gpd.GeoDataFrame,
    slope_raster,
    trees: gpd.GeoDataFrame,
    vbz_lines: gpd.GeoDataFrame,
    vbz_points: gpd.GeoDataFrame
):
    # Step 1: Filter parking types
    parking_filtered = parking_lots[
        ~parking_lots['parking'].isin(['underground', 'multi-storey'])
    ].copy()

    # Step 2: Process VBZ features
    vbz_line_buffers = vbz_lines.buffer(buffer_dist_vbz)
    vbz_point_buffers = vbz_points.buffer(buffer_dist_vbz)

    # Step 3: Buffer trees
    tree_buffers = trees.buffer(buffer_trees)

    # step 3.1: Buffer residential buildings
    buffer_buildings = building_footprints.geometry.apply(lambda x: x.buffer(buffer_dist_residential))

    # Step 4: Combine all buffers
    all_buffers = gpd.GeoSeries(
        list(vbz_line_buffers) + list(vbz_point_buffers) +
        list(tree_buffers) +
        list(buffer_buildings),
        crs=parking_lots.crs
    ).unary_union

    # Step 5: Zonal statistics for slope
    slope_stats = zonal_stats(
        parking_filtered,
        slope_raster.name,  # Use the raster filepath from the raster object
        stats=['mean', 'max'],
        nodata=slope_raster.nodata
    )

    # Add slope stats to parking lots
    parking_filtered['slope_mean'] = [s['mean'] for s in slope_stats]
    parking_filtered['slope_max'] = [s['max'] for s in slope_stats]

    # Step 6: Filter by slope
    slope_filtered = parking_filtered[parking_filtered['slope_mean'] <= max_slope]

    # Step 7: Spatial difference with buffers
    final_areas = slope_filtered.geometry.difference(all_buffers)

    # Step 8: Calculate areas and filter
    result = gpd.GeoDataFrame(geometry=final_areas, crs=parking_lots.crs)
    result['area'] = result.geometry.area
    suitable = result[result['area'] >= area_threshold]

    return suitable

def main():
    logging.info("Starting suitability analysis")
    
    # Load the raster
    logging.info("Loading slope raster")
    slope_raster = rasterio.open(snakemake.input.slope_raster)
    
    # Load datasets
    logging.info("Loading input datasets")
    tree_dataset = gpd.read_file(snakemake.input.trees).to_crs(slope_raster.crs)
    parking_lots = gpd.read_file(snakemake.input.parking_lots).to_crs(slope_raster.crs)
    rcps = gpd.read_file(snakemake.input.rcps).to_crs(slope_raster.crs)
    building_footprints = gpd.read_file(snakemake.input.buildings).to_crs(slope_raster.crs)
    
    # Filter buildings and turn them into 2d polygons
    logging.info("Processing building footprints")
    building_footprints = building_footprints[building_footprints['art_txt'].str.contains('wohn', case=False, na=False)]
    building_footprints['geometry'] = building_footprints.geometry.apply(
        lambda geom: loads(dumps(geom, output_dimension=2))
    )
    
    # List all layers in the vbz geopackage
    logging.info("Processing VBZ data")
    layers = fiona.listlayers(snakemake.input.vbz)
    line_layers = []
    point_layers = []
    
    # Loop through layers and separate those that contain only line or point geometries
    for layer in layers:
        gdf = gpd.read_file(snakemake.input.vbz, layer=layer)
        if gdf.geom_type.str.startswith("Line").all():
            line_layers.append(gdf)
        elif gdf.geom_type.str.startswith("Point").all():
            point_layers.append(gdf)
    
    # Merge all line layers into one GeoDataFrame
    vbz_lines = gpd.GeoDataFrame(
        pd.concat(line_layers, ignore_index=True), crs=line_layers[0].crs
    )
    
    # Merge all point layers into one GeoDataFrame
    vbz_points = gpd.GeoDataFrame(
        pd.concat(point_layers, ignore_index=True), crs=point_layers[0].crs
    )

    # Get parameters from snakemake
    buffer_dist_vbz = snakemake.params.get("buffer_dist_vbz", 2)
    buffer_trees = snakemake.params.get("buffer_trees", 2)
    max_slope = snakemake.params.get("max_slope", 5)
    area_threshold = snakemake.params.get("area_threshold", 16)
    buffer_buildings = snakemake.params.get("buffer_buildings", 14)
    
    # Apply suitability analysis
    logging.info("Running suitability analysis")
    suitable_areas = suitability_analysis(
        building_footprints=building_footprints,
        buffer_dist_residential=buffer_buildings,
        area_threshold=area_threshold,
        buffer_dist_vbz=buffer_dist_vbz,
        buffer_trees=buffer_trees,
        max_slope=max_slope,
        parking_lots=parking_lots,
        slope_raster=slope_raster,
        trees=tree_dataset,
        vbz_lines=vbz_lines,
        vbz_points=vbz_points
    )
    
    # Process existing RCP sites
    logging.info("Processing existing and potential sites")
    existing_sites = rcps.copy()
    existing_sites['ID'] = ['e_{}'.format(i+1) for i in existing_sites.index]
    existing_sites['status'] = 'open'
    existing_sites = existing_sites[['geometry', 'ID', 'status']]
    
    # Create a buffer around existing open sites
    buffer_union = existing_sites.geometry.buffer(125).unary_union
    
    # Filter out potential sites within the buffer
    filtered_potential_sites = suitable_areas[~suitable_areas.geometry.intersects(buffer_union)].copy()
    filtered_potential_sites['ID'] = ['p_{}'.format(i+1) for i in filtered_potential_sites.index]
    filtered_potential_sites['status'] = 'potential'
    filtered_potential_sites = filtered_potential_sites[['geometry', 'ID', 'status']]
    
    # Merge existing and potential sites
    merged_sites = pd.concat([existing_sites, filtered_potential_sites], ignore_index=True)
    
    # Make sure all geometries are points
    if not merged_sites['geometry'].geom_type.str.startswith('Point').all():
        merged_sites = gpd.GeoDataFrame(merged_sites, geometry=merged_sites.geometry.centroid, crs=merged_sites.crs)
    merged_sites['geometry'] = merged_sites.geometry.centroid
    
    # Save the results
    logging.info(f"Saving {len(merged_sites)} sites ({len(existing_sites)} existing, {len(filtered_potential_sites)} potential)")
    merged_sites.to_file(snakemake.output.potential_sites, driver='GPKG')
    
    # Create a map showing all sites
    logging.info("Creating map of suitable sites")
    
    # Ensure the data is in WGS84 for mapping
    merged_sites_wgs84 = merged_sites.to_crs('EPSG:4326')
    suitable_areas_wgs84 = suitable_areas.to_crs('EPSG:4326')
    
    # Create a base map centered on Zurich
    m = folium.Map(location=[47.3769, 8.5417], zoom_start=12)
    
    # Add all suitable areas as polygons (light grey)
    folium.GeoJson(
        suitable_areas_wgs84,
        name='All Suitable Areas',
        style_function=lambda x: {
            'fillColor': '#eeeeee',
            'color': '#aaaaaa',
            'weight': 1,
            'fillOpacity': 0.6
        },
        tooltip=folium.GeoJsonTooltip(fields=['area'], aliases=['Area (m²)'])
    ).add_to(m)
    
    # Add existing RCPs (blue)
    existing_sites_layer = folium.FeatureGroup(name='Existing RCPs')
    for idx, row in merged_sites_wgs84[merged_sites_wgs84['status'] == 'open'].iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=5,
            color='#1f77b4',
            fill=True,
            fill_color='#1f77b4',
            fill_opacity=0.8,
            popup=f"ID: {row['ID']}",
            tooltip=f"Existing RCP: {row['ID']}"
        ).add_to(existing_sites_layer)
    existing_sites_layer.add_to(m)
    
    # Add potential new sites (orange)
    potential_sites_layer = folium.FeatureGroup(name='Potential Sites')
    for idx, row in merged_sites_wgs84[merged_sites_wgs84['status'] == 'potential'].iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=3,
            color='#ff7f0e',
            fill=True,
            fill_color='#ff7f0e',
            fill_opacity=0.8,
            popup=f"ID: {row['ID']}",
            tooltip=f"Potential Site: {row['ID']}"
        ).add_to(potential_sites_layer)
    potential_sites_layer.add_to(m)
    
    # Add layer control
    folium.LayerControl().add_to(m)
    
    # Save the map
    m.save(snakemake.output.map)
    logging.info("Suitability analysis completed successfully")

if __name__ == "__main__":
    main()
