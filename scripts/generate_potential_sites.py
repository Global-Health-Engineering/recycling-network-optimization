#!/usr/bin/env python
import numpy as np
import geopandas as gpd
import fiona
import pandas as pd
import rasterio
from rasterstats import zonal_stats
import folium
from shapely.wkt import loads, dumps
from pathlib import Path

def suitability_analysis(
    buffer_dist_residential: float,
    building_footprints: gpd.GeoDataFrame,
    buffer_dist_vbz: float,
    area_threshold: float,
    buffer_trees: float,
    max_slope: float,
    parking_lots: gpd.GeoDataFrame,
    slope_raster: rasterio.DatasetReader,
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
        slope_raster.name,
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

# Main function
if __name__ == "__main__":
    # Get paths from Snakemake
    slope_raster_path = str(snakemake.input.slope_raster)
    
    # Load the raster
    slope_raster = rasterio.open(slope_raster_path)

    # Load the data and ensure CRS matches the raster
    tree_dataset = gpd.read_file(str(snakemake.input.trees)).to_crs(slope_raster.crs)
    parking_lots = gpd.read_file(str(snakemake.input.parking)).to_crs(slope_raster.crs)
    rcps = gpd.read_file(str(snakemake.input.rcps)).to_crs(slope_raster.crs)
    building_footprints = gpd.read_file(str(snakemake.input.buildings)).to_crs(slope_raster.crs)

    # Filter buildings and turn them into 2D polygons
    building_footprints = building_footprints[building_footprints['art_txt'].str.contains('wohn', case=False, na=False)]
    building_footprints['geometry'] = building_footprints.geometry.apply(
        lambda geom: loads(dumps(geom, output_dimension=2))
    )

    # List all layers in the vbz geopackage
    layers = fiona.listlayers(str(snakemake.input.vbz))
    line_layers = []
    point_layers = []

    # Loop through layers and separate those that contain only line or point geometries
    for layer in layers:
        gdf = gpd.read_file(str(snakemake.input.vbz), layer=layer)
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

    # Define parameters
    buffer_dist_vbz = snakemake.params.get("buffer_dist_vbz", 2)  # meters buffer around VBZ infrastructure
    buffer_trees = snakemake.params.get("buffer_trees", 2)  # meters buffer around trees
    max_slope = snakemake.params.get("max_slope", 5)  # Maximum slope in degrees
    area_threshold = snakemake.params.get("area_threshold", 16)  # Minimum area in square meters
    buffer_buildings = snakemake.params.get("buffer_buildings", 14)  # meters buffer around buildings

    # Apply suitability analysis
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

    # Assign unique IDs to existing sites and add status
    existing_sites = rcps.copy()
    existing_sites['ID'] = ['e_{}'.format(i+1) for i in existing_sites.index]
    existing_sites['status'] = 'open'
    existing_sites = existing_sites[['geometry', 'ID', 'status']]

    # Create a buffer of 125 meters around existing open sites
    buffer = existing_sites.geometry.buffer(125)

    # Merge all buffers into a single geometry
    buffer_union = buffer.unary_union

    # Filter out potential sites within the buffer
    filtered_potential_sites = suitable_areas[~suitable_areas.geometry.intersects(buffer_union)].copy()

    # Add unique IDs and status to filtered potential sites
    filtered_potential_sites['ID'] = ['p_{}'.format(i+1) for i in filtered_potential_sites.index]
    filtered_potential_sites['status'] = 'potential'
    filtered_potential_sites = filtered_potential_sites[['geometry', 'ID', 'status']]

    # Merge the existing and filtered potential sites
    merged_sites = pd.concat([existing_sites, filtered_potential_sites], ignore_index=True)

    # Make sure all data points are points
    if not merged_sites['geometry'].geom_type.str.startswith('Point').all():
        merged_sites = gpd.GeoDataFrame(merged_sites, geometry=merged_sites.geometry.centroid, crs=merged_sites.crs)

    merged_sites['geometry'] = merged_sites.geometry.centroid

    # Export the merged sites to a GeoPackage
    merged_sites.to_file(str(snakemake.output.sites), driver='GPKG')

    # Create an HTML map
    m = folium.Map(location=[47.3769, 8.5417], zoom_start=12)
    
    # Add suitable areas to the map
    folium.GeoJson(
        merged_sites,
        name='Suitable Areas',
        popup=folium.GeoJsonPopup(fields=['ID', 'status']),
    ).add_to(m)
    
    # Save the map to an HTML file
    m.save(str(snakemake.output.map))
