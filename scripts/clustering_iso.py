#!/usr/bin/env python3

import geopandas as gpd
import pandas as pd
from sklearn.cluster import DBSCAN
import folium
import branca.colormap as cm
from shapely.ops import unary_union
import sys
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(snakemake.log[0]),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def merge_isochrones_preserve_time(isochrones_gdf):
    """
    Merge isochrones preserving lower time values.

    Parameters:
    - isochrones_gdf: GeoDataFrame with isochrones and 'time' attribute.

    Returns:
    - GeoDataFrame with merged isochrones.
    """
    # Ensure CRS is EPSG:4326
    if isochrones_gdf.crs != "EPSG:4326":
        isochrones_gdf = isochrones_gdf.to_crs(epsg=4326)

    # Sort isochrones by 'time' ascending
    isochrones_sorted = isochrones_gdf.sort_values(by='time')

    merged_isochrones = gpd.GeoDataFrame(columns=isochrones_sorted.columns, crs="EPSG:4326")

    # Initialize an empty geometry for subtraction
    accumulated_geom = None

    for _, row in isochrones_sorted.iterrows():
        current_geom = row.geometry
        current_time = row['time']

        if accumulated_geom:
            remaining_geom = current_geom.difference(accumulated_geom)
        else:
            remaining_geom = current_geom

        if not remaining_geom.is_empty:
            new_row = row.copy()
            new_row.geometry = remaining_geom
            # Ensure the new_row GeoDataFrame has the correct CRS
            new_row = gpd.GeoDataFrame([new_row], crs="EPSG:4326")
            merged_isochrones = pd.concat([merged_isochrones, new_row], ignore_index=True)
            # Update accumulated geometry
            if accumulated_geom:
                accumulated_geom = unary_union([accumulated_geom, remaining_geom])
            else:
                accumulated_geom = remaining_geom
    return merged_isochrones

def main():
    # Load inputs from Snakemake
    logger.info("Loading input files")
    existing_isochrones = gpd.read_file(snakemake.input.isochrones)
    flats_pop = gpd.read_file(snakemake.input.flats)
    rcps = gpd.read_file(snakemake.input.rcps)
    potential_sites = gpd.read_file(snakemake.input.potential_sites)
    
    # Ensure consistent CRS
    rcps = rcps.to_crs('EPSG:4326')

    # Merge isochrones preserving time
    logger.info("Merging isochrones")
    merged_isochrones = merge_isochrones_preserve_time(existing_isochrones)
    merged_isochrones.to_file(snakemake.output.merged_isochrones, driver='GPKG')

    # Reproject flats_pop to match merged_isochrones CRS
    logger.info("Identifying flats outside isochrones")
    flats_pop_4326 = flats_pop.to_crs(merged_isochrones.crs)

    # Merge all isochrones into a single geometry
    iso_union = merged_isochrones.unary_union

    # Identify flats outside any isochrones
    flats_outside = flats_pop_4326[~flats_pop_4326.geometry.within(iso_union)]
    
    # Apply DBSCAN clustering
    logger.info("Applying DBSCAN clustering")
    # Convert to centroids and set up the data for clustering
    X = pd.DataFrame({
        'x': flats_outside.geometry.x,
        'y': flats_outside.geometry.y,
        'population': flats_outside['est_pop']
    })

    # Apply DBSCAN clustering with parameters from input
    db = DBSCAN(eps=snakemake.params.eps, min_samples=snakemake.params.min_samples).fit(X[['x', 'y']])
    X['cluster'] = db.labels_

    # Remove noise points
    clusters = X[X['cluster'] != -1]

    # Calculate cluster centers weighted by population
    logger.info("Calculating weighted cluster centers")
    cluster_centers = clusters.groupby('cluster').apply(
        lambda df: pd.Series({
            'x': (df['x'] * df['population']).sum() / df['population'].sum(),
            'y': (df['y'] * df['population']).sum() / df['population'].sum()
        })
    ).reset_index()

    # Create GeoDataFrame for new collection points
    cluster_centers_gdf = gpd.GeoDataFrame(
        cluster_centers,
        geometry=gpd.points_from_xy(cluster_centers['x'], cluster_centers['y']),
        crs="EPSG:4326"
    )
    
    # Find closest potential sites
    logger.info("Finding closest potential sites to cluster centers")
    potential_pot = potential_sites[potential_sites["status"] == "potential"].copy()

    # Reproject potential sites to EPSG:4326 if needed
    if potential_pot.crs != "EPSG:4326":
        potential_pot = potential_pot.to_crs("EPSG:4326")

    # For each cluster centre, find the closest potential location
    closest_locations = []
    for idx, centre in cluster_centers_gdf.iterrows():
        # Compute distances from this centre to all potential sites
        potential_pot['dist'] = potential_pot.geometry.distance(centre.geometry)
        # Get the potential site with the minimum distance
        min_idx = potential_pot['dist'].idxmin()
        min_loc = potential_pot.loc[min_idx]
        closest_locations.append({
            'potential_ID': min_loc['ID'],
            'geometry': min_loc.geometry
        })

    closest_locations_gdf = gpd.GeoDataFrame(closest_locations, geometry='geometry', crs="EPSG:4326")
    
    # Create a GeoDataFrame for existing RCPs
    logger.info("Creating output GeoDataFrame")
    existing = rcps.copy()
    existing['id'] = ['existing_' + str(i + 1) for i in range(len(existing))]

    # Create a GeoDataFrame for potential RCPs
    potentials = closest_locations_gdf.copy()
    potentials['id'] = ['pot_' + str(i + 1) for i in range(len(potentials))]

    # Combine the two groups and then select only the required columns
    rcp_summary = pd.concat([
        existing[['geometry', 'id']], 
        potentials[['geometry', 'id']]
    ]).reset_index(drop=True)

    # Export to file
    logger.info(f"Saving results to {snakemake.output.clustered_sites}")
    rcp_summary.to_file(snakemake.output.clustered_sites, driver='GPKG')
    
    # Generate map visualization
    logger.info("Generating map visualization")
    centroid = merged_isochrones.geometry.unary_union.centroid

    # Initialize the folium map centered on the centroid with specified tiles
    m = folium.Map(location=[centroid.y, centroid.x], zoom_start=12)

    # Convert 'time' column to numeric and convert seconds to minutes
    merged_isochrones['time'] = pd.to_numeric(merged_isochrones['time']) / 60

    # Define a viridis colormap based on time (in minutes)
    colormap = cm.linear.viridis.scale(
        merged_isochrones['time'].min(),
        merged_isochrones['time'].max()
    )
    colormap.caption = 'Walking Time (minutes)'
    colormap.add_to(m)

    folium.GeoJson(
        merged_isochrones,
        name='Merged Isochrones',
        style_function=lambda feature: {
            'fillColor': colormap(float(feature['properties']['time'])),
            'color': colormap(float(feature['properties']['time'])),
            'weight': 1,
            'fillOpacity': 0.5,
        },
        show=False
    ).add_to(m)

    # Add RCP dataset to the map with green markers
    rcp_layer = folium.FeatureGroup(name='RCP Locations')
    for _, row in rcps.iterrows():
        folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            popup=row['adresse'] if 'adresse' in row else '',
            icon=folium.Icon(color='green', icon='recycle', prefix='fa')
        ).add_to(rcp_layer)
    rcp_layer.add_to(m)

    # Add flats_outside as red CircleMarkers within a FeatureGroup
    flats_outside_layer = folium.FeatureGroup(name='Flats Outside', show=False)
    for _, row in flats_outside.iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=3,
            fill=True,
            color='red',
            fill_color='red',
            fill_opacity=0.6,
            popup=f'Population: {row.est_pop:.2f}'
        ).add_to(flats_outside_layer)
    flats_outside_layer.add_to(m)

    # Add new collection points to the map with blue + sign markers
    new_rcp_layer = folium.FeatureGroup(name='New RCP Locations', show=False)
    for _, point in potentials.iterrows():
        folium.Marker(
            location=[point.geometry.y, point.geometry.x],
            icon=folium.Icon(color='blue', icon='plus', prefix='fa')
        ).add_to(new_rcp_layer)
    new_rcp_layer.add_to(m)

    folium.LayerControl().add_to(m)

    # Save the map to an HTML file
    m.save(snakemake.output.html_map)
    
    logger.info("Clustering completed successfully")

if __name__ == "__main__":
    main()
