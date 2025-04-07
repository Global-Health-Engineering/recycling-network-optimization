#!/usr/bin/env python3

import geopandas as gpd
import pandas as pd
from sklearn.cluster import DBSCAN
import folium
import branca.colormap as cm
from shapely.ops import unary_union
import sys
import os
import logging

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import calculate_duration, generate_isochrone, merge_isochrones_preserve_time

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

# Get routing engine from params
ROUTING_ENGINE = snakemake.params.get('routing_engine', 'valhalla')
logger.info(f"Using {ROUTING_ENGINE} routing engine")

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
    logger.info(f"Found {len(flats_outside)} flats outside isochrones")

    # Convert to centroids and set up the data for clustering
    X = pd.DataFrame({
        'x': flats_outside.geometry.x,
        'y': flats_outside.geometry.y,
        'population': flats_outside['est_pop']
    })

    # Apply DBSCAN clustering
    logger.info("Applying DBSCAN clustering")
    db = DBSCAN(eps=snakemake.params.eps, min_samples=snakemake.params.min_samples).fit(X[['x', 'y']])
    X['cluster'] = db.labels_

    # Remove noise points
    clusters = X[X['cluster'] != -1]
    logger.info(f"Found {len(clusters['cluster'].unique())} clusters (excluding noise)")

    # Calculate cluster centers weighted by population
    logger.info("Calculating population-weighted cluster centers")
    cluster_centers = clusters.groupby('cluster').apply(
        lambda df: pd.Series({
            'x': (df['x'] * df['population']).sum() / df['population'].sum(),
            'y': (df['y'] * df['population']).sum() / df['population'].sum()
        })
    ).reset_index()

    # Create GeoDataFrame for cluster centers
    cluster_centers_gdf = gpd.GeoDataFrame(
        cluster_centers,
        geometry=gpd.points_from_xy(cluster_centers['x'], cluster_centers['y']),
        crs="EPSG:4326"
    )

    # Filter potential sites with status "potential"
    potential_pot = potential_sites[potential_sites["status"] == "potential"].copy()
    if potential_pot.crs != "EPSG:4326":
        potential_pot = potential_pot.to_crs("EPSG:4326")

    # For each cluster centre, find the closest potential location
    logger.info("Finding closest potential sites to cluster centers")
    closest_locations = []
    for idx, centre in cluster_centers_gdf.iterrows():
        min_duration = float('inf')
        min_site = None
        
        for _, pot_site in potential_pot.iterrows():
            # Calculate duration using utility function
            duration = calculate_duration(
                (centre.geometry.x, centre.geometry.y),
                (pot_site.geometry.x, pot_site.geometry.y)
            )
            
            if duration and duration < min_duration:
                min_duration = duration
                min_site = pot_site
        
        if min_site is not None:
            closest_locations.append({
                'potential_ID': min_site['ID'],
                'geometry': min_site.geometry,
                'duration_minutes': min_duration
            })

    closest_locations_gdf = gpd.GeoDataFrame(closest_locations, geometry='geometry', crs="EPSG:4326")
    logger.info(f"Found {len(closest_locations_gdf)} closest potential sites")

    # Create a GeoDataFrame for existing RCPs
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

    # Save the cluster centers to a file 
    logger.info(f"Saving cluster centers to {snakemake.output.cluster_centers}")
    cluster_centers_gdf.to_file(snakemake.output.cluster_centers, driver='GPKG')
    
    # Generate map visualization
    logger.info("Generating map visualization")
    centroid = merged_isochrones.geometry.unary_union.centroid

    # Initialize the folium map centered on the centroid with specified tiles
    m = folium.Map(location=[centroid.y, centroid.x], zoom_start=12, control_scale=True)

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

    # Add cluster centers to map
    cluster_centers_layer = folium.FeatureGroup(name='Cluster Centers')
    for _, row in cluster_centers_gdf.iterrows():
        folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            popup=f'Cluster {row.cluster}',
            icon=folium.Icon(color='purple', icon='bullseye', prefix='fa')
        ).add_to(cluster_centers_layer)
    cluster_centers_layer.add_to(m)

    # Add markers for each selected new site
    new_rcp_layer = folium.FeatureGroup(name='New RCP Locations')
    for _, row in closest_locations_gdf.iterrows():
        folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            popup=f'ID: {row["potential_ID"]}',
            icon=folium.Icon(color='blue', icon='map-marker', prefix='fa')
        ).add_to(new_rcp_layer)
    new_rcp_layer.add_to(m)

    # Add underserved flats as red CircleMarkers
    underserved_layer = folium.FeatureGroup(name='Underserved Flats', show=False)
    for _, row in flats_outside.iterrows():
        folium.CircleMarker(
            location=[row.geometry.y, row.geometry.x],
            radius=3,
            fill=True,
            color='orange',
            fill_color='orange',
            fill_opacity=0.6,
            popup=f'Population: {row.est_pop:.2f}'
        ).add_to(underserved_layer)
    underserved_layer.add_to(m)

    # Add layer control
    folium.LayerControl().add_to(m)
    legend_html = '''
    <div style="
        position: fixed;
        bottom: 50px;
        left: 50px;
        width: 150px;
        height: 110px;
        border:2px solid grey;
        z-index:9999;
        font-size:14px;
        background-color: white;
        opacity: 0.8;
        padding: 10px;
    ">
        <b>Legend</b><br>
        <i class="fa fa-bullseye" style="color: purple"></i>&nbsp;Cluster Center<br>
        <i class="fa fa-map-marker" style="color: blue"></i>&nbsp;New RCP Location<br>
        <i class="fa fa-recycle" style="color: green"></i>&nbsp;Existing RCP<br>
        <i class="fa fa-circle" style="color: orange"></i>&nbsp;Underserved Flats
    </div>
    '''
    m.get_root().html.add_child(folium.Element(legend_html))

    # Save the map
    m.save(snakemake.output.html_map)
    logger.info(f"Map saved to {snakemake.output.html_map}")

if __name__ == "__main__":
    main()
