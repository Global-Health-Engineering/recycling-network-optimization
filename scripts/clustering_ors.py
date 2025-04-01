import geopandas as gpd
import pandas as pd
import folium
import branca.colormap as cm
from sklearn.cluster import DBSCAN
import scripts.util as util
import sys
import os
import logging
import openrouteservice as ors
from openrouteservice.exceptions import ApiError

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import calculate_duration, initialize_ball_tree, find_nearest_rcp_duration, get_ors_client

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


# 1. Data Loading
def load_data():
    flats = gpd.read_file(snakemake.input['flats']).to_crs("EPSG:4326")
    rcps = gpd.read_file(snakemake.input['rcps']).to_crs("EPSG:4326")
    potential_sites = gpd.read_file(snakemake.input['potential_sites']).to_crs("EPSG:4326")
    return flats, rcps, potential_sites

# 2. Aggregating Flats by Building
def aggregate_flats(flats):
    return flats.groupby('egid').agg({'est_pop': 'sum', 'geometry': 'first'}).reset_index()

# 3. Compute Nearest RCP Durations
def compute_nearest_durations(buildings, rcps, route_url="http://localhost:8002/route"):

    # Prepare RCPs for nearest neighbor search
    tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps, 'poi_id')
    durations = buildings.copy()

    # Process each building to find its nearest RCP and the associated duration
    results = []
    for _, building in buildings.iterrows():
        # This utility call will use the correct routing engine based on the provided route_url
        nearest_id, duration = util.find_nearest_rcp_duration(
          building.geometry, tree, rcp_coords, rcp_ids, route_url
        )
        results.append((nearest_id, duration))

    durations['nearest_rcp_id'] = [r[0] for r in results]
    durations['duration'] = [r[1] for r in results]
    return gpd.GeoDataFrame(durations, geometry='geometry', crs="EPSG:4326")

# 4. Clustering Flats with High Population and Long Durations
def cluster_flats(flats_duration):
    iso_threshold = snakemake.params['iso_threshold']
    high_pop_unserved = flats_duration[
        (flats_duration['est_pop'] > 0) & 
        (flats_duration['duration'] >= iso_threshold)
    ].copy()
    high_pop_unserved = gpd.GeoDataFrame(high_pop_unserved, geometry='geometry', crs="EPSG:4326")
    
    coords = high_pop_unserved.geometry
    X = pd.DataFrame({
        'x': coords.x,
        'y': coords.y,
        'population': high_pop_unserved['est_pop']
    })
    db = DBSCAN(eps=snakemake.params.eps, min_samples=snakemake.params.min_samples).fit(X[['x', 'y']])
    X['cluster'] = db.labels_
    clusters = X[X['cluster'] != -1]

    # Calculate weighted cluster centers
    cluster_centers = clusters.groupby('cluster').apply(
        lambda df: pd.Series({
            'x': (df['x'] * df['population']).sum() / df['population'].sum(),
            'y': (df['y'] * df['population']).sum() / df['population'].sum()
        })
    ).reset_index()
    
    return gpd.GeoDataFrame(
        cluster_centers,
        geometry=gpd.points_from_xy(cluster_centers['x'], cluster_centers['y']),
        crs="EPSG:4326"
    ), high_pop_unserved

# 5. Find Closest Potential Locations for Each Cluster Centre
def find_closest_potential(cluster_centers, potential_sites, route_url="http://localhost:8002/route"):
    potential_pot = potential_sites[potential_sites["status"] == "potential"].copy()
    closest_locations = []
    
    for _, centre in cluster_centers.iterrows():
        centre_coords = (centre.geometry.x, centre.geometry.y)
        
        # Calculate durations for all potential sites
        durations = []
        for idx, site in potential_pot.iterrows():
            site_coords = (site.geometry.x, site.geometry.y)
            duration = util.calculate_duration(centre_coords, site_coords, route_url=route_url)
            durations.append(duration)
        
        potential_pot['duration'] = durations
        
        # If no duration could be calculated, skip this centre
        if potential_pot['duration'].isnull().all():
            continue
        
        # Select the potential site with the minimal duration
        min_idx = potential_pot['duration'].idxmin()
        min_loc = potential_pot.loc[min_idx]
        closest_locations.append({
            'potential_ID': min_loc['ID'],
            'geometry': min_loc.geometry,
            'duration': round(min_loc['duration'], 2) if min_loc['duration'] is not None else None
        })
    
    closest_locations_gdf = gpd.GeoDataFrame(closest_locations, geometry='geometry', crs="EPSG:4326")
    closest_locations_gdf['id'] = ['pot_' + str(i) for i in range(len(closest_locations_gdf))]
    return closest_locations_gdf

# 6. Merge Existing Recycling Points with New Potentials
def merge_locations(rcps, closest_locations_gdf):
    rcps['id'] = ['existing_' + str(i) for i in range(len(rcps))]
    merged_locations = pd.concat([closest_locations_gdf[['geometry', 'id']], rcps[['geometry', 'id']]])
    return merged_locations.reset_index(drop=True)

# 7. Build and Save the Folium Map
def build_map(rcps, closest_locations_gdf, flats_duration):
    colormap = cm.linear.viridis.scale(flats_duration['duration'].min(), 15)
    colormap.caption = 'Duration to Closest RCP (minutes)'
    outlier_color = "#ff0000"

    m = folium.Map(location=[47.3769, 8.5417], zoom_start=12, control_scale=True)

    existing_points = folium.FeatureGroup(name='Existing Recycling Points', show=False)
    new_recycling_points = folium.FeatureGroup(name='New Recycling Points', show=False)
    flats_layer = folium.FeatureGroup(name='Flats', show=False)

    # Add Existing RCP Markers
    for idx, row in rcps.iterrows():
        folium.Marker(
            location=[row['geometry'].y, row['geometry'].x],
            icon=folium.Icon(color='green', icon='recycle', prefix='fa'),
            popup=f"RCP ID: {idx}<br>Address: {row['adresse']}"
        ).add_to(existing_points)

    # Add New Potential Recycling Point Markers
    for _, row in closest_locations_gdf.iterrows():
        folium.Marker(
            location=[row.geometry.y, row.geometry.x],
            icon=folium.Icon(color='blue', icon='plus', prefix='fa'),
            popup=f"Potential Site ID: {row['potential_ID']}"
        ).add_to(new_recycling_points)

    # Add Flats with Duration-based Coloring
    for idx, row in flats_duration.iterrows():
        if not pd.isna(row['duration']):
            if row['duration'] > 15:
                color = outlier_color
                radius = 3
            else:
                color = colormap(row['duration'])
                radius = 1
            folium.CircleMarker(
                location=[row.geometry.y, row.geometry.x],
                radius=radius,
                color=color,
                fill=True,
                fill_color=color,
                fill_opacity=0.6,
                popup=f"Duration: {row['duration']:.2f} minutes, Population: {row['est_pop']:.2f}"
            ).add_to(m)

    # Assemble Map Layers
    flats_layer.add_to(m)
    existing_points.add_to(m)
    new_recycling_points.add_to(m)
    colormap.add_to(m)

    # Add Custom Legend
    legend_html = '''
    <div style="position: fixed; bottom: 50px; left: 50px; width: 180px; height: 130px; 
                border:2px solid grey; z-index:9999; font-size:14px; background-color:white;">
        &nbsp;<b>Legend</b><br>
        &nbsp;<i class="fa fa-circle fa-1x" style="color:red"></i>&nbsp; Duration > 15 min<br>
        &nbsp;<i class="fa fa-circle fa-1x" style="color:blue"></i>&nbsp; Flats (<= 15 min)<br>
        &nbsp;<i class="fa fa-recycle fa-1x" style="color:green"></i>&nbsp; Recycling Points<br>
        &nbsp;<i class="fa fa-plus fa-1x" style="color:blue"></i>&nbsp; New Recycling Points<br>
    </div>
    '''
    m.get_root().html.add_child(folium.Element(legend_html))
    folium.LayerControl().add_to(m)
    
    return m

def main():
    # Load and preprocess data
    flats, rcps, potential_sites = load_data()
    buildings = aggregate_flats(flats)
    
    # Determine route URL based on routing engine
    route_url = "http://localhost:8002/route"  # Default for Valhalla
    
    # Compute durations and export
    flats_duration = compute_nearest_durations(buildings, rcps, route_url)
    flats_duration.to_file(snakemake.output['flats_duration'], driver='GPKG')
    
    # Cluster flats and compute cluster centres
    cluster_centers, _ = cluster_flats(flats_duration)
    
    # Match the clusters to potential sites and merge with existing RCPs
    closest_locations_gdf = find_closest_potential(cluster_centers, potential_sites, route_url)
    merged_locations = merge_locations(rcps, closest_locations_gdf)
    merged_locations.to_file(snakemake.output['rcps_clustering_ors'], driver='GPKG')
    
    # Build and save the map
    m = build_map(rcps, closest_locations_gdf, flats_duration)
    m.save(snakemake.output['map_clustering_ors'])

if __name__ == "__main__":
    main()
