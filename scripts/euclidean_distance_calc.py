import geopandas as gpd
import pandas as pd
import numpy as np
import sys
import os
from snakemake.logging import logger
from shapely.geometry import Point
from scipy.spatial.distance import cdist
from sklearn.neighbors import BallTree

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Get file paths from snakemake
INPUT_FLATS = snakemake.input.flats
INPUT_RCPS1 = snakemake.input.rcps1
INPUT_RCPS2 = snakemake.input.rcps2
INPUT_RCPS3 = snakemake.input.rcps3

OUTPUT1 = snakemake.output[0]  # flats_duration_clustering_iso_euclidean.gpkg
OUTPUT2 = snakemake.output[1]  # flats_duration_clustering_ors_euclidean.gpkg
OUTPUT3 = snakemake.output[2]  # flats_duration_opt_euclidean.gpkg

# Get walking speed from params (km/h)
WALKING_SPEED_KMH = snakemake.params.walking_speed_kmh

logger.info(f"Started calculating Euclidean distances with walking speed: {WALKING_SPEED_KMH} km/h")

def calculate_euclidean_distance(point1, point2, crs_epsg=2056):
    """
    Calculate Euclidean distance between two points.
    
    Parameters:
    - point1: Shapely Point geometry
    - point2: Shapely Point geometry  
    - crs_epsg: EPSG code for projected coordinate system (default: 2056 for Swiss coordinates)
    
    Returns:
    - Distance in meters
    """
    # Create temporary GeoDataFrames to handle CRS transformation
    gdf1 = gpd.GeoDataFrame([1], geometry=[point1], crs="EPSG:4326")
    gdf2 = gpd.GeoDataFrame([1], geometry=[point2], crs="EPSG:4326")
    
    # Transform to projected CRS for accurate distance calculation
    gdf1_proj = gdf1.to_crs(f"EPSG:{crs_epsg}")
    gdf2_proj = gdf2.to_crs(f"EPSG:{crs_epsg}")
    
    # Calculate Euclidean distance
    point1_proj = gdf1_proj.geometry.iloc[0]
    point2_proj = gdf2_proj.geometry.iloc[0]
    
    distance = np.sqrt((point1_proj.x - point2_proj.x)**2 + (point1_proj.y - point2_proj.y)**2)
    return distance

def calculate_euclidean_distance_fast(point1_coords, point2_coords, crs_epsg=2056):
    """
    Fast Euclidean distance calculation using coordinate tuples.
    This avoids creating temporary GeoDataFrames for each calculation.
    
    Parameters:
    - point1_coords: Tuple of (longitude, latitude) for point1
    - point2_coords: Tuple of (longitude, latitude) for point2
    - crs_epsg: EPSG code for projected coordinate system (default: 2056 for Swiss coordinates)
    
    Returns:
    - Distance in meters
    """
    # Create points and transform once
    point1 = Point(point1_coords)
    point2 = Point(point2_coords)
    
    # Create temporary GeoDataFrames for transformation
    gdf1 = gpd.GeoDataFrame([1], geometry=[point1], crs="EPSG:4326")
    gdf2 = gpd.GeoDataFrame([1], geometry=[point2], crs="EPSG:4326")
    
    # Transform to projected CRS
    gdf1_proj = gdf1.to_crs(f"EPSG:{crs_epsg}")
    gdf2_proj = gdf2.to_crs(f"EPSG:{crs_epsg}")
    
    # Calculate distance
    p1_proj = gdf1_proj.geometry.iloc[0]
    p2_proj = gdf2_proj.geometry.iloc[0]
    
    distance = np.sqrt((p1_proj.x - p2_proj.x)**2 + (p1_proj.y - p2_proj.y)**2)
    return distance

def euclidean_distance_to_duration(distance_m, walking_speed_kmh):
    """
    Convert Euclidean distance to walking duration.
    
    Parameters:
    - distance_m: Distance in meters
    - walking_speed_kmh: Walking speed in km/h
    
    Returns:
    - Duration in minutes
    """
    distance_km = distance_m / 1000.0
    duration_hours = distance_km / walking_speed_kmh
    duration_minutes = duration_hours * 60.0
    return duration_minutes

def initialize_ball_tree_euclidean(rcps_gdf):
    """
    Initialize a BallTree for RCP coordinates - adapted from util.py
    
    Parameters:
    - rcps_gdf: GeoDataFrame containing RCPs
    
    Returns:
    - tree: BallTree object for nearest neighbor searches
    - rcp_coords: List of RCP coordinates as (longitude, latitude) tuples  
    - rcp_ids: List of RCP identifiers
    """
    # Determine the correct identifier column
    if 'poi_id' in rcps_gdf.columns:
        identifier_column = 'poi_id'
    elif 'ID' in rcps_gdf.columns:
        identifier_column = 'ID'
    elif 'id' in rcps_gdf.columns:
        identifier_column = 'id'
    else:
        identifier_column = rcps_gdf.columns[0]
        logger.warning(f"No standard ID column found, using {identifier_column}")
    
    rcp_coords = rcps_gdf.geometry.apply(lambda geom: (geom.x, geom.y)).tolist()
    rcp_ids = rcps_gdf[identifier_column].tolist()
    
    # Convert coordinates to radians for BallTree ([lat, lon])
    rcp_coords_rad = np.radians([coord[::-1] for coord in rcp_coords])  # [lat, lon]
    
    # Build BallTree for efficient nearest neighbor search
    tree = BallTree(rcp_coords_rad, metric='haversine')
    
    return tree, rcp_coords, rcp_ids

def find_nearest_rcp_euclidean(flat_geom, tree, rcp_coords, rcp_ids):
    """
    Find the nearest RCP using BallTree and calculate Euclidean walking duration.
    
    Parameters:
    - flat_geom: Shapely geometry of the flat
    - tree: BallTree instance with RCP coordinates
    - rcp_coords: List of RCP (longitude, latitude) tuples
    - rcp_ids: List of RCP identifiers
    
    Returns:
    - Tuple of (rcp_id, duration_min)
    """
    flat_coord = (flat_geom.x, flat_geom.y)
    flat_rad = np.radians([flat_coord[::-1]])  # [lat, lon] in radians
    
    # Find the nearest RCP using BallTree
    distance, index = tree.query(flat_rad, k=1)
    nearest_idx = index[0][0]
    
    # Get the nearest RCP coordinates and ID
    nearest_rcp_coord = rcp_coords[nearest_idx]
    nearest_rcp_id = rcp_ids[nearest_idx]
    
    # Calculate Euclidean distance using the faster coordinate-based function
    flat_coord = (flat_geom.x, flat_geom.y)
    distance = calculate_euclidean_distance_fast(flat_coord, nearest_rcp_coord)
    duration = euclidean_distance_to_duration(distance, WALKING_SPEED_KMH)
    
    return nearest_rcp_id, round(duration, 2)

# Load datasets
logger.info("Loading datasets...")
flats = gpd.read_file(INPUT_FLATS).to_crs(epsg=4326)
rcps1 = gpd.read_file(INPUT_RCPS1).to_crs(epsg=4326)
rcps2 = gpd.read_file(INPUT_RCPS2).to_crs(epsg=4326)
rcps3 = gpd.read_file(INPUT_RCPS3).to_crs(epsg=4326)

logger.info(f"Loaded {len(flats)} flats")
logger.info(f"Loaded {len(rcps1)} RCPs (clustering_iso)")
logger.info(f"Loaded {len(rcps2)} RCPs (clustering_ors)")
logger.info(f"Loaded {len(rcps3)} RCPs (optimization)")

# Create output directory if it doesn't exist
os.makedirs(os.path.dirname(OUTPUT1), exist_ok=True)

# Process each method
for rcps, output_path, method_name in [
    (rcps1, OUTPUT1, "clustering_iso"),
    (rcps2, OUTPUT2, "clustering_ors"),
    (rcps3, OUTPUT3, "optimization")
]:
    logger.info(f"Processing method: {method_name}")
    
    # Initialize BallTree for fast nearest neighbor search
    tree, rcp_coords, rcp_ids = initialize_ball_tree_euclidean(rcps)
    logger.info(f"Initialized BallTree with {len(rcp_coords)} RCPs for {method_name}")
    
    # Create a copy of flats dataframe for this method
    flats_method = flats.copy()
    
    # Calculate Euclidean duration to nearest RCP for each flat
    durations = []
    total_flats = len(flats_method)
    
    for idx, flat in flats_method.iterrows():
        nearest_id, duration = find_nearest_rcp_euclidean(flat.geometry, tree, rcp_coords, rcp_ids)
        durations.append({
            'flat_id': idx,
            'nearest_rcp': nearest_id,
            'duration_min': duration
        })
        
        if (idx + 1) % 500 == 0:
            logger.info(f"Processed {idx + 1}/{total_flats} flats for {method_name}")
    
    # Create duration dataframe and join with flats
    duration_df = pd.DataFrame(durations)
    flats_method = flats_method.merge(duration_df, left_index=True, right_on='flat_id', how='left')
    
    # Drop the redundant flat_id column
    flats_method = flats_method.drop('flat_id', axis=1)
    
    # Save to output file
    flats_method.to_file(output_path, driver='GPKG')
    logger.info(f"Saved Euclidean durations for {method_name} to {output_path}")
    
    # Log some statistics
    logger.info(f"Method {method_name} statistics:")
    logger.info(f"  Mean duration: {flats_method['duration_min'].mean():.2f} minutes")
    logger.info(f"  Median duration: {flats_method['duration_min'].median():.2f} minutes")
    logger.info(f"  Max duration: {flats_method['duration_min'].max():.2f} minutes")
    logger.info(f"  Min duration: {flats_method['duration_min'].min():.2f} minutes")

logger.info("All Euclidean distance calculations completed")
