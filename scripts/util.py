# util.py

import numpy as np
from sklearn.neighbors import BallTree
import requests
from shapely.geometry import Point
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def initialize_ball_tree(rcps, identifier_column):
    """
    Initialize a BallTree for RCP coordinates.

    Parameters:
    - rcps: GeoDataFrame containing RCPs with 'geometry' and 'poi_id' columns.

    Returns:
    - tree: BallTree object for nearest neighbor searches.
    - rcp_coords: List of RCP coordinates as (longitude, latitude) tuples.
    - rcp_ids: List of RCP identifiers.
    """
    rcp_coords = rcps.geometry.apply(lambda geom: (geom.x, geom.y)).tolist()
    rcp_ids = rcps[identifier_column].tolist()

    # Convert coordinates to radians for BallTree ([lat, lon])
    rcp_coords_rad = np.radians([coord[::-1] for coord in rcp_coords])  # [lat, lon]

    # Build BallTree for efficient nearest neighbor search
    tree = BallTree(rcp_coords_rad, metric='haversine')

    logger.info("BallTree initialized with RCP coordinates.")
    return tree, rcp_coords, rcp_ids

def calculate_duration_valhalla(origin, destination, valhalla_url="http://localhost:8002/route"):
    """
    Calculate walking duration between origin and destination using Valhalla.

    Parameters:
    - origin: Tuple of (longitude, latitude)
    - destination: Tuple of (longitude, latitude)
    - valhalla_url: URL for the Valhalla routing service

    Returns:
    - Duration in minutes if successful, else None
    """
    try:
        # Prepare Valhalla request
        valhalla_params = {
            "locations": [
                {"lat": origin[1], "lon": origin[0]},
                {"lat": destination[1], "lon": destination[0]}
            ],
            "costing": "pedestrian",
            "directions_options": {
                "units": "kilometers"
            },
            "format": "json"
        }
        
        # Make the request to Valhalla
        response = requests.post(valhalla_url, json=valhalla_params)
        
        if response.status_code == 200:
            data = response.json()
            duration_seconds = data["trip"]["legs"][0]["summary"]["time"]
            duration_minutes = duration_seconds / 60  # Convert to minutes
            return duration_minutes
        else:
            logger.error(f"Valhalla API error: Status code {response.status_code}")
            logger.error(response.text)
            return None
            
    except Exception as e:
        logger.error(f"Unexpected error for origin {origin} to destination {destination}: {e}")
        return None

def find_nearest_rcp_duration(flat_geom, tree, rcp_coords, rcp_ids, valhalla_url="http://localhost:8002/route", radius=5000):
    """
    Find the nearest RCP within a specified radius and calculate walking duration.

    Parameters:
    - flat_geom: Shapely geometry of the flat
    - tree: BallTree instance with RCP coordinates
    - rcp_coords: List of RCP (longitude, latitude) tuples
    - rcp_ids: List of RCP identifiers
    - valhalla_url: URL for the Valhalla routing service
    - radius: Search radius in meters (default: 5000)

    Returns:
    - Tuple of (rcp_id, duration_min)
    """
    flat_coord = (flat_geom.x, flat_geom.y)
    flat_rad = np.radians([flat_coord[::-1]])  # [lat, lon] in radians

    # Query for the top 5 nearest RCPs
    distance, index = tree.query(flat_rad, k=5)

    for dist, idx in zip(distance[0], index[0]):
        actual_distance = dist * 6371000  # Earth radius in meters
        if actual_distance <= radius:
            rcp_coord = rcp_coords[idx]
            duration = calculate_duration_valhalla(flat_coord, rcp_coord, valhalla_url)
            if duration is not None:
                return rcp_ids[idx], round(duration, 2)
    return None, None

def generate_isochrone_valhalla(coord, time_limit, valhalla_url="http://localhost:8002/isochrone"):
    """
    Generate an isochrone using Valhalla for a given coordinate and time limit.
    
    Parameters:
    - coord: Tuple of (longitude, latitude)
    - time_limit: Time in seconds
    - valhalla_url: URL for the Valhalla isochrone service
    
    Returns:
    - GeoJSON isochrone if successful, else None
    """
    try:
        # Prepare Valhalla isochrone request
        valhalla_params = {
            "locations": [{"lat": coord[1], "lon": coord[0]}],
            "costing": "pedestrian",
            "contours": [{"time": time_limit/60}],
            "polygons": True,
            "denoise": 0.5,
            "generalize": 50,
            "format": "json"
        }
        
        # Make the request to Valhalla
        response = requests.post(valhalla_url, json=valhalla_params)
        
        if response.status_code == 200:
            return response.json()
        else:
            logger.error(f"Valhalla Isochrone API error: Status code {response.status_code}")
            logger.error(response.text)
            return None
            
    except Exception as e:
        logger.error(f"Unexpected error generating isochrone for {coord}: {e}")
        return None

# For backward compatibility, keep functions with old names but use the new implementation
def calculate_duration(origin, destination, client=None, valhalla_url="http://localhost:8002/route"):
    """
    Backward compatibility function that uses Valhalla instead of ORS.
    """
    return calculate_duration_valhalla(origin, destination, valhalla_url)

# Placeholder to maintain compatibility with the merge_isochrones_preserve_time function
import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union

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

if __name__ == "__main__":
    def main():
        """
        Main function to test utility functions.
        """
        import geopandas as gpd

        # Example paths (update these paths as needed)
        rcp_shapefile_path = '/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp'

        # Load RCPs GeoDataFrame
        rcps = gpd.read_file(rcp_shapefile_path)
        rcps = rcps.to_crs("EPSG:4326")  # Ensure CRS is WGS84

        # Initialize BallTree
        tree, rcp_coords, rcp_ids = initialize_ball_tree(rcps)

        # Example flat geometry (replace with actual data)
        example_flat_geom = Point(8.5417, 47.3769)  # Longitude, Latitude for Zurich

        # Find nearest RCP and duration
        rcp_id, duration = find_nearest_rcp_duration(example_flat_geom, tree, rcp_coords, rcp_ids)
        if rcp_id and duration:
            logger.info(f"Nearest RCP ID: {rcp_id}, Duration: {duration} minutes")
        else:
            logger.info("No RCP found within the specified radius.")

    main()