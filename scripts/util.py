# util.py

import numpy as np
from sklearn.neighbors import BallTree
import requests
from shapely.geometry import Point, shape
import logging
import os
import yaml
import openrouteservice as ors
from openrouteservice.exceptions import ApiError
import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load configuration
def load_config():
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 
                              "config", "config.yaml")
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

config = load_config()
ROUTING_ENGINE = config.get("routing_engine", "valhalla")

# Initialize ORS client
def get_ors_client():
    """
    Initialize and return an OpenRouteService client
    """
    try:
        # Check if API key is available in config or environment
        api_key = config.get("ors", {}).get("api_key", "")
        if not api_key:
            api_key = os.environ.get('ORS_API_KEY', "")
            
        # Set the base URL (use default if not specified)
        base_url = config.get("ors", {}).get("route_url", "http://localhost:8080/ors")
        if base_url.startswith("'") and base_url.endswith("'"):
            base_url = base_url[1:-1]  # Remove quotes
            
        # If API key is provided, use public ORS instance
        if api_key:
            client = ors.Client(key=api_key)
        else:
            # Otherwise use local instance
            client = ors.Client(base_url=base_url)
        return client
    except Exception as e:
        logger.error(f"Failed to initialize ORS client: {e}")
        return None

# Get URLs based on selected routing engine
def get_route_url():
    if ROUTING_ENGINE == "valhalla":
        return config["valhalla"]["route_url"]
    else:
        return config["ors"]["route_url"]

def get_isochrone_url():
    if ROUTING_ENGINE == "valhalla":
        return config["valhalla"]["isochrone_url"]
    else:
        return config["ors"]["isochrone_url"]

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
    identifier_column = 'poi_id' if 'poi_id' in rcps.columns else ('ID' if 'ID' in rcps.columns else 'id')

    rcp_coords = rcps.geometry.apply(lambda geom: (geom.x, geom.y)).tolist()
    rcp_ids = rcps[identifier_column].tolist()

    # Convert coordinates to radians for BallTree ([lat, lon])
    rcp_coords_rad = np.radians([coord[::-1] for coord in rcp_coords])  # [lat, lon]

    # Build BallTree for efficient nearest neighbor search
    tree = BallTree(rcp_coords_rad, metric='haversine')

    logger.info("BallTree initialized with RCP coordinates.")
    return tree, rcp_coords, rcp_ids

def calculate_duration_valhalla(origin, destination, valhalla_url=None):
    """
    Calculate walking duration between origin and destination using Valhalla.

    Parameters:
    - origin: Tuple of (longitude, latitude)
    - destination: Tuple of (longitude, latitude)
    - valhalla_url: URL for the Valhalla routing service

    Returns:
    - Duration in minutes if successful, else None
    """
    if valhalla_url is None:
        valhalla_url = get_route_url()
        
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

def calculate_duration_ors(origin, destination, ors_url=None, api_key=None):
    """
    Calculate walking duration between two points using OpenRouteService Python client.
    
    Args:
        origin: (lon, lat) tuple or Point object
        destination: (lon, lat) tuple or Point object
        ors_url: URL for ORS API (used only if not using client)
        api_key: ORS API key (used only if not using client)
    
    Returns:
        Duration in minutes or None if calculation failed
    """
    try:
        # Convert Point objects to coordinate tuples if needed
        if isinstance(origin, Point):
            origin = (origin.x, origin.y)
        if isinstance(destination, Point):
            destination = (destination.x, destination.y)
        
        # Get ORS client
        client = get_ors_client()
        if client is None:
            logger.error("Failed to initialize ORS client")
            return None
            
        # Request directions using the Python client
        route = client.directions(
            coordinates=[origin, destination],
            profile='foot-walking',
            format='geojson',
            units='m',
            optimize_waypoints=False
        )
        
        # Extract and return duration in minutes
        duration_seconds = route['features'][0]['properties']['summary']['duration']
        return duration_seconds / 60  # Convert to minutes
        
    except ApiError as e:
        logger.error(f"ORS API error: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error calculating ORS duration: {e}")
        return None

def calculate_duration(origin, destination, client=None, route_url=None, api_key=None):
    """
    Calculate walking duration using the configured routing engine.
    """
    if ROUTING_ENGINE == "valhalla":
        return calculate_duration_valhalla(origin, destination, route_url)
    else:
        return calculate_duration_ors(origin, destination, route_url, api_key)

def find_nearest_rcp_duration(flat_geom, tree, rcp_coords, rcp_ids, valhalla_url="http://localhost:8002/route", radius=5000):
    """
    Find the nearest recycling point and calculate duration.
    Compatible with both Valhalla and ORS.
    """
    try:
        # Get flat coordinates
        if isinstance(flat_geom, Point):
            flat_coord = (flat_geom.x, flat_geom.y)
        else:
            flat_coord = (flat_geom.x, flat_geom.y)
        
        # Query the ball tree for nearest points
        distance, index = tree.query([[flat_coord[1], flat_coord[0]]], k=5)
        
        for dist, idx in zip(distance[0], index[0]):
            actual_distance = dist * 6371000  # Earth radius in meters
            if actual_distance <= radius:
                rcp_coord = rcp_coords[idx]
                duration = calculate_duration(flat_coord, rcp_coord, valhalla_url)
                if duration is not None:
                    return rcp_ids[idx], round(duration, 2)
        return None, None
        
    except Exception as e:
        logger.error(f"Error finding nearest RCP: {e}")
        return None, None

def generate_isochrone_valhalla(coord, time_limit, valhalla_url=None):
    """
    Generate an isochrone using Valhalla for a given coordinate and time limit.
    
    Parameters:
    - coord: Tuple of (longitude, latitude)
    - time_limit: Time in seconds
    - valhalla_url: URL for the Valhalla isochrone service
    
    Returns:
    - GeoJSON isochrone if successful, else None
    """
    if valhalla_url is None:
        valhalla_url = get_isochrone_url()
        
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

def generate_isochrone_ors(coord, time_limit, ors_url=None, api_key=None):
    """
    Generate an isochrone using OpenRouteService Python client.
    
    Args:
        coord: (lon, lat) tuple
        time_limit: time in seconds
        ors_url: URL for ORS API (used only if not using client)
        api_key: ORS API key (used only if not using client)
    
    Returns:
        GeoJSON response or None if generation failed
    """
    try:
        # Get ORS client
        client = get_ors_client()
        if client is None:
            logger.error("Failed to initialize ORS client")
            return None
        
        # Generate isochrone using the Python client
        isochrone = client.isochrones(
            locations=[[coord[0], coord[1]]],
            profile='foot-walking',
            range=[time_limit],  # time in seconds
            attributes=['area', 'reachfactor'],
            units='m',
            location_type='start',
        )
        
        return isochrone
        
    except ApiError as e:
        logger.error(f"ORS Isochrone API error: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error generating isochrone for {coord}: {e}")
        return None

def generate_isochrone(coord, time_limit, isochrone_url=None, api_key=None):
    """
    Generate an isochrone using the configured routing engine.
    """
    if ROUTING_ENGINE == "valhalla":
        return generate_isochrone_valhalla(coord, time_limit, isochrone_url)
    else:
        return generate_isochrone_ors(coord, time_limit, isochrone_url, api_key)

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