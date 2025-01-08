# util.py

import numpy as np
from sklearn.neighbors import BallTree
import openrouteservice
from shapely.geometry import Point
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def initialize_ball_tree(rcps, identifier_column='poi_id'):
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

def calculate_duration(origin, destination, client):
    """
    Calculate walking duration between origin and destination using ORS.

    Parameters:
    - origin: Tuple of (longitude, latitude)
    - destination: Tuple of (longitude, latitude)
    - client: OpenRouteService client instance

    Returns:
    - Duration in minutes if successful, else None
    """
    try:
        route = client.directions(
            coordinates=[origin, destination],
            profile='foot-walking',
            format='geojson'
        )
        duration_seconds = route['features'][0]['properties']['segments'][0]['duration']
        duration_minutes = duration_seconds / 60  # Convert to minutes
        return duration_minutes
    except openrouteservice.exceptions.ApiError as e:
        logger.error(f"ORS API error for origin {origin} to destination {destination}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error for origin {origin} to destination {destination}: {e}")
    return None

def find_nearest_rcp_duration(flat_geom, tree, rcp_coords, rcp_ids, client, radius=5000):
    """
    Find the nearest RCP within a specified radius and calculate walking duration.

    Parameters:
    - flat_geom: Shapely geometry of the flat
    - tree: BallTree instance with RCP coordinates
    - rcp_coords: List of RCP (longitude, latitude) tuples
    - rcp_ids: List of RCP identifiers
    - client: ORS client instance
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
            duration = calculate_duration(flat_coord, rcp_coord, client)
            if duration is not None:
                return rcp_ids[idx], round(duration, 2)
    return None, None


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

        # Initialize ORS client
        ors_key = os.getenv('ORS_API_KEY')
        if not ors_key:
            logger.error("OpenRouteService API key not found. Please set the 'ORS_API_KEY' environment variable.")
            return
        client = openrouteservice.Client(key=ors_key)

        # Example flat geometry (replace with actual data)
        example_flat_geom = Point(8.5417, 47.3769)  # Longitude, Latitude for Zurich

        # Find nearest RCP and duration
        rcp_id, duration = find_nearest_rcp_duration(example_flat_geom, tree, rcp_coords, rcp_ids, client)
        if rcp_id and duration:
            logger.info(f"Nearest RCP ID: {rcp_id}, Duration: {duration} minutes")
        else:
            logger.info("No RCP found within the specified radius.")

    main()