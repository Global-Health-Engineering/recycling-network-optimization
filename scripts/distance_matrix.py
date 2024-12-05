import openrouteservice
import os
from snakemake.logging import logger
import sys
import geopandas as gpd
import pandas as pd

# Constants
DUMP_COORDS = [8.512281878574365, 47.38447647508825]  # [longitude, latitude]
TRUCK_GARAGE_COORDS = [8.575500, 47.414889]  # [longitude, latitude]
DEPOT_COORDS = [DUMP_COORDS, TRUCK_GARAGE_COORDS]

def get_ors_client():
    """Initialize OpenRouteService client"""
    api_key = os.getenv("ORS_API_KEY")
    if not api_key:
        logger.error("OpenRouteService API key not found.")
        sys.exit(1)
    return openrouteservice.Client(key=api_key)

def calculate_distance_matrix(source_coords, depot_coords):
    """Calculate distance matrix between source coordinates and depot coordinates"""
    try:
        # Convert source coordinates from Points to list format
        source_coords_list = [[p.x, p.y] for p in source_coords]
        
        # depot_coords is already in the correct format
        depot_coords_list = depot_coords  # Already [longitude, latitude] pairs
        
        # Get distances from ORS
        distances = []
        for coord in source_coords_list:
            response = client.directions(
                coordinates=[coord, depot_coords_list[0]],  # To dump
                profile='driving-car',
                format='geojson'
            )
            dump_dist = response['features'][0]['properties']['segments'][0]['distance']
            
            response = client.directions(
                coordinates=[coord, depot_coords_list[1]],  # To garage
                profile='driving-car',
                format='geojson'
            )
            garage_dist = response['features'][0]['properties']['segments'][0]['distance']
            
            distances.append([dump_dist, garage_dist])

        # Create DataFrame
        df = pd.DataFrame(distances, columns=['Distance_to_Dump', 'Distance_to_Garage'])
        return df
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {e}")
        return None

# Main execution
try:
    os.chdir("/home/silas/projects/msc_thesis")
    logger.info("Changed working directory.")
    
    # Read and ensure CRS matches EPSG:4326
    gdf = gpd.read_file(snakemake.input[0])
    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
        logger.info("Converted CRS to EPSG:4326.")
    
    source_coords = gdf.geometry.tolist()
    
    # Initialize client
    client = get_ors_client()
    
    # Calculate matrix
    matrix = calculate_distance_matrix(source_coords, DEPOT_COORDS)
    if matrix is not None:
        matrix.to_csv(snakemake.output[0], index=False)
        logger.info("Calculated and saved distance matrix.")
    else:
        logger.error("Failed to calculate distance matrix.")
finally:
    logger.info("Script finished.")