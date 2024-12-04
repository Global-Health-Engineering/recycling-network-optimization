import openrouteservice
import os
from snakemake.logging import logger
import sys
import pandas as pd
import geopandas as gpd

# Constants
DUMP_COORDS = [8.512281878574365, 47.38447647508825]  # [longitude, latitude]
TRUCK_GARAGE_COORDS = [8.575500, 47.414889]  # [longitude, latitude]
source_coords = gpd.read_file(snakemake.input[0], header=None).geometry.tolist()[0]

def get_ors_client():
    """Initialize OpenRouteService client"""
    api_key = os.getenv("ORS_API_KEY")
    if not api_key:
        logger.error("OpenRouteService API key not found.")
        sys.exit(1)
    return openrouteservice.Client(key=api_key)

def calculate_distance_matrix(source_coords):
    """
    Calculate distance matrix using OpenRouteService API
    
    Args:
        source_coords (list): List of coordinates [lon, lat]
    
    Returns:
        dict: Distance matrix from OpenRouteService or None if error occurs
    """
    try:
        client = get_ors_client()
        matrix = client.distance_matrix(
            locations=[source_coords, DUMP_COORDS, TRUCK_GARAGE_COORDS],
            restrictions = {
            'profile':'driving-hgv',
            'restrictions': {
                'height': 4.0,
                'width': 2.5,
                'length': 17,
                'weight': 48,
                'axleload': 10.0,
                'hazmat': False,
                'hazmat_tunnel': False,
            }},
            metrics=['distance']
        )
        return matrix
    except openrouteservice.exceptions.ApiError as e:
        logger.error(f"API error: {e}")
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
    return None

# Execute script
try:
    os.chdir("/home/silas/projects/msc_thesis")
    logger.info("Changed working directory.")

    # Calculate distance matrix
    matrix = calculate_distance_matrix(source_coords)
    matrix.to_csv(snakemake.output[0])
    if matrix:
        logger.info("Calculated distance matrix.")
    else:
        logger.error("Failed to calculate distance matrix.")
finally:
    logger.info("Script finished.")
    