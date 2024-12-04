import openrouteservice
import os
from snakemake.logging import logger
import sys
import geopandas as gpd
import pandas as pd

# Constants
DUMP_COORDS = [8.512281878574365, 47.38447647508825]  # [longitude, latitude]
TRUCK_GARAGE_COORDS = [8.575500, 47.414889]  # [longitude, latitude]
source_coords = gpd.read_file(snakemake.input[0]).geometry.tolist()

def get_ors_client():
    """Initialize OpenRouteService client"""
    api_key = os.getenv("ORS_API_KEY")
    if not api_key:
        logger.error("OpenRouteService API key not found.")
        sys.exit(1)
    return openrouteservice.Client(key=api_key)

def calculate_distance_matrix(coordinates):
    try:
        # Convert shapely Point objects to [longitude, latitude] lists
        coords_list = [[p.x, p.y] for p in coordinates]
        matrix = client.distance_matrix(
            locations=coords_list,
            metrics=['distance'],
            profile='driving-car',
        )
        # Convert the matrix to a DataFrame for CSV export
        distances = matrix.get('distances', [])
        df = pd.DataFrame(distances)
        return df
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {e}")
        return None

# Initialize OpenRouteService client
client = get_ors_client()

# Execute script
try:
    os.chdir("/home/silas/projects/msc_thesis")
    logger.info("Changed working directory.")

    # Calculate distance matrix
    matrix = calculate_distance_matrix(source_coords)
    if matrix is not None:
        matrix.to_csv(snakemake.output[0], index=False)
        logger.info("Calculated and saved distance matrix.")
    else:
        logger.error("Failed to calculate distance matrix.")
finally:
    logger.info("Script finished.")