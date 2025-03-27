import geopandas as gpd
import pandas as pd
import numpy as np
import sys
import os
from snakemake.logging import logger

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import calculate_duration

# Get file paths from snakemake
POT_LOC_PATH = snakemake.input.potential_locations
DEMAND_POINTS_PATH = snakemake.input.demand_points
OUTPUT_MATRIX_PATH = snakemake.output.matrix_walking

# Get routing engine from params
ROUTING_ENGINE = snakemake.params.get('routing_engine', 'valhalla')

logger.info(f"Starting distance matrix calculation using {ROUTING_ENGINE}")

# Load data
potential_locations = gpd.read_file(POT_LOC_PATH).to_crs(epsg=4326)
demand_points = gpd.read_file(DEMAND_POINTS_PATH).to_crs(epsg=4326)

# Create empty dataframe for the distance matrix
facility_ids = potential_locations['ID'].tolist()
distance_matrix = pd.DataFrame(
    index=demand_points.index.astype(str),
    columns=facility_ids,
    dtype=float
)

# Calculate walking distances from each demand point to each facility
total_calculations = len(demand_points) * len(potential_locations)
current = 0

logger.info(f"Calculating {total_calculations} walking distances...")

for dp_idx, dp_row in demand_points.iterrows():
    dp_coord = (dp_row.geometry.x, dp_row.geometry.y)
    
    for fac_idx, fac_row in potential_locations.iterrows():
        fac_coord = (fac_row.geometry.x, fac_row.geometry.y)
        
        # Calculate walking duration using utility function
        duration = calculate_duration(dp_coord, fac_coord)
        
        # Store in distance matrix
        distance_matrix.loc[str(dp_idx), fac_row['ID']] = duration
        
        # Log progress
        current += 1
        if current % 100 == 0:
            logger.info(f"Processed {current}/{total_calculations} distances")

# Save distance matrix to CSV
distance_matrix.to_csv(OUTPUT_MATRIX_PATH)
logger.info(f"Distance matrix saved to {OUTPUT_MATRIX_PATH}")