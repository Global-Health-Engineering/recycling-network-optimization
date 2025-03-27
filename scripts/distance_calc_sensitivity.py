import geopandas as gpd
import pandas as pd
import numpy as np
import sys
import os
from snakemake.logging import logger

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import initialize_ball_tree, find_nearest_rcp_duration, calculate_duration

# Get file paths from snakemake
INPUT_FLATS = snakemake.input.flats
INPUT_RCPS = snakemake.input.rcps
OUTPUT_PATH = snakemake.output.duration

# Get routing engine from params
ROUTING_ENGINE = snakemake.params.get('routing_engine', 'valhalla')

logger.info(f"Started calculating distances for sensitivity analysis using {ROUTING_ENGINE}")

# Load datasets
flats = gpd.read_file(INPUT_FLATS).to_crs(epsg=4326)
rcps = gpd.read_file(INPUT_RCPS).to_crs(epsg=4326)

# Initialize BallTree for fast nearest neighbor search
tree, rcp_coords, rcp_ids = initialize_ball_tree(rcps, 'id')

# Create a copy of flats dataframe
flats_output = flats.copy()

# Calculate duration to nearest RCP for each flat
durations = []
for idx, flat in flats_output.iterrows():
    nearest_id, duration = find_nearest_rcp_duration(flat.geometry, tree, rcp_coords, rcp_ids)
    durations.append({
        'flat_id': idx,
        'nearest_rcp': nearest_id,
        'duration_min': duration
    })
    
    if idx % 1000 == 0:
        logger.info(f"Processed {idx}/{len(flats_output)} flats")

# Create duration dataframe and join with flats
duration_df = pd.DataFrame(durations)
flats_output = flats_output.merge(duration_df, left_index=True, right_on='flat_id', how='left')

# Save to output file
flats_output.to_file(OUTPUT_PATH, driver='GPKG')
logger.info(f"Saved durations to {OUTPUT_PATH}")