import geopandas as gpd
import pandas as pd
import sys
import os
from snakemake.logging import logger

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import initialize_ball_tree, find_nearest_rcp_duration

# Get file paths from snakemake
INPUT_FLATS = snakemake.input.flats
INPUT_RCPS1 = snakemake.input.rcps1
INPUT_RCPS2 = snakemake.input.rcps2
INPUT_RCPS3 = snakemake.input.rcps3
OUTPUT1 = snakemake.output.output1
OUTPUT2 = snakemake.output.output2
OUTPUT3 = snakemake.output.output3

# Get routing engine from params
ROUTING_ENGINE = snakemake.params.get('routing_engine', 'valhalla')

logger.info(f"Started calculating distances using {ROUTING_ENGINE} routing engine")

# Load datasets
flats = gpd.read_file(INPUT_FLATS).to_crs(epsg=4326)
rcps1 = gpd.read_file(INPUT_RCPS1).to_crs(epsg=4326)
rcps2 = gpd.read_file(INPUT_RCPS2).to_crs(epsg=4326)
rcps3 = gpd.read_file(INPUT_RCPS3).to_crs(epsg=4326)

# Process each method
for rcps, output_path, method_name in [
    (rcps1, OUTPUT1, "clustering_iso"),
    (rcps2, OUTPUT2, "clustering_ors"),
    (rcps3, OUTPUT3, "optimization")
]:
    logger.info(f"Processing method: {method_name}")
    
    # Initialize BallTree for fast nearest neighbor search
    tree, rcp_coords, rcp_ids = initialize_ball_tree(rcps, 'id')
    
    # Create a copy of flats dataframe for this method
    flats_method = flats.copy()
    
    # Calculate duration to nearest RCP for each flat
    durations = []
    for idx, flat in flats_method.iterrows():
        # No need to specify valhalla_url - handled by util.py
        nearest_id, duration = find_nearest_rcp_duration(flat.geometry, tree, rcp_coords, rcp_ids)
        durations.append({
            'flat_id': idx,
            'nearest_rcp': nearest_id,
            'duration_min': duration
        })
        
        if idx % 1000 == 0:
            logger.info(f"Processed {idx}/{len(flats_method)} flats for {method_name}")
    
    # Create duration dataframe and join with flats
    duration_df = pd.DataFrame(durations)
    flats_method = flats_method.merge(duration_df, left_index=True, right_on='flat_id', how='left')
    
    # Save to output file
    flats_method.to_file(output_path, driver='GPKG')
    logger.info(f"Saved durations for {method_name} to {output_path}")

logger.info("All distance calculations completed")