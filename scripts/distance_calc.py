import geopandas as gpd
import os
import pandas as pd
import openrouteservice as ors
import logging
import time
from snakemake.logging import logger
import sys
import scripts.util as util

start_time = time.perf_counter()

# Remove hardcoded parameters and paths
buffer_distance = snakemake.params.buffer_distance
n = snakemake.params.n if hasattr(snakemake.params, 'n') else None
FLATS_PATH = snakemake.input['flats']
RCPS_PATH = snakemake.input['rcps']
OUTPUT_PATH = snakemake.output[0]

# import datasets
flats_zh = gpd.read_file(FLATS_PATH)
rcps = gpd.read_file(RCPS_PATH)
# Initialize BallTree
tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps)

# Initialize ORS client
client = ors.Client(base_url='http://localhost:8080/ors')

# apply the find_nearest_rcp_duration function to each flat
flats_zh['nearest_rcp_id'], flats_zh['duration'] = zip(*flats_zh['geometry'].apply(func=lambda geom: util.find_nearest_rcp_duration(geom, tree, rcp_coords, rcp_ids, client)))


# Save the output
flats_zh.to_file(OUTPUT_PATH, driver='GPKG')