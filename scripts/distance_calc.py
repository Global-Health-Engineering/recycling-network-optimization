import geopandas as gpd
import os
import pandas as pd
import openrouteservice as ors
import logging
import time
import numpy as np
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

# make sure crs matches
flats_zh = flats_zh.to_crs("EPSG:4326")
rcps = rcps.to_crs("EPSG:4326")

# aggregate flats to buildings
buildings_zh = flats_zh.groupby('egid').agg({'est_pop': 'sum', 'geometry': 'first'}).reset_index()

# Initialize BallTree
tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps)

# Initialize ORS client
client = ors.Client(base_url='http://localhost:8080/ors')

# apply the find_nearest_rcp_duration function to each building
buildings_zh['nearest_rcp_id'], buildings_zh['duration'] = zip(*buildings_zh['geometry'].apply(lambda geom: util.find_nearest_rcp_duration(geom, tree, rcp_coords, rcp_ids, client)))

# convert to geoDataFrame
buildings_zh = gpd.GeoDataFrame(buildings_zh, geometry='geometry', crs="EPSG:4326")

# calculate impact measure
buildings_zh['impact'] = buildings_zh['est_pop'] * buildings_zh['duration']
buildings_zh['impact_log'] = np.log1p(buildings_zh['impact'])

# Save the output
buildings_zh.to_file(OUTPUT_PATH, driver='GPKG')