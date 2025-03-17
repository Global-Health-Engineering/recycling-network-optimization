import geopandas as gpd
import openrouteservice as ors
import time
import numpy as np
import scripts.util as util

start_time = time.perf_counter()
FLATS_PATH = snakemake.input['flats']

# Read flats dataset and aggregate to buildings
flats_zh = gpd.read_file(FLATS_PATH)
flats_zh = flats_zh.to_crs("EPSG:4326")
buildings_zh = flats_zh.groupby('egid').agg({'est_pop': 'sum', 'geometry': 'first'}).reset_index()
buildings_agg = gpd.GeoDataFrame(buildings_zh, geometry='geometry', crs="EPSG:4326")

# List of rcp keys and corresponding output indices
rcp_keys = ['rcps1', 'rcps2', 'rcps3', 'rcps4']

for i, key in enumerate(rcp_keys):
    # Read the current rcp dataset and set CRS
    rcps = gpd.read_file(snakemake.input[key])
    rcps = rcps.to_crs("EPSG:4326")

    # Initialize BallTree for current rcps and ORS client
    if key == 'rcps1':
        id_column = 'poi_id'
    else:
        id_column = 'id' if 'id' in rcps.columns else 'ID'

    tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps, id_column)

    client = ors.Client(base_url='http://localhost:8080/ors')

    # Copy the aggregated buildings and compute nearest rcp duration
    buildings = buildings_agg.copy()
    buildings['nearest_rcp_id'], buildings['duration'] = zip(*buildings['geometry'].apply(
        lambda geom: util.find_nearest_rcp_duration(geom, tree, rcp_coords, rcp_ids, client)
    ))

    # Calculate impact measures
    buildings['impact'] = buildings['est_pop'] * buildings['duration']
    buildings['impact_log'] = np.log1p(buildings['impact'])

    # Save the output for current rcp dataset
    buildings.to_file(snakemake.output[i], driver='GPKG')