import geopandas as gpd
import requests
from snakemake.logging import logger
from shapely.ops import unary_union
from shapely.geometry import shape
import pandas as pd
import sys
import os

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import generate_isochrone, merge_isochrones_preserve_time

# Obtain paths from snakemake
INPUT_FLATS = snakemake.input['flats']
INPUT_RCPS = snakemake.input['rcps']
OUTPUT_GPKG_all = snakemake.output['iso_all']
OUTPUT_GPKG_merged = snakemake.output['iso_merged']

# Get routing engine from params
ROUTING_ENGINE = snakemake.params.get('routing_engine', 'valhalla')

TIME_LIMITS = [300, 600]  # in seconds

flats = gpd.read_file(INPUT_FLATS).to_crs("EPSG:4326")
rcps = gpd.read_file(INPUT_RCPS).to_crs("EPSG:4326")

logger.info(f"Started the script with {ROUTING_ENGINE} routing engine.")

def count_flats_in_isochrones(flats, gdf):
    try:
        flats = flats.to_crs("EPSG:4326")
        logger.info("Converted flats to EPSG:4326.")

        # Use 'within' predicate and ensure geometries are valid
        gdf['geometry'] = gdf['geometry'].buffer(0)  # Clean up any invalid geometries
        flats_in_isochrones = gpd.sjoin(
            flats, gdf, how='inner', predicate='within')
        est_pop = flats_in_isochrones.groupby('index_right')['est_pop'].sum()
        gdf['est_pop'] = est_pop
        gdf['est_pop'] = gdf['est_pop'].fillna(0).astype(int)

        logger.info("Added population estimation to isochrones.")
        logger.info(f"Total population across all isochrones: {gdf['est_pop'].sum()}")
    except Exception as e:
        logger.error(f"An error occurred while counting flats in isochrones: {e}")

# generate 1-10 min isochrones without population estimation
all_isochrones = []

times = [60, 120, 180, 240, 300, 360, 420, 480, 540, 600]

for time in times:
    for _, row in rcps.iterrows():
        isochrone = generate_isochrone([row.geometry.x, row.geometry.y], time)
        if isochrone is not None:
            # Extract features from the isochrone response
            for feature in isochrone.get('features', []):
                all_isochrones.append({
                    'geometry': shape(feature['geometry']),
                    'poi_id': row['poi_id'],
                    'time': time
                })

if all_isochrones:
    all_isochrones_gdf = gpd.GeoDataFrame(all_isochrones, crs="EPSG:4326")
    all_isochrones_gdf.to_file(OUTPUT_GPKG_all, driver="GPKG")
    logger.info("Generated all isochrones.")
else:
    logger.error(f"No isochrones generated - check if {ROUTING_ENGINE} isochrone service is running")

# Merge isochrones
if all_isochrones:
    merged_isochrones_gdf = merge_isochrones_preserve_time(all_isochrones_gdf)
    merged_isochrones_gdf.to_file(OUTPUT_GPKG_merged, driver="GPKG")
    logger.info("Merged isochrones saved to file.")
else:
    logger.error("Could not merge isochrones as none were generated")
