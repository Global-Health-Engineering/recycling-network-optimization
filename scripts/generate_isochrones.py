import geopandas as gpd
import sys
import os
import requests
from snakemake.logging import logger
from shapely.ops import unary_union
from shapely.geometry import shape
import pandas as pd

# Obtain paths from snakemake
INPUT_FLATS = snakemake.input['flats']
INPUT_RCPS = snakemake.input['rcps']
OUTPUT_GPKG_all= snakemake.output['iso_all']
OUTPUT_GPKG_merged = snakemake.output['iso_merged']

TIME_LIMITS = [300, 600]  # in seconds
VALHALLA_ISOCHRONE_URL = "http://localhost:8002/isochrone"

flats = gpd.read_file(INPUT_FLATS).to_crs("EPSG:4326")
rcps = gpd.read_file(INPUT_RCPS).to_crs("EPSG:4326")

logger.info("Started the script.")

def generate_isochrones(locations, time_limit, valhalla_url=VALHALLA_ISOCHRONE_URL):
    try:
        # Prepare Valhalla isochrone request
        valhalla_params = {
            "locations": [{"lat": locations[1], "lon": locations[0]}],
            "costing": "pedestrian",
            "contours": [{"time": time_limit}],
            "polygons": True,
            "denoise": 0.5,
            "generalize": 50,
            "format": "json"
        }
        
        # Make the request to Valhalla
        response = requests.post(valhalla_url, json=valhalla_params)
        
        if response.status_code == 200:
            return response.json()
        else:
            logger.error(f"Valhalla Isochrone API error: Status code {response.status_code}")
            logger.error(response.text)
            return None
            
    except Exception as e:
        logger.error(f"Unexpected error generating isochrone for {locations}: {e}")
        return None

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

def merge_isochrones_preserve_time(isochrones_gdf):
    """
    Merge isochrones preserving lower time values.

    Parameters:
    - isochrones_gdf: GeoDataFrame with isochrones and 'time' attribute.

    Returns:
    - GeoDataFrame with merged isochrones.
    """
    # Ensure CRS is EPSG:4326
    if isochrones_gdf.crs != "EPSG:4326":
        isochrones_gdf = isochrones_gdf.to_crs(epsg=4326)

    # Sort isochrones by 'time' ascending
    isochrones_sorted = isochrones_gdf.sort_values(by='time')

    merged_isochrones = gpd.GeoDataFrame(columns=isochrones_sorted.columns, crs="EPSG:4326")

    # Initialize an empty geometry for subtraction
    accumulated_geom = None

    for _, row in isochrones_sorted.iterrows():
        current_geom = row.geometry
        current_time = row['time']

        if accumulated_geom:
            remaining_geom = current_geom.difference(accumulated_geom)
        else:
            remaining_geom = current_geom

        if not remaining_geom.is_empty:
            new_row = row.copy()
            new_row.geometry = remaining_geom
            # Ensure the new_row GeoDataFrame has the correct CRS
            new_row = gpd.GeoDataFrame([new_row], crs="EPSG:4326")
            merged_isochrones = pd.concat([merged_isochrones, new_row], ignore_index=True)
            # Update accumulated geometry
            if accumulated_geom:
                accumulated_geom = unary_union([accumulated_geom, remaining_geom])
            else:
                accumulated_geom = remaining_geom
                logger.info(f"Processed isochrone with time {current_time} minutes.")

    return merged_isochrones

# generate 1-10 min isochrones without population estimation
all_isochrones = []

times = [60, 120, 180, 240, 300, 360, 420, 480, 540, 600]

for time in times:
    for _, row in rcps.iterrows():
        isochrone = generate_isochrones([row.geometry.x, row.geometry.y], time)
        if isochrone is not None:
            # Extract features from the Valhalla isochrone response
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
    logger.error("No isochrones generated - check if Valhalla isochrone service is running")

# Merge isochrones
if all_isochrones:
    merged_isochrones_gdf = merge_isochrones_preserve_time(all_isochrones_gdf)
    merged_isochrones_gdf.to_file(OUTPUT_GPKG_merged, driver="GPKG")
    logger.info("Merged isochrones saved to file.")
else:
    logger.error("Could not merge isochrones as none were generated")
