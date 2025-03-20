import geopandas as gpd
import openrouteservice
import sys
import openrouteservice as ors
import os
from shapely.geometry import shape
from snakemake.logging import logger
from shapely.ops import unary_union
import pandas as pd

# Obtain paths from snakemake
INPUT_FLATS = snakemake.input['flats']
INPUT_RCPS = snakemake.input['rcps']
OUTPUT_GPKG_all= snakemake.output['iso_all']
OUTPUT_GPKG_merged = snakemake.output['iso_merged']

TIME_LIMITS = [300, 600]  # in seconds

client = ors.Client(base_url='http://localhost:8080/ors')
flats = gpd.read_file(INPUT_FLATS).to_crs("EPSG:4326")
rcps = gpd.read_file(INPUT_RCPS).to_crs("EPSG:4326")

logger.info("Started the script.")

def generate_isochrones(client, locations, time_limit):
    try:
        params = {
            "locations": [locations],
            "range": [time_limit],
            "range_type": "time",
            "units": "m",
            "location_type": "start",
            "smoothing": 0.3,
            "profile": "foot-walking",
        }
        isochrones = client.isochrones(**params)
        return isochrones
    except exceptions.ApiError as e:
        logger.error(f"API error for location {locations}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error for location {locations}: {e}")
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

def generate_and_save_isochrones(client, rcps, time_limit, output_path):
    iso = []
    for _, row in rcps.iterrows():
        lon, lat = row.geometry.x, row.geometry.y
        isochrone = generate_isochrones(client, [lon, lat], time_limit)
        if isochrone:
            for feature in isochrone['features']:
                iso.append({
                    'geometry': shape(feature['geometry']),
                    'poi_id': row['poi_id']
                })
        logger.info(f"Processed recycling point {row['poi_id']} for {time_limit//60} min.")
    
    if iso:
        gdf = gpd.GeoDataFrame(iso, crs="EPSG:4326")
        

        total_pop = flats['est_pop'].sum()
        logger.info(f"Imported flat data, population estimation: {total_pop}.")
        
        # Count flats within isochrones
        count_flats_in_isochrones(flats, gdf)
        
        gdf.to_file(output_path, driver="GPKG")
        logger.info(f"Saved {time_limit//60}-min isochrones with flat counts to {output_path}.")

# generate 1-10 min isochrones without population estimation
all_isochrones = []

times = [60, 120, 180, 240, 300, 360, 420, 480, 540, 600]

for time in times:
    for _, row in rcps.iterrows():
        isochrone = generate_isochrones(client, [row.geometry.x, row.geometry.y], time)
        if isochrone is not None:
            all_isochrones.extend([
                {
                    'geometry': shape(feature['geometry']),
                    'poi_id': row['poi_id'],
                    'time': time
                }
                for feature in isochrone['features']
            ])
if all_isochrones:
    all_isochrones_gdf = gpd.GeoDataFrame(all_isochrones, crs="EPSG:4326")
    all_isochrones_gdf.to_file(OUTPUT_GPKG_all, driver="GPKG")
    logger.info("Generated all isochrones.")
else:
    logger.error("unexpected error")

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
                logger.info(f"Processed isochrone with time {current_time} minutes.")

    return merged_isochrones

# Merge isochrones
merged_isochrones_gdf = merge_isochrones_preserve_time(all_isochrones_gdf)
merged_isochrones_gdf.to_file(OUTPUT_GPKG_merged, driver="GPKG")
