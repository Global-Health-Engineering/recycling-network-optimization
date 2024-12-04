import geopandas as gpd
import openrouteservice
from openrouteservice import exceptions
import sys
import os
from shapely.geometry import shape
from snakemake.logging import logger

# Obtain paths from snakemake
INPUT_FLATS = snakemake.input['flats']
INPUT_RCPS = snakemake.input['rcps']
OUTPUT_GPKG_5 = snakemake.output['iso_5min']
OUTPUT_GPKG_10 = snakemake.output['iso_10min']

# Get API key from environment variable or snakemake params
API_KEY = os.getenv("ORS_API_KEY")
if not API_KEY:
    logger.error("OpenRouteService API key not found.")
    sys.exit(1)

TIME_LIMITS = [300, 600]  # in seconds
n = 10  # number of recycling points to process

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
        existing_flats = flats[flats['wstatlang'] == "Bestehend"]
        logger.info("Selected existing flats.")

        existing_flats = existing_flats.to_crs("EPSG:4326")
        logger.info("Converted flats to EPSG:4326.")

        # Use 'within' predicate and ensure geometries are valid
        gdf['geometry'] = gdf['geometry'].buffer(0)  # Clean up any invalid geometries
        existing_flats = existing_flats.drop(columns=['index_left', 'index_right'], errors='ignore')
        gdf = gdf.drop(columns=['index_left', 'index_right'], errors='ignore')
        flats_in_isochrones = gpd.sjoin(existing_flats, gdf, how='inner', predicate='within')
        
        est_pop = flats_in_isochrones.groupby('index_right')['est_pop'].sum()
        gdf['est_pop'] = est_pop
        gdf['est_pop'] = gdf['est_pop'].fillna(0).astype(int)
        gdf['poi_id'] = gdf.index
        
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
                properties = feature['properties']
                cleaned_properties = {k: v for k, v in properties.items() if not isinstance(v, list)}
                iso.append({
                    'geometry': shape(feature['geometry']),
                    'poi_id': row['poi_id'],
                    'properties': cleaned_properties,
                    'est_pop': 0
                })
        logger.info(f"Processed recycling point {row['poi_id']} for {time_limit//60} min.")
    
    if iso:
        gdf = gpd.GeoDataFrame(iso, crs="EPSG:4326")
        
        flats = gpd.read_file(INPUT_FLATS)
        if 'est_pop' not in flats.columns:
            flats['est_pop'] = 0
        total_pop = flats['est_pop'].sum()
        logger.info(f"Imported flat data, population estimation: {total_pop}.")
        
        # Count flats within isochrones
        count_flats_in_isochrones(flats, gdf)
        
        gdf.to_file(output_path, driver="GPKG")
        logger.info(f"Saved {time_limit//60}-min isochrones with flat counts to {output_path}.")

try:
    client = openrouteservice.Client(key=API_KEY)
    rcps = gpd.read_file(INPUT_RCPS)
    logger.info("Imported datasets.")
    rcps = rcps.to_crs(epsg=4326)

    # Generate 5-minute isochrones
    generate_and_save_isochrones(client, rcps, TIME_LIMITS[0], OUTPUT_GPKG_5)

    # Generate 10-minute isochrones
    generate_and_save_isochrones(client, rcps, TIME_LIMITS[1], OUTPUT_GPKG_10)

except Exception as e:
    logger.critical(f"An unexpected error occurred: {e}")
    sys.exit(1)
