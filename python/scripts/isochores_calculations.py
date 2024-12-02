import geopandas as gpd
import openrouteservice
from openrouteservice import exceptions
import logging
import sys
import os
from shapely.geometry import shape

API_KEY = os.getenv("ORS_API_KEY")
INPUT_FLATS = "/home/silas/projects/msc_thesis/data/derived_data/flats_population.gpkg"
INPUT_RCPS = '/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp'
OUTPUT_GPKG_5 = "/home/silas/projects/msc_thesis/data/derived_data/isochores_5min.gpkg"
OUTPUT_GPKG_10 = "/home/silas/projects/msc_thesis/data/derived_data/isochores_10min.gpkg"
TIME_LIMITS = [300, 600]  # in seconds
n = 10 # number of recycling points to process

# Configure logging
logging.basicConfig(
    filename='isochores_calculations.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

logging.info("Started the script.")


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
        logging.error(f"API error for location {locations}: {e}")
    except Exception as e:
        logging.error(f"Unexpected error for location {locations}: {e}")
    return None

def count_flats_in_isochrones(flats, gdf):
    try:
        existing_flats = flats[flats['wstatlang'] == "Bestehend"]
        logging.info("Selected existing flats.")

        existing_flats = existing_flats.to_crs("EPSG:4326")
        logging.info("Converted flats to EPSG:4326.")

        # Use 'within' predicate and ensure geometries are valid
        gdf = gdf.buffer(0)  # Clean up any invalid geometries
        flats_in_isochrones = gpd.sjoin(existing_flats, gdf, how='inner', predicate='within')
        
        # Log diagnostic information
        logging.info(f"Total flats: {len(existing_flats)}")
        logging.info(f"Flats within isochrones: {len(flats_in_isochrones)}")
        
        est_pop = flats_in_isochrones.groupby('index_right')['est_pop'].sum()
        gdf['est_pop'] = est_pop
        gdf['est_pop'] = gdf['est_pop'].fillna(0).astype(int)
        gdf['poi_id'] = gdf.index
        
        logging.info("Added population estimation to isochrones.")
        logging.info(f"Total population across all isochrones: {gdf['est_pop'].sum()}")
    except Exception as e:
        logging.error(f"An error occurred while counting flats in isochrones: {e}")

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
                    'properties': cleaned_properties
                })
        logging.info(f"Processed recycling point {row['poi_id']} for {time_limit//60} min.")
    
    if iso:
        gdf = gpd.GeoDataFrame(iso, crs="EPSG:4326")
        
        # Import flat data
        flats = gpd.read_file(INPUT_FLATS)
        total_pop = flats['est_pop'].sum()
        logging.info(f"Imported flat data, population estimation: {total_pop}.")

        # Count flats within isochrones
        count_flats_in_isochrones(flats, gdf)
    
        gdf.to_file(output_path, driver="GPKG")
        logging.info(f"Saved {time_limit//60}-min isochrones with flat counts to {output_path}.")

def main():
    try:
        client = openrouteservice.Client(key=API_KEY)
        rcps = gpd.read_file(INPUT_RCPS)
        logging.info("Imported datasets.")
        rcps = rcps.to_crs(epsg=4326)

        # Generate 5-minute isochrones
        generate_and_save_isochrones(client, rcps, TIME_LIMITS[0], OUTPUT_GPKG_5)

        # Generate 10-minute isochrones
        generate_and_save_isochrones(client, rcps, TIME_LIMITS[1], OUTPUT_GPKG_10)

    except Exception as e:
        logging.critical(f"An unexpected error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()