import geopandas as gpd
import openrouteservice
from openrouteservice import exceptions
import logging
import sys
import os
from shapely.geometry import shape

API_KEY = os.getenv("ORS_API_KEY") 
INPUT_FLATS ="/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp"
INPUT_RCPS = '/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp'
OUTPUT_GPKG_5 = "/home/silas/projects/msc_thesis/data/derived_data/isochores_5min.gpkg"
OUTPUT_GPKG_10 = "/home/silas/projects/msc_thesis/data/derived_data/isochores_10min.gpkg"
TIME_LIMITS = [300, 600]  # in seconds
n=10
# Removed os.chdir to use absolute paths

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

def main():
    try:
        client = openrouteservice.Client(key=API_KEY)
        rcps = gpd.read_file("/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp")
        logging.info("Imported datasets.")
        rcps = rcps.to_crs(epsg=4326)[1:n]
        
        iso_5 = []
        
        for _, row in rcps.iterrows():
            lon, lat = row.geometry.x, row.geometry.y
            isochrone_5 = generate_isochrones(client, [lon, lat], TIME_LIMITS[0])
            if isochrone_5:
                for feature in isochrone_5['features']:
                    properties = feature['properties']
                    # Entfernen von Feldern mit Listentypen
                    cleaned_properties = {k: v for k, v in properties.items() if not isinstance(v, list)}
                    iso_5.append({
                        'geometry': shape(feature['geometry']),
                        'properties': cleaned_properties
                    })
            logging.info(f"Processed recycling point {row['poi_id']}.")
        if iso_5:
            gdf = gpd.GeoDataFrame(iso_5, crs="EPSG:4326")
            gdf.to_file(OUTPUT_GPKG_5, driver="GPKG")
            logging.info(f"Saved 5-min isochrones to {OUTPUT_GPKG_5}.")

    except Exception as e:
        logging.critical(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()