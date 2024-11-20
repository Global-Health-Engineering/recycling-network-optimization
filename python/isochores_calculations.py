import geopandas as gpd
import openrouteservice
from openrouteservice import exceptions
import logging
import sys

# Parameters
API_KEY = "your_openrouteservice_api_key"
INPUT_RCPS = "./data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
OUTPUT_GPKG_5 = "./data/derived_data/isochores_5min.gpkg"
OUTPUT_GPKG_10 = "./data/derived_data/isochores_10min.gpkg"
TIME_LIMITS = [300, 600]  # in seconds

# Configure logging
logging.basicConfig(
    filename='isochores_calculations.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def generate_isochrones(client, locations, time_limit):
    try:
        params = {
            "locations": [locations],
            "range": [time_limit],
            "units": "m",
            "location_type": "start",
            "smoothing": True
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
        rcps = gpd.read_file(INPUT_RCPS)
        rcps = rcps.to_crs(epsg=4326)
        
        iso_5 = []
        iso_10 = []
        
        for idx, row in rcps.iterrows():
            lon, lat = row.geometry.x, row.geometry.y
            isochrone_5 = generate_isochrones(client, [lon, lat], TIME_LIMITS[0])
            isochrone_10 = generate_isochrones(client, [lon, lat], TIME_LIMITS[1])
            
            if isochrone_5:
                iso_5.append(isochrone_5)
            if isochrone_10:
                iso_10.append(isochrone_10)
        
        if iso_5:
            gpd.GeoDataFrame.from_features(iso_5).to_file(OUTPUT_GPKG_5, driver="GPKG")
            logging.info("Saved 5-minute isochores.")
        
        if iso_10:
            gpd.GeoDataFrame.from_features(iso_10).to_file(OUTPUT_GPKG_10, driver="GPKG")
            logging.info("Saved 10-minute isochores.")
        
    except Exception as e:
        logging.critical(f"An unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()