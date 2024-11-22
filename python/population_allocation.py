import logging
import geopandas as gpd
import numpy as np

# Paths
POPULATION_POLYGON_PATH = '/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp'
OUTPUT_PATH = '/home/silas/projects/msc_thesis/data/derived_data/population_allocation.gpkg'
FLAT_DB_PATH = "/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp"

# Configure logging
logging.basicConfig(
    filename='population_allocation.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def allocate_population(flat_db_path, population_polygon_path):
    try:
        logging.info("Loading flat database.")
        flats = gpd.read_file(flat_db_path)

        logging.info("Loading population polygons.")
        population = gpd.read_file(population_polygon_path)

        logging.info("Ensuring CRS match.")
        if flats.crs != population.crs:
            flats = flats.to_crs(population.crs)

        logging.info("Spatial joining flats with population polygons.")
        flats = gpd.sjoin(flats, population, how='left', op='within')

        logging.info("Allocating population to flats based on 'wazim'.")
        # grouped variable removed as it is not used

        est_pop = np.zeros(len(flats))

        for idx, row in flats.iterrows():
            # cell_pop variable removed as it is not used
            wazim = row['wazim']
            est_pop[idx] = wazim  # Allocate population based on number of rooms in the flats

            #logging.info(f"Flat ID {row['egid']}: Allocated population {est_pop[idx]}")

        flats['est_pop'] = est_pop
        logging.info("Population allocation completed.")
        flats.to_file(OUTPUT_PATH)
        logging.info(f"Allocated population saved to {OUTPUT_PATH}")

    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    allocate_population(FLAT_DB_PATH, POPULATION_POLYGON_PATH)
