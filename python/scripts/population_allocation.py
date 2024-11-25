import logging
import geopandas as gpd
import numpy as np

# Paths
POPULATION_POLYGON_PATH = '/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp'
OUTPUT_PATH = '/home/silas/projects/msc_thesis/data/derived_data/flats_w/_pop.gpkg'
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
        # -9999 values are nodata/0 in all colunns
        logging.info("Replacing -9999 values with 0.")
        population = population.replace(-999, 0)
        logging.info("Ensuring CRS match.")
        if flats.crs != population.crs:
            flats = flats.to_crs(population.crs)

        logging.info("Spatial joining flats with population polygons.")
        flats = gpd.sjoin(flats, population, how='left', predicate='within')
        logging.info(f"Spatial join completed, length of flats: {len(flats)}")
        flats.query('wstatlang == "Bestehend"', inplace=True)
        logging.info("Allocating population to flats based on 'wazim'.")
        est_pop = np.zeros(len(flats))

        for i, (_, row) in enumerate(flats.iterrows()):
            # cell_pop variable removed as it is not used
            # Calculate total population and total wazim in the current cell
            cell_pop = row['PERS_N']  
            cell_wazim = flats[flats['index_right'] == row['index_right']]['wazim'].sum()
            
            # Compute average number of people per wazim
            avg_people_per_wazim = cell_pop / cell_wazim if cell_wazim != 0 else 0
            
            # Allocate population based on average people per wazim
            est_pop[i] = row['wazim'] * avg_people_per_wazim

            #logging.info(f"Flat ID {row['egid']}: Allocated population {est_pop[idx]}")

        #global average population per wazim
        avg_pop_per_wazim = est_pop.sum() / flats['wazim'].sum()
        logging.info(f"Average population per wazim: {avg_pop_per_wazim}")
        # Replace NaN values with the average population per wazim
        est_pop[np.isnan(est_pop)] = avg_pop_per_wazim * flats['wazim'][np.isnan(est_pop)]

        flats['est_pop'] = est_pop
        logging.info("Population allocation completed.")
        flats.to_file(OUTPUT_PATH)
        logging.info(f"Allocated population saved to {OUTPUT_PATH}")
        logging.info(f"Total population allocated: {flats['est_pop'].sum()}")

    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    allocate_population(FLAT_DB_PATH, POPULATION_POLYGON_PATH)
