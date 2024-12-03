import geopandas as gpd
import numpy as np
from snakemake.logging import logger

# Remove hardcoded paths and configure logging using snakemake variables
flat_db_path = snakemake.input['flats']
population_polygon_path = snakemake.input['population']
output_path = snakemake.output[0]

def allocate_population(flat_db_path, population_polygon_path, output_path):
    try:
        logger.info("Loading flat database.")
        flats = gpd.read_file(flat_db_path)
        logger.info("Loading population polygons.")
        population = gpd.read_file(population_polygon_path)
        # -9999 values are nodata/0 in all columns
        logger.info("Replacing -9999 values with 0.")
        population = population.replace(-999, 0)
        logger.info("Ensuring CRS match.")
        if flats.crs != population.crs:
            flats = flats.to_crs(population.crs)

        logger.info("Spatial joining flats with population polygons.")
        flats = gpd.sjoin(flats, population, how='left', predicate='within')
        logger.info(f"Spatial join completed, length of flats: {len(flats)}")
        flats.query('wstatlang == "Bestehend"', inplace=True)
        logger.info("Allocating population to flats based on 'wazim'.")
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

            #logger.info(f"Flat ID {row['egid']}: Allocated population {est_pop[idx]}")

        # Global average population per wazim
        avg_pop_per_wazim = est_pop.sum() / flats['wazim'].sum()
        logger.info(f"Average population per wazim: {avg_pop_per_wazim}")
        # Replace NaN values with the average population per wazim
        est_pop[np.isnan(est_pop)] = avg_pop_per_wazim * flats['wazim'][np.isnan(est_pop)]

        flats['est_pop'] = est_pop
        logger.info("Population allocation completed.")
        flats.to_file(output_path)
        logger.info(f"Allocated population saved to {output_path}")
        logger.info(f"Total population allocated: {flats['est_pop'].sum()}")

    except Exception as e:
        logger.error(f"An error occurred: {e}")

if __name__ == "__main__":
    allocate_population(flat_db_path, population_polygon_path, output_path)
