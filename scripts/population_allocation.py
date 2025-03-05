import geopandas as gpd
import numpy as np
from snakemake.logging import logger
import pandas as pd

# Remove hardcoded paths and configure logging using snakemake variables
flat_db_path = snakemake.input['flats']
population_polygon_path = snakemake.input['population']
output_path = snakemake.output[0]

# buffer around planned flats to exclude existing flats which will be demolished
exclusion_buffer = snakemake.params.get('exclusion_buffer') 

def allocate_population(flat_db_path, population_polygon_path, output_path):
    try:
        logger.info("Loading flat database.")
        flats = gpd.read_file(flat_db_path)

        # Group flats by building to speed up spatial join
        flats=flats.groupby('egid').agg({'geometry': 'first', 'wazim': 'sum', 'wstatlang':'first'}).reset_index()
        flats = gpd.GeoDataFrame(flats, crs='EPSG:2056', geometry=flats['geometry'])
        population = gpd.read_file(population_polygon_path)

        # -9999 values are nodata/0 in all columns
        population = population.replace(-999, 0)

        # Ensure CRS match
        if flats.crs != population.crs:
            flats = flats.to_crs(population.crs)

        flats = gpd.sjoin(flats, population, how='left', predicate='intersects')
        logger.info(f"Spatial join completed, length of flats: {len(flats)}")

        # Categories which should always be included
        status_categories = ["Bestehend", "Freigegeben"]

        # Add additional categories based on parameters
        if snakemake.params.get('flats_under_construction', True):
            status_categories.append("Im Bau")
            
        if snakemake.params.get('flats_in_planning', True):
            status_categories.append("Bewilligt")
 
        # Create the filter query
        status_filter = " or ".join([f'wstatlang == "{status}"' for status in status_categories])
        flats.query(status_filter, inplace=True)

        # allocate population to flats based on wazim
        if snakemake.params.get('flats_in_planning', False):
            planned_flats = flats[flats["wstatlang"] == "Bewilligt"]
            if not planned_flats.empty:
                planned_buffer = planned_flats.copy()
                planned_buffer["geometry"] = planned_buffer.geometry.buffer(exclusion_buffer)
                existing_flats = flats[flats["wstatlang"] == "Bestehend"]
                join_df = gpd.sjoin(existing_flats, planned_buffer[["geometry"]], how="left", predicate="within", 
                    lsuffix="_left", rsuffix="_planned")
                # Find the correct index column (might be 'index_planned' now)
                index_cols = [col for col in join_df.columns if 'index_planned' in col.lower()]
                if index_cols:
                    right_index_col = index_cols[0]
                    indices_to_remove = join_df[join_df[right_index_col].notnull()].index
                    flats = flats.drop(indices_to_remove)
                else:
                    logger.warning("No suitable index column found in spatial join result. Skipping removal step.")

        logger.info("Filtered flats based on construction status, allocating population to flats based on 'wazim'.") 
        est_pop = np.zeros(len(flats))

        # Create Series for population instead of numpy array so indices are preserved
        est_pop = pd.Series(np.zeros(len(flats)), index=flats.index)

        # Step 1: Allocate population for 'Bestehend' flats using cell-specific averages
        existing_mask = flats["wstatlang"] == "Bestehend"
        for i, row in flats[existing_mask].iterrows():
            flats_in_cell = flats[(flats["index_right"] == row["index_right"]) & existing_mask]
            total_wazim = flats_in_cell["wazim"].sum()
            avg_people_per_wazim = row["PERS_N"] / total_wazim if total_wazim != 0 else 0
            est_pop[i] = row["wazim"] * avg_people_per_wazim

        # Step 2: Compute global average based on 'Bestehend' flats only
        total_population_existing = est_pop[existing_mask].sum()
        total_wazim_existing = flats.loc[existing_mask, "wazim"].sum()
        avg_pop_per_wazim = total_population_existing / total_wazim_existing if total_wazim_existing != 0 else 0

        # Allocate population for the other flats using the global average
        non_existing_mask = flats["wstatlang"] != "Bestehend"
        for i, row in flats[non_existing_mask].iterrows():
            est_pop[i] = row["wazim"] * avg_pop_per_wazim
        logger.info(f"Average population per wazim: {avg_pop_per_wazim}")

        # Replace NaN and 0 values with the average population per wazim
        est_pop[(np.isnan(est_pop)) | (est_pop == 0)] = avg_pop_per_wazim * flats['wazim'][(np.isnan(est_pop)) | (est_pop == 0)]

        flats['est_pop'] = est_pop
        logger.info("Population allocation completed.")
        # select only relevant columns
        flats = flats[['egid', 'est_pop', 'geometry']]
        flats.to_file(output_path)
        logger.info(f"Allocated population saved to {output_path}")
        logger.info(f"Total population allocated: {flats['est_pop'].sum()}")

    except Exception as e:
        logger.error(f"An error occurred: {e}")

if __name__ == "__main__":
    allocate_population(flat_db_path, population_polygon_path, output_path)
