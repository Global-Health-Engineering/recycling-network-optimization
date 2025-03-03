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
        population = population.replace(-9999, 0)
        logger.info("Ensuring CRS match.")
        if flats.crs != population.crs:
            flats = flats.to_crs(population.crs)

        logger.info("Spatial joining flats with population polygons.")
        flats = gpd.sjoin(flats, population, how='left', predicate='within')
        logger.info(f"Spatial join completed, length of flats: {len(flats)}")

        # Categories which should always be included
        status_categories = ["Bestehend", "Freigegeben"]

        # Add additional categories based on parameters
        if snakemake.params.get('flats_under_construction', False):
            status_categories.append("Im Bau")
            
        if snakemake.params.get('flats_in_planning', False):
            status_categories.append("Bewilligt")
 
        # Create the filter query
        status_filter = " or ".join([f'wstatlang == "{status}"' for status in status_categories])
        flats.query(status_filter, inplace=True)

        logger.info("Filtered flats based on construction status, \n allocating population to flats based on 'wazim'.") 

        est_pop = np.zeros(len(flats))

        # Precompute the cell-level average people per wazim,
        # using the flats with status "Bestehend" when available.
        cell_avg = {}
        for cell, group in flats.groupby("index_right"):
            bestaend = group[group["wstatlang"] == "Bestehend"]
            if not bestaend.empty:
                total_wazim = bestaend["wazim"].sum()
                # Assuming the cell population (PERS_N) is the same for every flat in the cell
                cell_pop = bestaend.iloc[0]["PERS_N"]
                cell_avg[cell] = cell_pop / total_wazim if total_wazim != 0 else 0
            else:
                # For cells with no "Bestehend" flat, fall back on overall cell average
                total_wazim = group["wazim"].sum(ata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp, /home/silas/rcp_project/rcp_project/data/derived_data/rcps_clustering_ors.gpkg, /home/silas/rcp_project/rcp_project/data/derived_data/rcps_clustering_iso.gpkg
    output: /home/silas/rcp_project/rcp_project/data/derived_data/flats_duration_current.gpkg, /home/silas/rcp_project/rcp_project/data/derived_data/flats_duration_clustering_iso.gpkg, /home/silas/rcp_project/rcp_project/data/derived_data/flats_duration_clustering_ors.gpkg
    log: logs/distance_calc.log (check log file(s) for error details)
    conda-env: /home/silas/rcp_project/rcp_project/.snakemake/conda/d7f122aa3ec4f1aca2098dd7ccd0ea25_

Complete log: .snakemake/log/2025-03-03T155327.625258.snakemake.log
WorkflowError:
At least one job did not complete suc)
                cell_pop = group.iloc[0]["PERS_N"]
                cell_avg[cell] = cell_pop / total_wazim if total_wazim != 0 else 0

        # Allocate population per flat based on status
        for i, (_, row) in enumerate(flats.iterrows()):
            cell_average = cell_avg.get(row["index_right"], 0)
            if row["wstatlang"] == "Bestehend":
                # For 'Bestehend' flats, compute the cell average from only 'Bestehend' flats
                flats_in_cell = flats[(flats["index_right"] == row["index_right"]) &
                                      (flats["wstatlang"] == "Bestehend")]
                total_wazim = flats_in_cell["wazim"].sum()
                avg_people_per_wazim = row["PERS_N"] / total_wazim if total_wazim != 0 else 0
                est_pop[i] = row["wazim"] * avg_people_per_wazim
            else:
                # For all other statuses, use the precomputed cell average
                est_pop[i] = row["wazim"] * cell_average

        # Global average population per wazim
        avg_pop_per_wazim = est_pop.sum() / flats['wazim'].sum()
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
