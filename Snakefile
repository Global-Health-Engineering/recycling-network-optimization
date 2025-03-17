# Snakefile

# Data paths
RAW_DATA = "/home/silas/rcp_project/rcp_project/data/raw_data"
DERIVED_DATA = "/home/silas/rcp_project/rcp_project/data/derived_data"
PLOTS_PATH = "/home/silas/rcp_project/rcp_project/data/plots"

# Define the cluster numbers we want to generate - keep in main file as it's used in multiple places
CLUSTERS = [10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500]

# Update the all rule to include all outputs
rule all:
    input:
        DERIVED_DATA + "/isochrones_all.gpkg",
        DERIVED_DATA + "/iso_merged.gpkg",
        DERIVED_DATA + "/distance_matrix_walking.csv",
        DERIVED_DATA + "/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/flats_duration_opt.gpkg",
        DERIVED_DATA + "/merged_isochrones.gpkg",
        PLOTS_PATH + "/map_clustering_iso.html"

# Main workflow rules
include: "rules/data_preparation.smk"
include: "rules/optimization_rules.smk"

# Sensitivity analysis rules
include: "rules/sensitivity_rules.smk"