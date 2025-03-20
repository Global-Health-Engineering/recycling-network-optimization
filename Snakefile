# Snakefile

# Data paths
RAW_DATA = "/home/silas/rcp_project/rcp_project/data/raw_data"
DERIVED_DATA = "/home/silas/rcp_project/rcp_project/data/derived_data"
PLOTS_PATH = "/home/silas/rcp_project/rcp_project/data/plots"

# Update the all rule to include all outputs
rule all:
    input:
        DERIVED_DATA + "/workflow/isochrones_all.gpkg",
        DERIVED_DATA + "/workflow/iso_merged.gpkg",
        DERIVED_DATA + "/workflow/distance_matrix_walking.csv",
        DERIVED_DATA + "/workflow/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_opt.gpkg",
        DERIVED_DATA + "/workflow/merged_isochrones.gpkg",

# Main workflow rules
include: "rules/data_preparation.smk"
include: "rules/optimisation_rules.smk"

# Analysis for number of facilities
include: "rules/p_analysis_rules.smk"

# Sensitivity analysis rules
include: "rules/sensitivity_rules.smk"