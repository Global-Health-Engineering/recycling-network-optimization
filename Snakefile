# Snakefile

# Load configuration
configfile: "config/config.yaml"

# Data paths
RAW_DATA = "/home/silas/rcp_project/rcp_project/data/raw_data"
DERIVED_DATA = "/home/silas/rcp_project/rcp_project/data/derived_data"
PLOTS_PATH = "/home/silas/rcp_project/rcp_project/data/plots"

# Set routing engine from config
ROUTING_ENGINE = config.get("routing_engine", "valhalla")  # Default to valhalla if not specified

# Update the all rule to include all outputs
rule all:
    input:
        DERIVED_DATA + "/workflow/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_opt.gpkg",
        PLOTS_PATH + "/workflow/method_comparison_table.tex",
        DERIVED_DATA + "/workflow/method_comparison.csv",

# Main workflow rules
include: "rules/data_preparation.smk"
include: "rules/optimisation_rules.smk"

# Analysis for number of facilities
include: "rules/p_analysis_rules.smk"

# Sensitivity analysis rules
include: "rules/sensitivity_rules.smk"