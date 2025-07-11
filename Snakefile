# Snakefile

# Load configuration
configfile: "config/config.yaml"

# Data paths
RAW_DATA = "data/raw_data"
DERIVED_DATA = "data/derived_data"
PLOTS_PATH = "data/plots"

# Set routing engine from config
ROUTING_ENGINE = config.get("routing_engine", "valhalla")  # Default to valhalla if not specified

# Update the all rule to include all outputs
rule all:
    input:
        DERIVED_DATA + "/workflow/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_opt.gpkg",
        PLOTS_PATH + "/workflow/method_comparison_table.tex",
        PLOTS_PATH + "/workflow/map_clustering_iso.html",
        DERIVED_DATA + "/workflow/method_comparison.csv",

# Main workflow rules
include: "rules/data_preparation.smk"
include: "rules/optimisation_rules.smk"

# Analysis for number of facilities
include: "rules/p_analysis_rules.smk"

# Sensitivity analysis rules
include: "rules/sensitivity_rules.smk"