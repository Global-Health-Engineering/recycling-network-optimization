# Euclidean distance analysis rules
# This file contains rules for calculating Euclidean walking durations as an alternative to routing engine calculations

rule euclidean_distances_to_rcp:
    input:
        flats = DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps_current = RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        rcps1 = DERIVED_DATA + "/workflow/rcps_clustering_ors.gpkg",
        rcps2 = DERIVED_DATA + "/workflow/rcps_clustering_iso.gpkg",
        rcps3 = DERIVED_DATA + "/workflow/rcps_optimisation.gpkg"
    output:
        DERIVED_DATA + "/euclidean_analysis/flats_duration_current_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/flats_duration_clustering_iso_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/flats_duration_clustering_ors_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/flats_duration_opt_euclidean.gpkg"
    params:
        walking_speed_kmh=config.get("euclidean_analysis", {}).get("walking_speed_kmh", 5.0)  # Default 5 km/h walking speed
    log:
        "logs/euclidean_distance_calc.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/euclidean_distance_calc.py"

# Rule to run all euclidean analysis
rule euclidean_analysis_all:
    input:
        DERIVED_DATA + "/euclidean_analysis/flats_duration_current_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/flats_duration_clustering_iso_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/flats_duration_clustering_ors_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/flats_duration_opt_euclidean.gpkg",
        DERIVED_DATA + "/euclidean_analysis/euclidean_routing_comparison.tex"

# Rule to generate comprehensive comparison table
rule euclidean_routing_comparison:
    input:
        # Euclidean results
        current_euclidean = DERIVED_DATA + "/euclidean_analysis/flats_duration_current_euclidean.gpkg",
        iso_euclidean = DERIVED_DATA + "/euclidean_analysis/flats_duration_clustering_iso_euclidean.gpkg",
        ors_euclidean = DERIVED_DATA + "/euclidean_analysis/flats_duration_clustering_ors_euclidean.gpkg",
        opt_euclidean = DERIVED_DATA + "/euclidean_analysis/flats_duration_opt_euclidean.gpkg",
        # Routing results (from main workflow)
        current_routing = DERIVED_DATA + "/workflow/flats_duration_current.gpkg",
        iso_routing = DERIVED_DATA + "/workflow/flats_duration_clustering_iso.gpkg",
        ors_routing = DERIVED_DATA + "/workflow/flats_duration_clustering_ors.gpkg",
        opt_routing = DERIVED_DATA + "/workflow/flats_duration_opt.gpkg"
    output:
        DERIVED_DATA + "/euclidean_analysis/euclidean_routing_comparison.tex"
    params:
        data_dir = DERIVED_DATA
    log:
        "logs/euclidean_routing_comparison.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/generate_euclidean_comparison.py"
