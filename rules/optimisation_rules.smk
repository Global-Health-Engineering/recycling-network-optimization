# Optimization and clustering rules

rule linear_optimisation:
    input:
        demand_points=DERIVED_DATA + "/workflow/kmeans_clusters.gpkg",
        potential_sites=DERIVED_DATA + "/workflow/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/workflow/distance_matrix.csv",
    conda:
        "../envs/solver_env.yaml"
    output:
        sites=DERIVED_DATA + "/workflow/rcps_optimisation.gpkg"
    params:
        num_facilities=config["optimization"]["linear_optimisation"]["num_facilities"],
        routing_engine=ROUTING_ENGINE  # Pass the routing engine parameter
    log:
        "logs/linear_optimization.log"
    script:
        "../scripts/linear_optimization.py"

rule clustering_ors:
    input:
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_sites=DERIVED_DATA + "/workflow/all_pot_sites.gpkg"
    output:
        rcps_clustering_ors=DERIVED_DATA + "/workflow/rcps_clustering_ors.gpkg",
        flats_duration=DERIVED_DATA + "/workflow/flats_duration_current.gpkg",
        map_clustering_ors=PLOTS_PATH + "/workflow/map_clustering_ors.html"
    params:
        eps=config["optimization"]["clustering"]["dbscan"]["eps"],
        min_samples=config["optimization"]["clustering"]["dbscan"]["min_samples"],
        iso_threshold=config["optimization"]["clustering"]["clustering_ors"]["iso_threshold"],
        routing_engine=ROUTING_ENGINE
    log:
        "logs/spatial_clustering.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/clustering_ors.py"

rule clustering_isochrones:
    input:
        isochrones=DERIVED_DATA + "/workflow/iso_merged.gpkg",
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_sites=DERIVED_DATA + "/workflow/all_pot_sites.gpkg"
    output:
        merged_isochrones=DERIVED_DATA + "/workflow/merged_isochrones.gpkg",
        clustered_sites=DERIVED_DATA + "/workflow/rcps_clustering_iso.gpkg",
        html_map=PLOTS_PATH + "/workflow/map_clustering.html"
    params:
        eps=config["optimization"]["clustering"]["dbscan"]["eps"],
        min_samples=config["optimization"]["clustering"]["dbscan"]["min_samples"],
        routing_engine=ROUTING_ENGINE
    log:
        "logs/spatial_clustering.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/clustering_iso.py"

rule compare_methods:
    input:
        duration_current=DERIVED_DATA + "/workflow/flats_duration_current.gpkg",
        duration_iso=DERIVED_DATA + "/workflow/flats_duration_clustering_iso.gpkg",
        duration_ors=DERIVED_DATA + "/workflow/flats_duration_clustering_ors.gpkg",
        duration_opt=DERIVED_DATA + "/workflow/flats_duration_opt.gpkg",
        rcp_opt=DERIVED_DATA + "/workflow/rcps_optimisation.gpkg",
        rcp_iso=DERIVED_DATA + "/workflow/rcps_clustering_iso.gpkg",
        rcp_ors=DERIVED_DATA + "/workflow/rcps_clustering_ors.gpkg",
        rcp_current=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        comparison_csv=DERIVED_DATA + "/workflow/method_comparison.csv",
        comparison_plot=PLOTS_PATH + "/workflow/method_comparison.png",
        efficiency_plot=PLOTS_PATH + "/workflow/people_brought_in_per_rcp.png",
        latex_table=PLOTS_PATH + "/workflow/method_comparison_table.tex"
    log:
        "logs/method_comparison.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/method_analysis.py"