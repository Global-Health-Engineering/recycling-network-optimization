# Optimization and clustering rules

rule linear_optimisation:
    input:
        demand_points=DERIVED_DATA + "/kmeans_clusters.gpkg",
        potential_sites=DERIVED_DATA + "/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/distance_matrix_walking.csv",
        flats=DERIVED_DATA + "/flats_population.gpkg"
    conda:
        "envs/solver_env.yaml"
    output:
        sites=DERIVED_DATA + "/rcps_optimisation.gpkg"
    params:
        num_facilities=10,   # Number of NEW facilities to open
        pop_limit=1500       # Maximum population outside 10-minute radius
    log:
        "logs/linear_optimization.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/linear_optimization.py"

rule clustering_ors:
    input:
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_sites=DERIVED_DATA + "/all_pot_sites.gpkg"
    output:
        rcps_clustering_ors=DERIVED_DATA + "/rcps_clustering_ors.gpkg",
        flats_duration=DERIVED_DATA + "/flats_duration_current.gpkg",
        map_clustering_ors=PLOTS_PATH + "/map_clustering_ors.html"
    params:
        eps=0.005,           # DBSCAN epsilon parameter for clustering
        min_samples=20,      # DBSCAN minimum samples parameter
        iso_threshold=10     # Threshold for underserved buildings
    log:
        "logs/spatial_clustering.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/clustering_ors.py"

rule clustering_isochrones:
    input:
        isochrones=DERIVED_DATA + "/iso_merged.gpkg",
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_sites=DERIVED_DATA + "/all_pot_sites.gpkg"
    output:
        merged_isochrones=DERIVED_DATA + "/merged_isochrones.gpkg",
        clustered_sites=DERIVED_DATA + "/rcps_clustering_iso.gpkg",
        html_map=PLOTS_PATH + "/map_clustering.html"
    params:
        eps=0.005,           # DBSCAN epsilon parameter for clustering
        min_samples=20       # DBSCAN minimum samples parameter
    log:
        "logs/spatial_clustering.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/clustering_iso.py"