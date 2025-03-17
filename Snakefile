# Snakefile

# Data paths
RAW_DATA = "/home/silas/rcp_project/rcp_project/data/raw_data"
DERIVED_DATA = "/home/silas/rcp_project/rcp_project/data/derived_data"
PLOTS_PATH = "/home/silas/rcp_project/rcp_project/data/plots"

# Define the cluster numbers we want to generate
CLUSTERS = [10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500]

# Configuration
#configfile: "config.yaml"

# target rule for sensitivity analysis
rule run_sensitivity_analysis:
    input:
        # Generate sensitivity outputs for all cluster numbers
        expand(DERIVED_DATA + "/sensitivity_clusters/flats_duration_{n}.gpkg", n=CLUSTERS)

# Update the all rule to include all cluster variants
rule all:
    input:
        DERIVED_DATA + "/isochrones_all.gpkg",
        DERIVED_DATA + "/iso_merged.gpkg",
        expand(DERIVED_DATA + "sensitivity_clusters/kmeans_clusters_{n}.gpkg", n=CLUSTERS),
        expand(PLOTS_PATH + "sensitivity_clusters/html/kmeans_clusters_{n}.html", n=CLUSTERS),
        DERIVED_DATA + "/distance_matrix_trucks.csv",
        DERIVED_DATA + "/distance_matrix_walking.csv",
        DERIVED_DATA + "/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/flats_duration_opt.gpkg",
        DERIVED_DATA + "/merged_isochrones.gpkg",
        PLOTS_PATH + "/map_clustering_iso.html"

rule allocate_population:
    input:
        flats=RAW_DATA + "/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp",
        population=RAW_DATA + "/geodata_stadt_Zuerich/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp"
    params:
        flats_under_construction=True,
        flats_in_planning=True,
        exclusion_buffer=5
    output:
        DERIVED_DATA + "/flats_population.gpkg"
    log:
        "logs/population_allocation.log"
    conda: "envs/geo_env.yaml"
    script:
        "scripts/population_allocation.py"

rule generate_isochrones:
    input:
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        iso_5min=DERIVED_DATA + "/isochrones_5min.gpkg",
        iso_10min=DERIVED_DATA + "/isochrones_10min.gpkg",
        iso_all=DERIVED_DATA + "/isochrones_all.gpkg",
        iso_merged=DERIVED_DATA + "/iso_merged.gpkg"
    log:
        "logs/isochores_calculations.log"
    conda: "envs/geo_env.yaml"
    script:
        "scripts/generate_isochrones.py"

rule calculate_distance_matrices:
    input:
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_locations=DERIVED_DATA + "/all_pot_sites.gpkg",
        demand_points=DERIVED_DATA + "/kmeans_clusters.gpkg"
    output:
        matrix_trucks=DERIVED_DATA + "/distance_matrix_trucks.csv",
        matrix_walking=DERIVED_DATA + "/distance_matrix_walking.csv"
    log:
        "logs/distance_matrix.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/distance_matrix.py"

rule calculate_distances_to_rcp:
    input:
        flats = DERIVED_DATA + "/flats_population.gpkg",
        rcps1 = DERIVED_DATA+ "/rcps_clustering_ors.gpkg",
        rcps2 = DERIVED_DATA + "/rcps_clustering_iso.gpkg",
        rcps3 = DERIVED_DATA + "/rcps_optimisation.gpkg"
    output:
        DERIVED_DATA + "/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/flats_duration_opt.gpkg"
    log:
        "logs/distance_calc.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/distance_calc.py"

rule generate_demand_points:
    input:
        flats=DERIVED_DATA + "/flats_duration_current.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        gpkg=DERIVED_DATA + "sensitivity_clusters/kmeans_clusters_{n_clusters}.gpkg",
        html_map=PLOTS_PATH + "sensitivity_clusters/kmeans_clusters_{n_clusters}.html"
    params:
        n_clusters=lambda wildcards: int(wildcards.n_clusters)
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/demand_points.py"

rule linear_optimisation:
    input:
        demand_points=DERIVED_DATA + "/kmeans_clusters.gpkg",
        potential_sites=DERIVED_DATA + "/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/distance_matrix_walking.csv",
        flats=DERIVED_DATA + "/flats_population.gpkg"
    conda:
        "envs/solver_env.yaml"
    output:
        sites=DERIVED_DATA + "/rcps_optimisation.gpkg",
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
        min_samples=20,       # DBSCAN minimum samples parameter
        iso_threshold=10    # Threshold for underserved buildings
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

rule sensitivity_distance_matrices:
    input:
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_locations=DERIVED_DATA + "/all_pot_sites.gpkg",
        demand_points=DERIVED_DATA + "sensitivity_clusters/kmeans_clusters_{n_clusters}.gpkg"
    output:
        matrix_walking=DERIVED_DATA + "/sensitivity_clusters/distance_matrix_{n_clusters}.csv"
    log:
        "logs/sensitivity/distance_matrix_{n_clusters}.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/distance_matrix.py"
        
rule sensitivity_linear_optimisation:
    input:
        demand_points=DERIVED_DATA + "sensitivity_clusters/kmeans_clusters_{n_clusters}.gpkg",
        potential_sites=DERIVED_DATA + "/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/sensitivity_clusters/distance_matrix_{n_clusters}.csv",
        flats=DERIVED_DATA + "/flats_population.gpkg"
    output:
        sites=DERIVED_DATA + "/sensitivity_clusters/rcps_optimisation_{n_clusters}.gpkg"
    params:
        num_facilities=10,
        pop_limit=1500
    log:
        "logs/sensitivity/linear_optimization_{n_clusters}.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/linear_optimization.py"

rule sensitivity_distance_calculation:
    input:
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=DERIVED_DATA + "/sensitivity_clusters/rcps_optimisation_{n_clusters}.gpkg"
    output:
        duration=DERIVED_DATA + "/sensitivity_clusters/flats_duration_{n_clusters}.gpkg"
    log:
        "logs/sensitivity/distance_calc_{n_clusters}.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/distance_calc_sensitivity.py"   