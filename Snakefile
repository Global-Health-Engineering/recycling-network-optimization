# Snakefile

# Data paths
RAW_DATA = "/home/silas/rcp_project/rcp_project/data/raw_data"
DERIVED_DATA = "/home/silas/rcp_project/rcp_project/data/derived_data"
PLOTS_PATH = "/home/silas/rcp_project/rcp_project/data/plots"

# Configuration
#configfile: "config.yaml"

rule all:
    input:
        DERIVED_DATA + "/isochrones_5min.gpkg",
        DERIVED_DATA + "/isochrones_10min.gpkg",
        DERIVED_DATA + "/isochrones_all.gpkg",
        DERIVED_DATA + "/iso_merged.gpkg",
        DERIVED_DATA + "/kmeans_clusters.gpkg",
        PLOTS_PATH + "/kmeans_clusters.html",
        DERIVED_DATA + "/distance_matrix_trucks.csv",
        DERIVED_DATA + "/distance_matrix_walking.csv",
        DERIVED_DATA + "/flats_duration_current.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_ors.gpkg"
        #DERIVED_DATA + "/flats_duration_optimisation_1.gpkg"

rule allocate_population:
    input:
        flats=RAW_DATA + "/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp",
        population=RAW_DATA + "/geodata_stadt_Zuerich/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp"
    params:
        flats_under_construction=True,
        flats_in_planning=True
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
        rcps1 = RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        rcps2 = DERIVED_DATA+ "/rcps_clustering_ors.gpkg",
        rcps3 = DERIVED_DATA + "/rcps_clustering_iso.gpkg" 
        #rcps4 = DERIVED_DATA + "/rcps_optimisation_1.gpkg"
    output:
        DERIVED_DATA + "/flats_duration_current.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/flats_duration_clustering_ors.gpkg"
        #DERIVED_DATA + "/flats_duration_optimisation_1.gpkg"
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
        gpkg=DERIVED_DATA + "/kmeans_clusters.gpkg",
        html_map=PLOTS_PATH + "/kmeans_clusters.html"
    params:
        n_clusters=3000
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/demand_points.py"
        #DERIVED_DATA + "/flats_duration_optimisation_1.gpkg"