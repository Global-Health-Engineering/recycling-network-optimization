# Data preparation rules

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