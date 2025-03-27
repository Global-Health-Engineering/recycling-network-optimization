# Data preparation rules

rule allocate_population:
    input:
        flats=RAW_DATA + "/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp",
        population=RAW_DATA + "/geodata_stadt_Zuerich/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp"
    params:
        flats_under_construction=config["data_preparation"]["population_allocation"]["flats_under_construction"],
        flats_in_planning=config["data_preparation"]["population_allocation"]["flats_in_planning"],
        exclusion_buffer=config["data_preparation"]["population_allocation"]["exclusion_buffer"]
    output:
        DERIVED_DATA + "/workflow/flats_population.gpkg"
    log:
        "logs/population_allocation.log"
    conda: "../envs/geo_env.yaml"
    script:
        "../scripts/population_allocation.py"

rule generate_isochrones:
    input:
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        iso_all=DERIVED_DATA + "/workflow/isochrones_all.gpkg",
        iso_merged=DERIVED_DATA + "/workflow/iso_merged.gpkg"
    params:
        routing_engine=ROUTING_ENGINE  # Pass the routing engine parameter
    log:
        "logs/isochores_calculations.log"
    conda: "../envs/geo_env.yaml"
    script:
        "../scripts/generate_isochrones.py"

rule calculate_distance_matrices:
    input:
        potential_locations=DERIVED_DATA + "/workflow/all_pot_sites.gpkg",
        demand_points=DERIVED_DATA + "/workflow/kmeans_clusters.gpkg"
    output:
        matrix_walking=DERIVED_DATA + "/workflow/distance_matrix.csv"
    params:
        routing_engine=ROUTING_ENGINE
    log:
        "logs/distance_matrix.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_matrix.py"

rule calculate_distances_to_rcp:
    input:
        flats = DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps1 = DERIVED_DATA+ "/workflow/rcps_clustering_ors.gpkg",
        rcps2 = DERIVED_DATA + "/workflow/rcps_clustering_iso.gpkg",
        rcps3 = DERIVED_DATA + "/workflow/rcps_optimisation.gpkg"
    output:
        DERIVED_DATA + "/workflow/flats_duration_clustering_iso.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_clustering_ors.gpkg",
        DERIVED_DATA + "/workflow/flats_duration_opt.gpkg"
    params:
        routing_engine=ROUTING_ENGINE  # Pass the routing engine parameter
    log:
        "logs/distance_calc.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_calc.py"

rule generate_potential_sites:
    input:
        slope_raster=RAW_DATA + "/slope_zurich.tif",
        trees=RAW_DATA + "/geodata_stadt_Zuerich/trees/data/data.gpkg",
        parking_lots=RAW_DATA + "/osm_data/parking_lots_zurich.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        buildings=RAW_DATA + "/geodata_stadt_Zuerich/3d_buildings/data/data.gpkg",
        vbz=RAW_DATA + "/geodata_stadt_Zuerich/vbz/data/data.gpkg"
    output:
        potential_sites=DERIVED_DATA + "/workflow/all_pot_sites.gpkg",
        map=PLOTS_PATH + "/workflow/suitable_sites_map.html"
    params:
        buffer_dist_vbz=config["data_preparation"]["potential_sites"]["buffer_dist_vbz"],
        buffer_trees=config["data_preparation"]["potential_sites"]["buffer_trees"],
        max_slope=config["data_preparation"]["potential_sites"]["max_slope"],
        area_threshold=config["data_preparation"]["potential_sites"]["area_threshold"],
        buffer_buildings=config["data_preparation"]["potential_sites"]["buffer_buildings"]
    log:
        "logs/generate_potential_sites.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/suitability_analysis_nb.py"

rule generate_demand_points:
    input:
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        gpkg=DERIVED_DATA + "/workflow/kmeans_clusters.gpkg",
        html_map=PLOTS_PATH + "/workflow/kmeans_clusters.html"
    params:
        n_clusters=config["data_preparation"]["demand_points"]["n_clusters"]
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/demand_points.py"