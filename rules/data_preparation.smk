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
    conda: "../envs/geo_env.yaml"
    script:
        "../scripts/population_allocation.py"

rule generate_isochrones:
    input:
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        iso_all=DERIVED_DATA + "/isochrones_all.gpkg",
        iso_merged=DERIVED_DATA + "/iso_merged.gpkg"
    log:
        "logs/isochores_calculations.log"
    conda: "../envs/geo_env.yaml"
    script:
        "../scripts/generate_isochrones.py"

rule calculate_distance_matrices:
    input:
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_locations=DERIVED_DATA + "/all_pot_sites.gpkg",
        demand_points=DERIVED_DATA + "/kmeans_clusters.gpkg"
    output:
        matrix_walking=DERIVED_DATA + "/distance_matrix_walking.csv"
    log:
        "logs/distance_matrix.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_matrix.py"

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
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_calc.py"

rule generate_potential_sites:
    input:
        slope_raster=DERIVED_DATA + "/slope_zurich.tif",
        trees=RAW_DATA + "/geodata_stadt_Zuerich/trees/data/data.gpkg",
        flats=DERIVED_DATA + "/flats_duration_current.gpkg",
        parking=RAW_DATA + "/osm_data/parking_lots_zurich.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        buildings=RAW_DATA + "/geodata_stadt_Zuerich/3d_buildings/data/data.gpkg",
        vbz=RAW_DATA + "/geodata_stadt_Zuerich/vbz/data/data.gpkg"
    output:
        sites=DERIVED_DATA + "/all_pot_sites.gpkg",
        map=PLOTS_PATH + "/suitable_sites_map.html"
    params:
        buffer_dist_vbz=2,         # meters buffer around VBZ infrastructure
        buffer_trees=2,            # meters buffer around trees
        max_slope=5,               # Maximum slope in degrees
        area_threshold=16,         # Minimum area in square meters
        buffer_buildings=14        # meters buffer around buildings
    log:
        "logs/site_selection/potential_sites.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/generate_potential_sites.py"

rule generate_demand_points:
    input:
        flats=DERIVED_DATA + "/flats_duration_current.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        gpkg=DERIVED_DATA + "/kmeans_clusters.gpkg",
        html_map=PLOTS_PATH + "/kmeans_clusters.html"
    params:
        n_clusters=1000  # Default value
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/demand_points.py"