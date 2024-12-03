# Snakefile

# Data paths
RAW_DATA = "/home/silas/projects/msc_thesis/data/raw_data"
DERIVED_DATA = "/home/silas/projects/msc_thesis/data/derived_data"
PLOTS_PATH = "/home/silas/projects/msc_thesis/docs/reports/thesis/figures"

# Configuration
#configfile: "config.yaml"

rule all:
    input:
        DERIVED_DATA + "/isochores_5min.gpkg",
        DERIVED_DATA + "/isochores_10min.gpkg",
        DERIVED_DATA + "/flats_subset_with_rcp.shp"

rule population_allocation:
    input:
        flats=RAW_DATA + "/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp",
        population=RAW_DATA + "/geodata_stadt_Zuerich/Raumliche_Bevolkerungsstatistik_-OGD/BEVOELKERUNG_HA_F.shp"
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
        iso_5min=DERIVED_DATA + "/isochores_5min.gpkg",
        iso_10min=DERIVED_DATA + "/isochores_10min.gpkg"
    log:
        "logs/isochores_calculations.log"
    conda: "envs/geo_env.yaml"
    script:
        "scripts/isochores_calculations.py"

rule calculate_distances_to_rcp:
    input:
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        DERIVED_DATA + "/flats_subset_with_rcp.shp"
    params:
        n=100 
    log:
        "logs/distance_calc.log"
    conda:
        "envs/geo_env.yaml"
    script:
        "scripts/distance_calc.py"