# Sensitivity analysis rules
CLUSTERS = config["sensitivity_analysis"]["clusters"]

# Target rule for sensitivity analysis
rule run_sensitivity_analysis:
    input:
        PLOTS_PATH + "/sensitivity_analysis/comparison_plot.png",
        DERIVED_DATA + "/sensitivity_analysis/summary_metrics.csv",
        expand(DERIVED_DATA + "/sensitivity_analysis/flats_duration_{n}.gpkg", n=CLUSTERS)


# Generate demand points for sensitivity analysis
rule sensitivity_n_demand_points:
    input:
        flats=DERIVED_DATA + "/workflow/flats_duration_current.gpkg",
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp"
    output:
        gpkg=DERIVED_DATA + "/sensitivity_analysis/kmeans_clusters_{n_clusters}.gpkg",
        html_map=PLOTS_PATH + "/sensitivity_analysis/html/kmeans_clusters_{n_clusters}.html"
    params:
        n_clusters=lambda wildcards: int(wildcards.n_clusters)
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/demand_points.py"

rule sensitivity_distance_matrices:
    input:
        rcps=RAW_DATA + "/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp",
        potential_locations=DERIVED_DATA + "/workflow/all_pot_sites.gpkg",
        demand_points=DERIVED_DATA + "/sensitivity_analysis/kmeans_clusters_{n_clusters}.gpkg"
    output:
        matrix_walking=DERIVED_DATA + "/sensitivity_analysis/distance_matrix_{n_clusters}.csv"
    log:
        "logs/sensitivity/distance_matrix_{n_clusters}.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_matrix.py"
        
rule sensitivity_linear_optimisation:
    input:
        demand_points=DERIVED_DATA + "/sensitivity_analysis/kmeans_clusters_{n_clusters}.gpkg",
        potential_sites=DERIVED_DATA + "/workflow/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/sensitivity_analysis/distance_matrix_{n_clusters}.csv",
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg"
    output:
        sites=DERIVED_DATA + "/sensitivity_analysis/rcps_optimisation_{n_clusters}.gpkg"
    params:
        num_facilities=config["sensitivity_analysis"]["linear_optimisation"]["num_facilities"],
    log:
        "logs/sensitivity/linear_optimization_{n_clusters}.log"
    conda:
        "../envs/solver_env.yaml"
    script:
        "../scripts/linear_optimization.py"

rule sensitivity_distance_calculation:
    input:
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps=DERIVED_DATA + "/sensitivity_analysis/rcps_optimisation_{n_clusters}.gpkg"
    output:
        duration=DERIVED_DATA + "/sensitivity_analysis/flats_duration_{n_clusters}.gpkg"
    log:
        "logs/sensitivity/distance_calc_{n_clusters}.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_calc_sensitivity.py"

# Analysis rule for sensitivity results
rule analyze_sensitivity_results:
    input:
        duration_files=expand(DERIVED_DATA + "/sensitivity_analysis/flats_duration_{n}.gpkg", n=CLUSTERS)
    output:
        summary=DERIVED_DATA + "/sensitivity_analysis/summary_metrics.csv",
        plot=PLOTS_PATH + "/sensitivity_analysis/comparison_plot.png"
    log:
        "logs/sensitivity/analysis.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/sensitivity_analysis.py"