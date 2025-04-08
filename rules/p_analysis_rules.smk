# Number of facilites (p-values) to test
P_VALUES = range(
    config["p_analysis"]["range"]["start"],
    config["p_analysis"]["range"]["end"],
    config["p_analysis"]["range"]["step"]
)

# Target rule for p-analysis
rule run_p_analysis:
    input:
        PLOTS_PATH + "/p-analysis/metrics_comparison_plot.png",
        PLOTS_PATH + "/p-analysis/optimality_gap_plot.png",
        DERIVED_DATA + "/p-analysis/summary_metrics.csv"

rule p_analysis_optimisation:
    input:
        demand_points=DERIVED_DATA + "/workflow/kmeans_clusters.gpkg", 
        potential_sites=DERIVED_DATA + "/workflow/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/workflow/distance_matrix.csv",
    output:
        sites=DERIVED_DATA + "/p-analysis/rcps_optimisation_{p}.gpkg",
        optimality_gap=DERIVED_DATA + "/p-analysis/optimality_gap_{p}.txt"
    params:
        num_facilities=lambda wildcards: int(wildcards.p),
        routing_engine=ROUTING_ENGINE
    log:
        "logs/p-analysis/linear_optimization_{p}.log"
    conda:
        "../envs/solver_env.yaml"
    script:
        "../scripts/linear_optimization.py"

rule p_analysis_distance_calculation:
    input:
        flats=DERIVED_DATA + "/workflow/flats_population.gpkg",
        rcps=DERIVED_DATA + "/p-analysis/rcps_optimisation_{p}.gpkg"
    output:
        duration=DERIVED_DATA + "/p-analysis/flats_duration_p_{p}.gpkg"
    params:
        routing_engine=ROUTING_ENGINE  # Pass the routing engine parameter
    log:
        "logs/p-analysis/distance_calc_{p}.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_calc_sensitivity.py"

rule analyse_p_results:
    input:
        duration_files=expand(DERIVED_DATA + "/p-analysis/flats_duration_p_{p}.gpkg", p=P_VALUES),
        optimality_gap_files=expand(DERIVED_DATA + "/p-analysis/optimality_gap_{p}.txt", p=P_VALUES)
    output:
        summary=DERIVED_DATA + "/p-analysis/summary_metrics.csv",
        metrics_plot=PLOTS_PATH + "/p-analysis/metrics_comparison_plot.png",
        optimality_plot=PLOTS_PATH + "/p-analysis/optimality_gap_plot.png"
    log:
        "logs/p-analysis/summary.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/analyse_p_results.py"