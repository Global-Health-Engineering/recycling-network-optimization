


# Number of facilites to test, 1-15
P_VALUES= range(0,15, 1)

# Target rule for p-analysis

rule run_p_analysis:
    input:
        PLOTS_PATH + "/p-analysis/p_comparison_plot.png",
        DERIVED_DATA + "/p-analysis/summary_metrics.csv",

rule p_analysis_optimisation:
    input:
        demand_points=DERIVED_DATA + 
        potential_sites=DERIVED_DATA + "/all_pot_sites.gpkg",
        distance_matrix=DERIVED_DATA + "/sensitivity_clusters/distance_matrix_{n_clusters}.csv",
        flats=DERIVED_DATA + "/flats_population.gpkg"
    output:
        sites=DERIVED_DATA + "//rcps_optimisation_{p}.gpkg"
    params:
        num_facilities= lambda wildcards: int(wildcards.p),
        pop_limit=3000
    log:
        "logs/sensitivity/linear_optimization_{p}.log"
    conda:
        "../envs/solver_env.yaml"
    script:
        "../scripts/linear_optimization.py"

rule p_analysis_distance_calculation:
    input:
        flats=DERIVED_DATA + "/flats_population.gpkg",
        rcps=DERIVED_DATA + "/sensitivity_clusters/rcps_optimisation_{p}.gpkg"
    output:
        duration=DERIVED_DATA + "/p-analysis/flats_duration_p_{p}.gpkg"
    log:
        "logs/sensitivity/distance_matrix_{p}.log"
    conda:
        "../envs/geo_env.yaml"
    script:
        "../scripts/distance_calc_sensitivity.py"

rule analyse_p_results:
    input:
        duration_files=expand(DERIVED_DATA + "/p-analysis/flats_duration_p_{p}.gpkg", p=P_VALUES)
    output:
        summary=DERIVED_DATA + "/p-analysis/summary_metrics.csv",
        plot=PLOTS_PATH + "/p-analysis/p_comparison_plot.png"
    log:
        "logs/p-analysis/summary.log"
    conda:
        "../envs/geo_env.yaml"