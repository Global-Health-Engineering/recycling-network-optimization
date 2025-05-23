# Derived Data

This directory contains processed and intermediate data produced by the analysis workflow.

## Structure

- `workflow/` - Contains the main outputs from each stage of the analysis
  - `flats_population.gpkg` - Buildings with estimated population counts
  - `iso_merged.gpkg` - Merged walking time isochrones
  - `all_pot_sites.gpkg` - Suitable sites for new RCPs
  - `kmeans_clusters.gpkg` - Demand point clusters
  - `distance_matrix.csv` - Walking duration calculations
  - `rcps_clustering_ors.gpkg` - Routing-based clustering results
  - `rcps_clustering_iso.gpkg` - Isochrone-based clustering results
  - `rcps_optimisation.gpkg` - Linear optimization results
  - `flats_duration_current.gpkg` - Current walking durations
  - `flats_duration_opt.gpkg` - Optimized walking durations
  - `method_comparison.csv` - Comparison of optimization methods

- `sensitivity_analysis/` - Contains outputs from sensitivity analysis
  - `kmeans_clusters_{n}.gpkg` - Demand point clusters with varying numbers
  - `distance_matrix_{n}.csv` - Walking duration matrices with different cluster counts
  - `rcps_optimisation_{n}.gpkg` - Optimization results with varying parameters
  - `optimality_gap_{n}.txt` - Optimization quality metrics

- `p_analysis/` - Contains outputs from the p-analysis (facility count analysis)
  - Results from varying the number of facilities (p-value)
  - Analysis of diminishing returns in coverage improvement

## Note

All files in this directory (except README.md) are ignored by git and will be generated by running the Snakemake workflow.
