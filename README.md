<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15320179.svg)](https://doi.org/10.5281/zenodo.15320179)
<!-- badges: end -->

<h1> Algorithms and data for optimizing the recycling network in the city of Zurich, Switzerland </h1>

<b>Contributors</b>  
- Silas Schweizer <a href="https://orcid.org/0009-0000-1284-7063">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /> 0009-0000-1284-7063
</a> *author, developer, maintainer*  
- Jakub Tkaczuk <a href="https://orcid.org/0000-0001-7997-9423">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /> 0000-0001-7997-9423
</a> *supervisor, maintainer*  
- Elizabeth Tilley <a href="https://orcid.org/0000-0002-2095-9724">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /> 0000-0002-2095-9724
</a> *supervisor*  
- Nicolas Seemann-Ricard <a href="https://orcid.org/0000-0002-0945-7475">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /> 0000-0002-0945-7475
</a> *technical support*  

<br>
<p align="middle"> 
<img src="img/ETH_GHE_logo_negative.svg" width=600>
<br><br>
This repository compliments the openly-accessible master’s thesis, available on the<br \>  
<a href="">ETH Research Collection</a>.
</p>

## Project Overview
This project aims to optimize the placement of Recycling Collection Points (RCPs) in Zurich, Switzerland. Using spatial analysis and optimization techniques, we identify optimal locations for new RCPs to improve accessibility for residents while minimizing the total number of facilities needed. The workflow uses a routing engine to calculate actual walking durations using the street network rather than relying on airline distances.


## Directory Structure
```
rcp_project/
├── config/                # Configuration files
│   └── config.yaml        # Main config with analysis parameters
├── data/                  # Data directory
│   ├── derived_data/      # Processed datasets
│   │   └── workflow/      # Intermediate and final analysis outputs
│   ├── plots/             # Generated visualizations
│   └── raw_data/          # Original unprocessed data
│       └── geodata_stadt_Zuerich/          # Geographic data from Zurich
│           ├── ssz.gwr_stzh_wohnungen.shp  # Building statistics
│           ├── BEVOELKERUNG_HA_F.shp       # Population data
│           ├── Gemeindegrenzen_-OGD.gpkg   # Municipal boundaries
│           └── stzh.poi_sammelstelle_view.shp # Existing RCPs
├── envs/                  # Conda environment definitions
│   ├── geo_env.yaml       # Environment for geospatial analysis
│   ├── snake_env.yaml     # Environment for running Snakemake
│   └── solver_env.yaml    # Environment for optimization solvers
├── img/                   # Images for documentation
├── logs/                  # Logs from workflow execution
├── rules/                 # Snakemake workflow rules
│   ├── data_preparation.smk  # Rules for data preparation
│   ├── optimisation_rules.smk # Rules for optimization
│   ├── p_analysis_rules.smk   # Rules for p-analysis analysis
│   └── sensitivity_rules.smk # Rules for sensitivity analysis
├── scripts/               # Python scripts used in workflow
└── Snakefile              # Main workflow definition
```

## Data Description

### Raw Data
- **Building Statistics**: Detailed information about buildings in Zurich (`ssz.gwr_stzh_wohnungen.shp`)
- **Population Data**: Spatial distribution of population in hectare cells (`BEVOELKERUNG_HA_F.shp`)
- **Municipal Boundaries**: Administrative boundaries of Zurich (`Gemeindegrenzen_-OGD.gpkg`)
- **Existing RCPs**: Current recycling collection points (`stzh.poi_sammelstelle_view.shp`)

### Derived Data
- **Population-Allocated Buildings**: Buildings with estimated population counts (`flats_population.gpkg`)
- **Suitable RCP Sites**: Potential locations for new RCPs (`all_pot_sites.gpkg`)
- **Isochrones**: Walking time service areas around RCPs (`merged_isochrones.gpkg`, `isochrones_all.gpkg`)
- **Walking Duration**: Buildings with calculated walking time to nearest RCP (`flats_duration_current.gpkg`)
- **Clustering Results**: Results from clustering analyses (`rcps_clustering_ors.gpkg`, `rcps_clustering_iso.gpkg`)
- **Optimization Results**: Solutions from optimization models (`rcps_optimisation.gpkg`)

## Comprehensive Workflow

The project workflow consists of the following main stages:

### 1. Data Preparation

1. **Population Allocation**
   - Allocate population to residential buildings based on building size and type
   - Include new buildings (under construction and approved) for future planning
   ```bash
   snakemake --use-conda data/derived_data/workflow/flats_population.gpkg
   ```

2. **Isochrone Generation**
   - Generate walking time service areas (isochrones) around existing RCPs
   - Merge isochrones to identify areas with adequate service coverage
   ```bash
   snakemake --use-conda data/derived_data/workflow/iso_merged.gpkg
   ```

3. **Suitability Analysis**
   - Identify suitable sites for new RCPs based on land use, accessibility, and distance to existing RCPs
   - Apply constraints (no placement in certain zones, minimum distance from buildings)
   ```bash
   snakemake --use-conda data/derived_data/workflow/all_pot_sites.gpkg
   ```

4. **Demand Point Generation**
   - Generate representative demand points using k-means clustering on building locations
   - Weight clusters by population to better represent demand distribution
   ```bash
   snakemake --use-conda data/derived_data/workflow/kmeans_clusters.gpkg
   ```

5. **Distance Matrix Calculation**
   - Calculate walking durations between potential RCP sites and demand points
   - Uses routing services (Valhalla or OpenRouteService) for accurate walking times
   ```bash
   snakemake --use-conda data/derived_data/workflow/distance_matrix.csv
   ```

### 2. Optimization Methods

1. **Clustering-Based Optimization with Routing**
   - Identify underserved buildings (>10 min walking time to nearest RCP)
   - Cluster underserved buildings to identify potential new RCP locations
   - Select optimal placement of new RCPs within each cluster
   ```bash
   snakemake --use-conda data/derived_data/workflow/rcps_clustering_ors.gpkg
   ```

2. **Isochrone-Based Optimization**
   - Identify areas outside existing service coverage (10-min walking isochrones)
   - Find clusters of underserved buildings
   - Place new RCPs to maximize coverage of underserved population
   ```bash
   snakemake --use-conda data/derived_data/workflow/rcps_clustering_iso.gpkg
   ```

3. **Linear Optimization**
   - Formulate and solve a facility location problem to minimize the number of new RCPs
   - Consider coverage constraints to ensure all buildings are served within threshold time
   ```bash
   snakemake --use-conda data/derived_data/workflow/rcps_optimisation.gpkg
   ```

### 3. Evaluation and Analysis

1. **Walking Duration Calculations**
   - Calculate walking durations from each building to nearest RCP for each solution
   - Compare solutions by coverage percentage and population impact
   ```bash
   snakemake --use-conda data/derived_data/workflow/flats_duration_opt.gpkg
   ```

2. **Method Comparison**
   - Compare different optimization approaches (clustering vs. linear optimization)
   - Evaluate solutions based on coverage metrics and efficiency
   ```bash
   snakemake --use-conda data/derived_data/workflow/method_comparison.csv
   ```

3. **Sensitivity Analysis**
   - Analyze impact of clustering parameters on solution quality
   - Determine optimal number of clusters for best coverage/cost ratio
   ```bash
   snakemake --use-conda run_sensitivity_analysis
   ```

4. **Facility Count Analysis (p-analysis)**
   - Analyze impact of varying the number of new facilities (p-value)
   - Identify diminishing returns in coverage improvement
   ```bash
   snakemake --use-conda run_p_analysis
   ```

### 4. Visualization

1. **Suitability Maps**
   - Visualize suitable locations for new RCPs
   - Highlight existing and potential RCP sites

2. **Walking Duration Maps**
   - Color-coded visualization of walking time to nearest RCP
   - Identify underserved areas (>10 min walking time)

3. **Isochrone Maps**
   - Visualize service areas based on walking time
   - Show coverage of existing and proposed RCPs

4. **Comparison Visualizations**
   - Compare different optimization methods
   - Visualize coverage improvements and population impact

## Installation and Usage

### Prerequisites
- Conda/Miniconda
- Git
- A local instance of [Valhalla](https://github.com/valhalla/valhalla/) or [OpenRouteService](https://openrouteservice.org/); an ORS API key works too but is too slow for a real analysis
- The Gurobi solver speeds things up a bit but is optional

### Setup
1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd rcp_project
   ```

2. Create and activate the base environment:
   ```bash
   conda env create -f envs/snake_env.yaml
   conda activate snake_env
   ```

### Running the Workflow
- Complete workflow:
  ```bash
  snakemake --use-conda
  ```

- Specific target:
  ```bash
  snakemake --use-conda <target-file>
  ```

- Sensitivity analysis:
  ```bash
  snakemake --use-conda run_sensitivity_analysis
  ```

- p-Analysis (number of facilities):
  ```bash
  snakemake --use-conda run_p_analysis
  ```

### Configuration
Adjust parameters in `config/config.yaml`:
- `routing_engine`: Choose routing service (valhalla or ors)
- `data_preparation`: Parameters for data preparation steps
- `optimization`: Parameters for the optimization models
- `p_values`: Values for p-analysis (number of facilities)

## Results and Outputs

The workflow produces the following key outputs:

1. **Optimized RCP Locations**
   - Proposed locations for new RCPs under different optimization strategies
   - Files: `rcps_clustering_ors.gpkg`, `rcps_clustering_iso.gpkg`, `rcps_optimisation.gpkg`

2. **Coverage Analysis**
   - Buildings with calculated walking time to nearest RCP
   - Population covered within different time thresholds
   - Files: `flats_duration_current.gpkg`, etc.

3. **Method Comparison**
   - Comparison metrics between different optimization approaches
   - Efficiency analysis (population served per new RCP)
   - Files: `method_comparison.csv`, `method_comparison_table.tex`

4. **Visualization Maps**
   - Interactive HTML maps showing analysis results
   - Files in `data/plots/` directory

## Routing Engines
The project supports two routing engines:

1. **Valhalla** (default)
   - High-performance open-source routing engine
   - Requires local installation, see config file

2. **OpenRouteService (ORS)**
   - Alternative routing service with API access
   - Requires API key configuration or local instance, see config file

To switch routing engines, modify the `routing_engine` parameter in `config/config.yaml`.

> **IMPORTANT CAVEAT**: Currently, the routing engine selection in the configuration file does not work properly. The main branch only works with Valhalla, not with ORS. If you need to use OpenRouteService, please use the `generate_ors_results` branch instead.
>
> This project was developed as a master thesis at the Global Health Engineering group in collaboration with ERZ (Entsorgung + Recycling Zürich). The routing engine issues are expected to be fixed when the work is published in a future paper.
