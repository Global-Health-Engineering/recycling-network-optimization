# Recycling-network-optimization

Masters thesis project by Silas Schweizer, supervised by Dr. Jakub Tkaczuk and suüpported by Nicolas Seeman. 

# Goal
- Test different methods of optimising the recycling network of Zurich.
- Identify optimal locations for new recycling collection points (RCPs)
- Maximize population coverage and minimize walking distances
- Compare different approaches to network optimization

# Methods
The workflow uses the open source ORS routing algorithm to calculate realistic walking distances to the closest recycling collection point. Based on this, a linear optimisation algorithm selects the most suitable locations for p additional recycling collection points.

The project implements and compares several approaches:
- Isochrone-based clustering
- Spatial clustering based on route calculation
- P-median linear optimization
- Sensitivity analysis with varying cluster sizes

# Data Sources
- Building footprints and population data from Stadt Zürich
- Public transport infrastructure (VBZ)
- Elevation models for slope analysis
- Tree and data by Grün Stadt Zürich
- Existing recycling collection points

# Project Structure
- `rules/`: Snakemake workflow rules
- `scripts/`: Python scripts for spatial analysis and optimization
- `envs/`: Conda environment configurations
- `data/`: Input and output data (raw and derived)

# Installation
1. Clone the repository
2. Install dependencies using conda: `conda env create -f envs/snake_env.yaml`
3. Activate the environment: `conda activate snake_env.yaml`

# Usage
Run the complete workflow:
```bash
snakemake --use-conda
```

Run the sensitivity analysis for the number of k-means clusters:

```bash
snakemake --use-conda run_sensitivity_analysis
```

Run the analysis for the number of recycling facilites to be installed and the expected impact:

```bash
snakemake --use-conda run_p_analysis
```
# Timeline
The model is currently under development, a working version and documentation will be released in May 2025.

