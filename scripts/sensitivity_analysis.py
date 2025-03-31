import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.ticker as mticker
import yaml
import os
import sys

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import load_config

# Load config
config = load_config()

# Get cluster numbers from config
CLUSTERS = config["sensitivity_analysis"]["clusters"]

# Get input and output files from snakemake
duration_files = snakemake.input.duration_files
output_summary = snakemake.output.summary
output_plot = snakemake.output.plot

# Analyze results
results = []

# First pass: calculate metrics for each file
for duration_file in duration_files:
    # Extract cluster number from filename
    n = int(Path(duration_file).stem.split('_')[-1])
    
    durations = gpd.read_file(duration_file)
    
    # Calculate metrics
    avg_duration = (durations['duration_min'] * durations['est_pop']).sum() / durations['est_pop'].sum()
    pop_outside_10min = durations[durations['duration_min'] > 10]['est_pop'].sum()
    
    results.append({
        'clusters': n,
        'avg_duration': avg_duration,
        'pop_outside_10min': pop_outside_10min
    })

# Create results DataFrame
results_df = pd.DataFrame(results)

# Sort by cluster count for consistent plotting
results_df = results_df.sort_values('clusters')

# Plot results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(results_df['clusters'], results_df['avg_duration'], marker='o')
ax1.set_title(r'Average Walking Duration by $n_{clusters}$')
ax1.set_xlabel('Number of Demand Clusters')
ax1.set_ylabel('Weighted Average Duration (min:s)')

def format_min_sec(x, pos):
    minutes = int(x)
    seconds = int(round((x - minutes) * 60))
    return f"{minutes}:{seconds:02d}"

ax1.yaxis.set_major_formatter(mticker.FuncFormatter(format_min_sec))

ax2.plot(results_df['clusters'], results_df['pop_outside_10min'], marker='o')
ax2.set_title(r'Population Outside 10min by $n_clusters$')
ax2.set_xlabel('Number of Demand Clusters')
ax2.set_ylabel('Population Outside 10min')

plt.tight_layout()
plt.savefig(output_plot, dpi=350)

# Save results to CSV
results_df.to_csv(output_summary, index=False)
