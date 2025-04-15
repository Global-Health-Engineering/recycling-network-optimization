import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.ticker as mticker
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
optimality_gap_files = snakemake.input.optimality_gap_files
output_summary = snakemake.output.summary
output_plot = snakemake.output.plot
output_optimality_plot = snakemake.output.optimality_plot

# Analyze results
results = []

# First pass: calculate metrics for each file
for duration_file, gap_file in zip(duration_files, optimality_gap_files):
    # Extract cluster number from filename
    n = int(Path(duration_file).stem.split('_')[-1])
    
    durations = gpd.read_file(duration_file)
    
    # Calculate metrics
    avg_duration = (durations['duration_min'] * durations['est_pop']).sum() / durations['est_pop'].sum()
    pop_outside_10min = durations[durations['duration_min'] > 10]['est_pop'].sum()
    
    # Read optimality gap value from file
    try:
        with open(gap_file, 'r') as f:
            optimality_gap = float(f.read().strip())
    except Exception as e:
        print(f"Warning: Could not read optimality gap for n={n}: {e}")
        optimality_gap = np.nan
    
    results.append({
        'clusters': n,
        'avg_duration': avg_duration,
        'pop_outside_10min': pop_outside_10min,
        'optimality_gap': optimality_gap
    })

# Create results DataFrame
results_df = pd.DataFrame(results)

# Sort by cluster count for consistent plotting
results_df = results_df.sort_values('clusters')

# Plot results (first figure - metrics)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(results_df['clusters'], results_df['avg_duration'], marker='o')
ax1.set_title(r'Average Walking Duration by $n_{clusters}$')
ax1.set_xlabel('Number of Demand Clusters')
ax1.set_ylabel('Weighted Average Duration (min:s)')
ax1.grid(True, linestyle='--', alpha=0.7)

def format_min_sec(x, pos):
    minutes = int(x)
    seconds = int(round((x - minutes) * 60))
    return f"{minutes}:{seconds:02d}"

ax1.yaxis.set_major_formatter(mticker.FuncFormatter(format_min_sec))

ax2.plot(results_df['clusters'], results_df['pop_outside_10min'], marker='o', color='orangered')
ax2.set_title(r'Population Outside 10min by $n_{clusters}$')
ax2.set_xlabel('Number of Demand Clusters')
ax2.set_ylabel('Population Outside 10min')
ax2.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig(output_plot, dpi=350)

# Create separate plot for optimality gap (second figure)
fig2, ax3 = plt.subplots(figsize=(8, 6))

# Plot optimality gap
ax3.plot(results_df['clusters'], results_df['optimality_gap'], marker='o', linestyle='-', color='green')
ax3.set_title('Optimality Gap by Number of Clusters')
ax3.set_xlabel('Number of Demand Clusters')
ax3.set_ylabel('Optimality Gap')
ax3.grid(True, linestyle='--', alpha=0.7)

# Format y-axis as percentage if values are between 0-1
if results_df['optimality_gap'].max() <= 1.0:
    ax3.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))

plt.tight_layout()
plt.savefig(output_optimality_plot, dpi=350)

# Save results to CSV
results_df.to_csv(output_summary, index=False)
