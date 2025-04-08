import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.ticker as mticker
import json
from shapely.geometry import Point, mapping
from shapely.geometry import shape

# Get input and output files from snakemake
duration_files = snakemake.input.duration_files
optimality_gap_files = snakemake.input.optimality_gap_files  # Add optimality gap files
output_summary = snakemake.output.summary
output_plot = snakemake.output.plot

# Define p values from file names
p_values = [int(Path(file).stem.split('_')[-1]) for file in duration_files]

# Analyze results
results = []

# First pass: calculate metrics for each p value
for p, duration_file, gap_file in zip(p_values, duration_files, optimality_gap_files):
    durations = gpd.read_file(duration_file)
    
    # Calculate metrics
    avg_duration = (durations['duration_min'] * durations['est_pop']).sum() / durations['est_pop'].sum()
    pop_outside_10min = durations[durations['duration_min'] > 10]['est_pop'].sum()
    
    # Read optimality gap value from file
    try:
        with open(gap_file, 'r') as f:
            optimality_gap = float(f.read().strip())
    except Exception as e:
        print(f"Warning: Could not read optimality gap for p={p}: {e}")
        optimality_gap = np.nan
    
    # Store basic metrics
    results.append({
        'p_value': p,
        'avg_duration': avg_duration,
        'pop_outside_10min': pop_outside_10min,
        'optimality_gap': optimality_gap
    })
    
    # (Original site info collection removed)
    site_file = Path(str(duration_file).replace('flats_duration_p_', 'rcps_optimisation_'))
    if not site_file.exists():
        print(f"Warning: Site file not found for p={p}: {site_file}")

# Create results DataFrame with only the indicator metrics and p values
results_df = pd.DataFrame(results)

# Plot results
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))  # Change to 3 subplots

# Sort by p_value to ensure proper plotting order
results_df = results_df.sort_values('p_value')

# Plot average duration by p value
ax1.plot(results_df['p_value'], results_df['avg_duration'], marker='o', linestyle='-')
ax1.set_title('Average Walking Duration by Number of Facilities')
ax1.set_xlabel('Number of Facilities (p)')
ax1.set_ylabel('Weighted Average Duration (min:s)')
ax1.grid(True, linestyle='--', alpha=0.7)

# Format y-axis to show minutes:seconds
def format_min_sec(x, pos):
    minutes = int(x)
    seconds = int(round((x - minutes) * 60))
    return f"{minutes}:{seconds:02d}"

ax1.yaxis.set_major_formatter(mticker.FuncFormatter(format_min_sec))

# Plot population outside 10min walking distance
ax2.plot(results_df['p_value'], results_df['pop_outside_10min'], marker='o', linestyle='-', color='orangered')
ax2.set_title('Population Outside 10min Walking Distance')
ax2.set_xlabel('Number of Facilities (p)')
ax2.set_ylabel('Population Count')
ax2.grid(True, linestyle='--', alpha=0.7)

# Plot optimality gap
ax3.plot(results_df['p_value'], results_df['optimality_gap'], marker='o', linestyle='-', color='green')
ax3.set_title('Optimality Gap by Number of Facilities')
ax3.set_xlabel('Number of Facilities (p)')
ax3.set_ylabel('Optimality Gap')
ax3.grid(True, linestyle='--', alpha=0.7)

# Format y-axis as percentage if values are between 0-1
if results_df['optimality_gap'].max() <= 1.0:
    ax3.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))

plt.tight_layout()
plt.savefig(output_plot, dpi=350)

# Save summary results (only indicators and p value)
results_df.to_csv(output_summary, index=False)

print(f"Analysis complete. Results saved to {output_summary} and {output_plot}")
