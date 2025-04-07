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
output_summary = snakemake.output.summary
output_plot = snakemake.output.plot

# Define p values from file names
p_values = [int(Path(file).stem.split('_')[-1]) for file in duration_files]

# Analyze results
results = []

# First pass: calculate metrics for each p value
for p, duration_file in zip(p_values, duration_files):
    durations = gpd.read_file(duration_file)
    
    # Calculate metrics
    avg_duration = (durations['duration_min'] * durations['est_pop']).sum() / durations['est_pop'].sum()
    pop_outside_10min = durations[durations['duration_min'] > 10]['est_pop'].sum()
    
    # Store basic metrics
    results.append({
        'p_value': p,
        'avg_duration': avg_duration,
        'pop_outside_10min': pop_outside_10min
    })
    
    # (Original site info collection removed)
    site_file = Path(str(duration_file).replace('flats_duration_p_', 'rcps_optimisation_'))
    if not site_file.exists():
        print(f"Warning: Site file not found for p={p}: {site_file}")

# Create results DataFrame with only the indicator metrics and p values
results_df = pd.DataFrame(results)

# Plot results
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

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


plt.tight_layout()
plt.savefig(output_plot, dpi=350)

# Save summary results (only indicators and p value)
results_df.to_csv(output_summary, index=False)

print(f"Analysis complete. Results saved to {output_summary} and {output_plot}")
