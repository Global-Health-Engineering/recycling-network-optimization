import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.ticker as mticker


# Define cluster numbers and output path
CLUSTERS = [10, 20, 30, 40, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500]
DERIVED_DATA = Path("/home/silas/rcp_project/rcp_project/data/derived_data")
PLOTS_PATH = Path("/home/silas/rcp_project/rcp_project/data/plots")

# Analyze results
results = []
for n in CLUSTERS:
    duration_file = DERIVED_DATA / f"sensitivity_analysis/flats_duration_{n}.gpkg"
    if not duration_file.exists():
        print(f"Skipping {n} - file not found")
        continue
        
    durations = gpd.read_file(duration_file)
    
    # Calculate metrics
    avg_duration = (durations['duration'] * durations['est_pop']).sum() / durations['est_pop'].sum()
    pop_outside_10min = durations[durations['duration'] > 10]['est_pop'].sum()
    
    results.append({
        'clusters': n,
        'avg_duration': avg_duration,
        'pop_outside_10min': pop_outside_10min
    })

# Create results DataFrame
results_df = pd.DataFrame(results)

# Plot results

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(results_df['clusters'], results_df['avg_duration'], marker='o')
ax1.set_title('Average Walking Duration by Cluster Count')
ax1.set_xlabel('Number of Demand Clusters')
ax1.set_ylabel('Weighted Average Duration (min:s)')

def format_min_sec(x, pos):
    minutes = int(x)
    seconds = int(round((x - minutes) * 60))
    return f"{minutes}:{seconds:02d}"

ax1.yaxis.set_major_formatter(mticker.FuncFormatter(format_min_sec))

ax2.plot(results_df['clusters'], results_df['pop_outside_10min'], marker='o')
ax2.set_title('Population Outside 10min by Cluster Count')
ax2.set_xlabel('Number of Demand Clusters')
ax2.set_ylabel('Population Outside 10min')

plt.tight_layout()
plt.savefig(snakemake.output['plot'], dpi=350)
plt.show()

# Save results to CSV
results_df.to_csv(snakemake.output['summary'], index=False)
