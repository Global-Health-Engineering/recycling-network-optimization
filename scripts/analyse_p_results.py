import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import matplotlib.ticker as mticker
import json
from shapely.geometry import Point, mapping

# Get input and output files from snakemake
duration_files = snakemake.input.duration_files
output_summary = snakemake.output.summary
output_plot = snakemake.output.plot

# Define p values from file names
p_values = [int(Path(file).stem.split('_')[-1]) for file in duration_files]

# Analyze results
results = []
site_info = {}

# First pass: calculate metrics for each p value
for p, duration_file in zip(p_values, duration_files):
    durations = gpd.read_file(duration_file)
    
    # Calculate metrics
    avg_duration = (durations['duration'] * durations['est_pop']).sum() / durations['est_pop'].sum()
    pop_outside_10min = durations[durations['duration'] > 10]['est_pop'].sum()
    
    # Store basic metrics
    results.append({
        'p_value': p,
        'avg_duration': avg_duration,
        'pop_outside_10min': pop_outside_10min
    })
    
    # Get corresponding RCP site file
    site_file = str(duration_file).replace('flats_duration_p_', 'rcps_optimisation_')
    site_file = Path(site_file.replace('/p-analysis/', '/sensitivity_clusters/'))
    
    if site_file.exists():
        sites = gpd.read_file(site_file)
        
        # Store site information
        site_info[p] = {
            'count': len(sites),
            'sites': sites.to_json()
        }
    else:
        print(f"Warning: Site file not found for p={p}: {site_file}")

# Create results DataFrame
results_df = pd.DataFrame(results)

# Create a GeoDataFrame with all sites and their p values
all_sites = []
for p, info in site_info.items():
    if 'sites' in info:
        sites = gpd.read_file(json.dumps(info['sites']))
        sites['p_value'] = p
        all_sites.append(sites)

if all_sites:
    all_sites_gdf = pd.concat(all_sites, ignore_index=True)
    
    # Convert to regular dataframe with geometry as WKT
    all_sites_df = pd.DataFrame(all_sites_gdf.drop(columns='geometry'))
    all_sites_df['geometry'] = all_sites_gdf.geometry.apply(lambda x: x.wkt)
    
    # Merge with results
    results_df = results_df.merge(all_sites_df.groupby('p_value')['id'].apply(list).reset_index(),
                                on='p_value', how='left')
    results_df = results_df.merge(all_sites_df.groupby('p_value')['geometry'].apply(list).reset_index(),
                                 on='p_value', how='left')
else:
    results_df['id'] = None
    results_df['geometry'] = None

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

# Add marginal improvement percentages
for i in range(1, len(results_df)):
    prev_pop = results_df.iloc[i-1]['pop_outside_10min']
    curr_pop = results_df.iloc[i]['pop_outside_10min']
    
    if prev_pop > 0:  # Avoid division by zero
        pct_improvement = ((prev_pop - curr_pop) / prev_pop) * 100
        ax2.annotate(f"{pct_improvement:.1f}%", 
                   xy=(results_df.iloc[i]['p_value'], curr_pop),
                   xytext=(0, -15), 
                   textcoords='offset points',
                   ha='center', fontsize=8)

plt.tight_layout()
plt.savefig(output_plot, dpi=350)

# Save results to CSV - convert lists to strings for CSV compatibility
results_df['site_ids'] = results_df['id'].apply(lambda x: ','.join(map(str, x)) if isinstance(x, list) else '')
results_df['site_geometries'] = results_df['geometry'].apply(lambda x: '|'.join(map(str, x)) if isinstance(x, list) else '')
results_df = results_df.drop(columns=['id', 'geometry'])

results_df.to_csv(output_summary, index=False)

print(f"Analysis complete. Results saved to {output_summary} and {output_plot}")
