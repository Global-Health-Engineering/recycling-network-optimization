#!/usr/bin/env python3
"""
Script to generate comprehensive comparison between Euclidean distance and routing engine results.
This script is designed to be run as part of the Snakemake workflow and outputs a LaTeX table
comparing all RCP methods (current, iso, ors, opt).
"""

import geopandas as gpd
import pandas as pd
import numpy as np
import os
from snakemake.logging import logger

def format_time(minutes):
    """Format duration in minutes to mm:ss format."""
    total_seconds = int(minutes * 60)
    mins = total_seconds // 60
    secs = total_seconds % 60
    return f"{mins}:{secs:02d}"

def load_and_calculate_metrics(routing_file, euclidean_file, method_name):
    """Load files and calculate metrics for a specific method."""
    if not os.path.exists(routing_file):
        logger.warning(f"Routing file not found for {method_name}: {routing_file}")
        return None, None, None, None
    
    if not os.path.exists(euclidean_file):
        logger.warning(f"Euclidean file not found for {method_name}: {euclidean_file}")
        return None, None, None, None
    
    # Load data
    routing_data = gpd.read_file(routing_file)
    euclidean_data = gpd.read_file(euclidean_file)
    
    # Calculate metrics
    routing_mean, routing_outside = calculate_metrics(routing_data, method_name, data_type='routing')
    euclidean_mean, euclidean_outside = calculate_metrics(euclidean_data, method_name, data_type='euclidean')
    
    return routing_mean, routing_outside, euclidean_mean, euclidean_outside

def calculate_metrics(data, method_name, data_type='routing'):
    """
    Calculate population-weighted mean walking duration and population outside 10 min.
    
    Parameters:
    - data: GeoDataFrame with the data
    - method_name: Method name (current, clustering_iso, etc.)
    - data_type: Either 'routing' or 'euclidean'
    
    Returns:
    - Tuple of (pop_weighted_mean, pop_outside_10)
    """
    # Find population column (similar to method_analysis.py logic)
    pop_col = None
    for col in ['est_pop', 'pop_total', 'population', 'pop']:
        if col in data.columns:
            pop_col = col
            break
    
    if pop_col is None:
        raise ValueError("Population column not found in input data.")

    # Determine duration column based on method and data type
    if method_name == 'current' and data_type == 'routing':
        duration_col = 'duration'  # Current routing uses 'duration' column
    else:
        duration_col = 'duration_min'  # All euclidean and other routing use 'duration_min'
    
    if duration_col not in data.columns:
        raise ValueError(f"Duration column '{duration_col}' not found in {method_name} {data_type} data.")

    # Population-weighted mean walking duration
    total_pop = data[pop_col].sum()
    pop_weighted_mean = (data[duration_col] * data[pop_col]).sum() / total_pop
    
    # Number of people living outside 10 min
    pop_outside_10 = data.loc[data[duration_col] > 10, pop_col].sum()
    
    return pop_weighted_mean, pop_outside_10

# Get file paths from snakemake
data_dir = snakemake.params.data_dir
output_file = snakemake.output[0]

# Define file paths for all methods
methods = [
    ('current', 'current'),
    ('clustering_iso', 'clustering_iso'), 
    ('clustering_ors', 'clustering_ors'),
    ('opt', 'opt')
]

logger.info("Generating comprehensive Euclidean vs Routing comparison...")

# Collect results for all methods
results = []

for short_name, long_name in methods:
    # Construct file paths
    routing_file = f"{data_dir}/workflow/flats_duration_{long_name}.gpkg"
    euclidean_file = f"{data_dir}/euclidean_analysis/flats_duration_{short_name}_euclidean.gpkg"
    
    logger.info(f"Processing method: {short_name}")
    
    # Calculate metrics
    routing_mean, routing_outside, euclidean_mean, euclidean_outside = load_and_calculate_metrics(
        routing_file, euclidean_file, short_name
    )
    
    if routing_mean is not None:
        results.append({
            'method': short_name,
            'method_display': short_name.replace('_', ' ').title(),
            'routing_mean': routing_mean,
            'routing_outside': routing_outside,
            'euclidean_mean': euclidean_mean,
            'euclidean_outside': euclidean_outside
        })
        logger.info(f"Processed {short_name}: routing_mean={routing_mean:.2f}, euclidean_mean={euclidean_mean:.2f}")
    else:
        logger.warning(f"Skipped {short_name} due to missing files")

if not results:
    logger.error("No valid results found for any method")
    raise ValueError("No valid comparison data found")

# Generate comprehensive LaTeX table
latex_table = """% Comprehensive comparison of routing engine and Euclidean distance results
\\begin{table}[ht]
\\centering
\\begin{tabular}{lrrrr}
\\toprule
& \\multicolumn{2}{c}{Population outside 10 min} & \\multicolumn{2}{c}{Pop.-weighted mean duration [min:sec]} \\\\
\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}
Method & Routing & Euclidean & Routing & Euclidean \\\\
\\midrule
"""

for result in results:
    method_display = result['method_display']
    routing_outside = int(result['routing_outside'])
    euclidean_outside = int(result['euclidean_outside']) 
    routing_time = format_time(result['routing_mean'])
    euclidean_time = format_time(result['euclidean_mean'])
    
    latex_table += f"{method_display} & {routing_outside:,} & {euclidean_outside:,} & {routing_time} & {euclidean_time} \\\\\n"

latex_table += """\\bottomrule
\\end{tabular}
\\caption{Comparison of accessibility metrics between routing engine and Euclidean distance calculations for all RCP methods}
\\label{tab:euclidean_routing_comparison}
\\end{table}
"""

# Create output directory if needed
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Save to output file
with open(output_file, 'w') as f:
    f.write(latex_table)

logger.info(f"Comprehensive comparison table saved to: {output_file}")

# Also print to console for verification
print("\n" + "="*80)
print("COMPREHENSIVE EUCLIDEAN VS ROUTING COMPARISON")
print("="*80)
print(latex_table)
print("="*80)

logger.info("Euclidean comparison analysis completed successfully")
