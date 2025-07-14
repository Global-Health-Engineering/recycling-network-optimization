#!/usr/bin/env python3
"""
Script to compare Euclidean distance results with routing engine results.
This script helps analyze the differences between straight-line distance calculations
and actual routing-based duration calculations.
"""

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def load_comparison_data(routing_file, euclidean_file):
    """Load and merge routing and euclidean data for comparison."""
    
    # Load both datasets
    routing_data = gpd.read_file(routing_file)
    euclidean_data = gpd.read_file(euclidean_file)
    
    # Merge on geometry or index (assuming same flats in same order)
    # Using index merge since both should have the same flats
    comparison_df = pd.DataFrame({
        'routing_duration': routing_data['duration_min'],
        'euclidean_duration': euclidean_data['duration_min'],
        'routing_rcp': routing_data['nearest_rcp'],
        'euclidean_rcp': euclidean_data['nearest_rcp']
    })
    
    # Calculate differences
    comparison_df['duration_difference'] = comparison_df['routing_duration'] - comparison_df['euclidean_duration']
    comparison_df['duration_ratio'] = comparison_df['routing_duration'] / comparison_df['euclidean_duration']
    comparison_df['same_rcp'] = comparison_df['routing_rcp'] == comparison_df['euclidean_rcp']
    
    return comparison_df

def generate_comparison_plots(comparison_df, method_name, output_dir):
    """Generate comparison plots."""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style and colors
    plt.style.use('default')
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']  # Default matplotlib colors
    
    # 1. Scatter plot: Routing vs Euclidean durations
    plt.figure(figsize=(10, 8))
    plt.scatter(comparison_df['euclidean_duration'], comparison_df['routing_duration'], 
                alpha=0.6, s=20, color=colors[0])
    
    # Add diagonal line for perfect correlation
    max_val = max(comparison_df['euclidean_duration'].max(), comparison_df['routing_duration'].max())
    plt.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Perfect correlation')
    
    plt.xlabel('Euclidean Duration (minutes)')
    plt.ylabel('Routing Duration (minutes)')
    plt.title(f'Routing vs Euclidean Durations - {method_name}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/duration_scatter_{method_name}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Histogram of duration differences
    plt.figure(figsize=(10, 6))
    plt.hist(comparison_df['duration_difference'], bins=50, alpha=0.7, edgecolor='black', color=colors[1])
    plt.xlabel('Duration Difference (Routing - Euclidean) [minutes]')
    plt.ylabel('Frequency')
    plt.title(f'Distribution of Duration Differences - {method_name}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/duration_diff_hist_{method_name}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Box plot of duration ratios
    plt.figure(figsize=(8, 6))
    box_plot = plt.boxplot(comparison_df['duration_ratio'], vert=True, patch_artist=True)
    box_plot['boxes'][0].set_facecolor(colors[2])
    box_plot['boxes'][0].set_alpha(0.7)
    plt.ylabel('Duration Ratio (Routing / Euclidean)')
    plt.title(f'Distribution of Duration Ratios - {method_name}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/duration_ratio_box_{method_name}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Combined comparison plot with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Subplot 1: Scatter plot
    ax1.scatter(comparison_df['euclidean_duration'], comparison_df['routing_duration'], 
                alpha=0.6, s=15, color=colors[0])
    max_val = max(comparison_df['euclidean_duration'].max(), comparison_df['routing_duration'].max())
    ax1.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='Perfect correlation')
    ax1.set_xlabel('Euclidean Duration (minutes)')
    ax1.set_ylabel('Routing Duration (minutes)')
    ax1.set_title('Routing vs Euclidean Durations')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Subplot 2: Duration differences histogram
    ax2.hist(comparison_df['duration_difference'], bins=30, alpha=0.7, edgecolor='black', color=colors[1])
    ax2.set_xlabel('Duration Difference (minutes)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Distribution of Duration Differences')
    ax2.grid(True, alpha=0.3)
    
    # Subplot 3: Duration ratios histogram
    ax3.hist(comparison_df['duration_ratio'], bins=30, alpha=0.7, edgecolor='black', color=colors[2])
    ax3.set_xlabel('Duration Ratio (Routing / Euclidean)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Distribution of Duration Ratios')
    ax3.grid(True, alpha=0.3)
    
    # Subplot 4: RCP selection agreement
    same_rcp_pct = comparison_df['same_rcp'].mean() * 100
    different_rcp_pct = 100 - same_rcp_pct
    ax4.pie([same_rcp_pct, different_rcp_pct], 
            labels=[f'Same RCP\n({same_rcp_pct:.1f}%)', f'Different RCP\n({different_rcp_pct:.1f}%)'],
            colors=[colors[3], colors[4]], autopct='%1.1f%%', startangle=90)
    ax4.set_title('RCP Selection Agreement')
    
    plt.suptitle(f'Comparison Summary - {method_name}', fontsize=16, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.savefig(f'{output_dir}/comparison_summary_{method_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

def print_comparison_statistics(comparison_df, method_name):
    """Print statistical summary of the comparison."""
    
    print(f"\n=== Comparison Statistics - {method_name} ===")
    print(f"Number of flats: {len(comparison_df)}")
    print(f"Same RCP selected: {comparison_df['same_rcp'].sum()} ({comparison_df['same_rcp'].mean()*100:.1f}%)")
    
    print(f"\nDuration Statistics:")
    print(f"  Routing - Mean: {comparison_df['routing_duration'].mean():.2f} min, Median: {comparison_df['routing_duration'].median():.2f} min")
    print(f"  Euclidean - Mean: {comparison_df['euclidean_duration'].mean():.2f} min, Median: {comparison_df['euclidean_duration'].median():.2f} min")
    
    print(f"\nDuration Differences (Routing - Euclidean):")
    print(f"  Mean: {comparison_df['duration_difference'].mean():.2f} min")
    print(f"  Median: {comparison_df['duration_difference'].median():.2f} min")
    print(f"  Std: {comparison_df['duration_difference'].std():.2f} min")
    print(f"  Min: {comparison_df['duration_difference'].min():.2f} min")
    print(f"  Max: {comparison_df['duration_difference'].max():.2f} min")
    
    print(f"\nDuration Ratios (Routing / Euclidean):")
    print(f"  Mean: {comparison_df['duration_ratio'].mean():.2f}")
    print(f"  Median: {comparison_df['duration_ratio'].median():.2f}")
    print(f"  Std: {comparison_df['duration_ratio'].std():.2f}")
    
    # Correlation
    correlation = comparison_df['routing_duration'].corr(comparison_df['euclidean_duration'])
    print(f"\nCorrelation: {correlation:.3f}")

def generate_metrics_and_latex(comparison_df, routing_data, euclidean_data, method_name, output_dir):
    """
    Generate metrics and output a LaTeX table comparing routing and euclidean results.
    Metrics:
    - Number of people living outside 10 minutes duration to an RCP
    - Population-weighted mean walking duration
    """
    # Try to get population column name
    pop_col = None
    for candidate in ['population', 'pop', 'POP', 'Population', 'est_pop']:
        if candidate in routing_data.columns:
            pop_col = candidate
            break
    if pop_col is None:
        raise ValueError("Population column not found in input data.")

    # Add population to comparison_df
    comparison_df['population'] = routing_data[pop_col]

    # Number of people living outside 10 min (routing)
    pop_outside_10_routing = comparison_df.loc[comparison_df['routing_duration'] > 10, 'population'].sum()
    # Number of people living outside 10 min (euclidean)
    pop_outside_10_euclidean = comparison_df.loc[comparison_df['euclidean_duration'] > 10, 'population'].sum()

    # Population-weighted mean walking duration (routing)
    pop_weighted_mean_routing = np.average(comparison_df['routing_duration'], weights=comparison_df['population'])
    # Population-weighted mean walking duration (euclidean)
    pop_weighted_mean_euclidean = np.average(comparison_df['euclidean_duration'], weights=comparison_df['population'])

    # Total population
    total_population = comparison_df['population'].sum()

    # Output LaTeX table
    latex_table = f"""
% Comparison of routing and euclidean metrics for {method_name}
\begin{{table}}[ht]
\centering
\begin{{tabular}}{{lrr}}
\toprule
Metric & Routing Engine & Euclidean \\
\midrule
Population outside 10 min & {int(pop_outside_10_routing):,} & {int(pop_outside_10_euclidean):,} \\
Population-weighted mean duration [min] & {pop_weighted_mean_routing:.2f} & {pop_weighted_mean_euclidean:.2f} \\
Total population & {int(total_population):,} & {int(total_population):,} \\
\bottomrule
\end{{tabular}}
\caption{{Comparison of accessibility metrics for routing engine and Euclidean distance (method: {method_name})}}
\label{{tab:comparison_{method_name}}}
\end{{table}}
"""
    # Save to file
    os.makedirs(output_dir, exist_ok=True)
    latex_path = os.path.join(output_dir, f"comparison_metrics_{method_name}.tex")
    with open(latex_path, "w") as f:
        f.write(latex_table)
    print(f"LaTeX table saved to: {latex_path}")
    print(latex_table)

def main():
    parser = argparse.ArgumentParser(description='Compare routing engine and Euclidean distance results')
    parser.add_argument('--method', choices=['iso', 'ors', 'opt'], default='iso',
                        help='Method to compare (iso, ors, or opt)')
    parser.add_argument('--data-dir', default='data/derived_data',
                        help='Data directory path')
    parser.add_argument('--output-dir', default='data/plots/euclidean_comparison',
                        help='Output directory for plots')
    parser.add_argument('--no-plots', action='store_true',
                        help='Skip generating plots')
    
    args = parser.parse_args()
    
    # Construct file paths
    method_map = {
        'iso': 'clustering_iso',
        'ors': 'clustering_ors', 
        'opt': 'opt'
    }
    
    method_name = method_map[args.method]
    routing_file = f"{args.data_dir}/workflow/flats_duration_{method_name}.gpkg"
    euclidean_file = f"{args.data_dir}/euclidean_analysis/flats_duration_{method_name}_euclidean.gpkg"
    
    # Check if files exist
    if not os.path.exists(routing_file):
        print(f"Error: Routing file not found: {routing_file}")
        print("Please run the main workflow first.")
        return
    
    if not os.path.exists(euclidean_file):
        print(f"Error: Euclidean file not found: {euclidean_file}")
        print("Please run the Euclidean analysis first: snakemake euclidean_analysis_all")
        return
    
    # Load and compare data
    print(f"Loading data for method: {method_name}")
    comparison_df = load_comparison_data(routing_file, euclidean_file)
    routing_data = gpd.read_file(routing_file)
    euclidean_data = gpd.read_file(euclidean_file)

    # Print statistics
    print_comparison_statistics(comparison_df, method_name)

    # Generate metrics and LaTeX table
    generate_metrics_and_latex(comparison_df, routing_data, euclidean_data, method_name, args.output_dir)

    # Generate plots
    if not args.no_plots:
        print(f"\nGenerating comparison plots...")
        generate_comparison_plots(comparison_df, method_name, args.output_dir)
        print(f"Plots saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
