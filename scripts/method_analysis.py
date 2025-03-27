#!/usr/bin/env python
import geopandas as gpd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mdutils.mdutils import MdUtils
import logging
import sys

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(snakemake.log[0]), logging.StreamHandler(sys.stdout)]
)

def main():
    logging.info("Starting method analysis")
    
    # Load the data
    logging.info("Loading input datasets")
    duration_orig = gpd.read_file(snakemake.input.duration_current)
    duration_iso = gpd.read_file(snakemake.input.duration_iso)
    duration_ors = gpd.read_file(snakemake.input.duration_ors)
    duration_opt = gpd.read_file(snakemake.input.duration_opt)
    
    # Calculate population weighted duration for each dataset
    logging.info("Calculating weighted durations")
    orig_weighted = (duration_orig['duration'] * duration_orig['est_pop']).sum() / duration_orig['est_pop'].sum()
    iso_weighted = (duration_iso['duration_min'] * duration_iso['est_pop']).sum() / duration_iso['est_pop'].sum()
    ors_weighted = (duration_ors['duration_min'] * duration_ors['est_pop']).sum() / duration_ors['est_pop'].sum()
    opt_weighted = (duration_opt['duration_min'] * duration_opt['est_pop']).sum() / duration_opt['est_pop'].sum()
    
    # Calculate population outside 10-minute range for each dataset
    logging.info("Calculating populations outside network coverage")
    orig_outside = duration_orig[duration_orig['duration'] > 10]['est_pop'].sum()
    iso_outside = duration_iso[duration_iso['duration_min'] > 10]['est_pop'].sum()
    ors_outside = duration_ors[duration_ors['duration_min'] > 10]['est_pop'].sum()
    opt_outside = duration_opt[duration_opt['duration_min'] > 10]['est_pop'].sum()
    
    # Calculate total population
    total_pop = duration_orig['est_pop'].sum()
    
    # Calculate coverage percentage for each dataset
    orig_coverage = ((total_pop - orig_outside) / total_pop) * 100
    iso_coverage = ((total_pop - iso_outside) / total_pop) * 100
    ors_coverage = ((total_pop - ors_outside) / total_pop) * 100
    opt_coverage = ((total_pop - opt_outside) / total_pop) * 100
    
    # Function to format time in min:sec
    def format_time(minutes):
        mins = int(minutes)
        secs = int((minutes - mins) * 60)
        return f"{mins}:{secs:02d}"
    
    # Create a DataFrame with the comparison results
    # Load RPC files to count new RCP sites for each method
    rpc_iso = gpd.read_file(snakemake.input.rpc_iso)
    rpc_ors = gpd.read_file(snakemake.input.rpc_ors)
    rpc_opt = gpd.read_file(snakemake.input.rpc_opt)
    rcp_current = gpd.read_file(snakemake.input.rcp_current)

    additional_iso = len(rpc_iso)-len(rcp_current)
    additional_ors = len(rpc_ors)-len(rcp_current)
    additional_opt = len(rpc_opt)-len(rcp_current)

    comparison_df = pd.DataFrame({
        'Method': ['Current Situation', 'Iso+Clustering', 'ORS+Clustering', 'Linear Optimisation'],
        'Population Outside 10min': [orig_outside, iso_outside, ors_outside, opt_outside],
        'Coverage (%)': [orig_coverage, iso_coverage, ors_coverage, opt_coverage],
        'Average Walking Time (min)': [format_time(orig_weighted), format_time(iso_weighted), 
                                       format_time(ors_weighted), format_time(opt_weighted)],
        'Number of additional RCPs': [0, additional_iso, additional_ors, additional_opt]
    })
    
    # Save the comparison results to CSV
    logging.info("Saving comparison results to CSV")
    comparison_df.to_csv(snakemake.output.comparison_csv, index=False)
    
    # Font size parameters
    TITLE_SIZE = 19
    LABEL_SIZE = 14
    TICK_SIZE = 14
    
    # Set dark theme style
    plt.style.use('dark_background')
    
    # Create a figure with three subplots
    logging.info("Creating comparison plots")
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 8))
    fig.patch.set_facecolor('#1a1a1a')
    
    # Define blue color palette
    colors = ['#1f77b4', '#2b94d1', '#37b0ee']
    
    # Plot 1: People outside the network
    ax1.bar(comparison_df['Method'], comparison_df['Population Outside 10min'], color=colors)
    ax1.set_title('People Outside Network', pad=15, fontsize=TITLE_SIZE)
    ax1.set_ylabel('Number of People', fontsize=LABEL_SIZE)
    ax1.tick_params(axis='x', rotation=45, labelsize=TICK_SIZE)
    ax1.set_facecolor('#1a1a1a')
    
    # Plot 2: Walking time
    walking_times = [float(t.split(':')[0]) + float(t.split(':')[1])/60 for t in comparison_df['Average Walking Time (min)']]
    ax2.bar(comparison_df['Method'], walking_times, color=colors)
    ax2.set_title('Population Weighted Walking Time', pad=15, fontsize=TITLE_SIZE)
    ax2.set_ylabel('Minutes', fontsize=LABEL_SIZE)
    ax2.tick_params(axis='x', rotation=45, labelsize=TICK_SIZE)
    ax2.set_facecolor('#1a1a1a')
    
    # Plot 3: Number of New RCPs
    ax3.bar(comparison_df['Method'][1:], comparison_df['Number of additional RCPs'][1:], color=colors[:-1])
    ax3.set_title('Number of New RCPs', pad=15, fontsize=TITLE_SIZE)
    ax3.set_ylabel('Count', fontsize=LABEL_SIZE)
    ax3.tick_params(axis='x', rotation=45, labelsize=TICK_SIZE)
    ax3.set_facecolor('#1a1a1a')
    
    plt.tight_layout()
    
    # Save comparison plot
    logging.info("Saving method comparison plot")
    plt.savefig(snakemake.output.comparison_plot, dpi=400, bbox_inches='tight', facecolor='#1a1a1a')
    plt.close()
    
    # Calculate people brought in by comparing to current situation
    logging.info("Creating efficiency plot")
    people_brought_in = comparison_df.loc[0, 'Population Outside 10min'] - comparison_df.loc[1:, 'Population Outside 10min']
    new_rcps = comparison_df.loc[1:, 'Number of additional RCPs']
    
    # Calculate people per RCP
    people_per_rcp = people_brought_in / new_rcps
    
    # Create figure with dark theme
    plt.figure(figsize=(10, 6), facecolor='#1a1a1a')
    
    # Plot the data
    plt.bar(comparison_df['Method'][1:], people_per_rcp, color=colors[0:2])
    
    # Customize the plot
    plt.title('People Brought Into Network per New RCP', fontsize=TITLE_SIZE, color='white', pad=15)
    plt.ylabel('Number of People per RCP', fontsize=LABEL_SIZE, color='white')
    plt.xticks(rotation=45, fontsize=TICK_SIZE, color='white')
    plt.gca().set_facecolor('#1a1a1a')
    plt.gca().tick_params(colors='white')
    plt.gca().spines['bottom'].set_color('white')
    plt.gca().spines['top'].set_color('white')
    plt.gca().spines['left'].set_color('white')
    plt.gca().spines['right'].set_color('white')
    
    plt.tight_layout()
    
    # Save efficiency plot
    logging.info("Saving efficiency plot")
    plt.savefig(snakemake.output.efficiency_plot, dpi=400, bbox_inches='tight', facecolor='#1a1a1a')
    plt.close()
    
    # Create markdown file with comparison table
    logging.info("Creating markdown table")
    mdFile = MdUtils(file_name=snakemake.output.markdown_table)
    
    # Add title
    mdFile.write("\n## Method Comparison Results\n")
    
    # Create table header 
    headers = comparison_df.columns.tolist()
    table = [headers]
    
    # Add rows
    for _, row in comparison_df.iterrows():
        table.append([str(val) for val in row.values])
    
    # Create table with mdutils
    mdFile.new_table(columns=len(headers), rows=len(table), text=sum(table, []))
    
    # Create file
    mdFile.create_md_file()
    
    logging.info("Method analysis completed successfully")

if __name__ == "__main__":
    main()
