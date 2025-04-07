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
    rcp_iso = gpd.read_file(snakemake.input.rcp_iso)
    rcp_ors = gpd.read_file(snakemake.input.rcp_ors)
    rcp_opt = gpd.read_file(snakemake.input.rcp_opt)
    rcp_current = gpd.read_file(snakemake.input.rcp_current)

    additional_iso = len(rcp_iso)-len(rcp_current)
    additional_ors = len(rcp_ors)-len(rcp_current)
    additional_opt = len(rcp_opt)-len(rcp_current)

    comparison_df = pd.DataFrame({
        'Method': ['Current Situation', 'Iso+Clustering', 'ORS+Clustering', 'Linear Optimisation'],
        'Population Outside 10min': [orig_outside, iso_outside, ors_outside, opt_outside],
        'Coverage (%)': [orig_coverage, iso_coverage, ors_coverage, opt_coverage],
        'Average Walking Time (min)': [format_time(orig_weighted), format_time(iso_weighted), 
                                       format_time(ors_weighted), format_time(opt_weighted)],
        'Number of additional RCPs': [0, additional_iso, additional_ors, additional_opt]
    })
    
    # Create a new column for numeric walking time
    comparison_df['walking_time'] = comparison_df['Average Walking Time (min)'].apply(
        lambda x: float(x.split(':')[0]) + float(x.split(':')[1]) / 60.0
    )
    
    # Save the comparison results to CSV
    logging.info("Saving comparison results to CSV")
    comparison_df.to_csv(snakemake.output.comparison_csv, index=False)
    
    # Font size parameters
    TITLE_SIZE = 19
    LABEL_SIZE = 14
    TICK_SIZE = 14
    
    # Bright theme parameters
    background_color = 'white'
    text_color = 'black'
    colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']  # Blue, Orange, Green, Purple

    logging.info("Rewriting comparison plots with bright theme")
    plt.figure(figsize=(18, 6), facecolor=background_color)

    # Plot 1: People Outside Network
    plt.subplot(1, 3, 1)
    bars1 = plt.bar(comparison_df['Method'], comparison_df['Population Outside 10min'], color=colors)
    plt.title('People Outside Network', fontsize=TITLE_SIZE, color=text_color, pad=15)
    plt.ylabel('Number of People', fontsize=LABEL_SIZE, color=text_color)
    plt.xticks(rotation=45, fontsize=TICK_SIZE, color=text_color, ha='right')
    plt.gca().set_facecolor(background_color)
    plt.gca().tick_params(colors=text_color)
    for spine in plt.gca().spines.values():
        spine.set_color(text_color)

    # Plot 2: Walking Time
    plt.subplot(1, 3, 2)
    bars2 = plt.bar(comparison_df['Method'], comparison_df['walking_time'], color=colors)
    plt.title('Population Weighted Walking Time', fontsize=TITLE_SIZE, color=text_color, pad=15)
    plt.ylabel('Minutes', fontsize=LABEL_SIZE, color=text_color)
    plt.xticks(rotation=45, fontsize=TICK_SIZE, color=text_color, ha='right')
    plt.gca().set_facecolor(background_color)
    plt.gca().tick_params(colors=text_color)
    for spine in plt.gca().spines.values():
        spine.set_color(text_color)

    # Plot 3: Number of New RCPs
    plt.subplot(1, 3, 3)
    bars3 = plt.bar(comparison_df['Method'][1:], comparison_df['Number of additional RCPs'][1:], color=colors[1:])
    plt.title('Number of New RCPs', fontsize=TITLE_SIZE, color=text_color, pad=15)
    plt.ylabel('Count', fontsize=LABEL_SIZE, color=text_color)
    plt.xticks(rotation=45, fontsize=TICK_SIZE, color=text_color, ha='right')
    plt.gca().set_facecolor(background_color)
    plt.gca().tick_params(colors=text_color)
    for spine in plt.gca().spines.values():
        spine.set_color(text_color)

    plt.tight_layout()

    # Save the combined figure
    plt.savefig(snakemake.output.comparison_plot , dpi=400, bbox_inches='tight', facecolor=background_color)
    
    # Calculate people brought in by comparing to current situation
    logging.info("Creating efficiency plot")
    people_brought_in = comparison_df.loc[0, 'Population Outside 10min'] - comparison_df.loc[1:, 'Population Outside 10min']
    new_rcps = comparison_df.loc[1:, 'Number of additional RCPs']
    
    # Calculate people per RCP
    people_per_rcp = people_brought_in / new_rcps
    
    # Create figure with bright theme
    plt.figure(figsize=(10, 6), facecolor=background_color)
    
    # Plot the data
    plt.bar(comparison_df['Method'][1:], people_per_rcp, color=colors[1:])
    
    # Customize the plot
    plt.title('People Brought Into Network per New RCP', fontsize=TITLE_SIZE, color=text_color, pad=15)
    plt.ylabel('Number of People per RCP', fontsize=LABEL_SIZE, color=text_color)
    plt.xticks(rotation=45, fontsize=TICK_SIZE, color=text_color, ha='right')
    plt.gca().set_facecolor(background_color)
    plt.gca().tick_params(colors=text_color)
    for spine in plt.gca().spines.values():
        spine.set_color(text_color)
    
    plt.tight_layout()
    
    # Save efficiency plot
    logging.info("Saving efficiency plot")
    plt.savefig(snakemake.output.efficiency_plot, dpi=400, bbox_inches='tight', facecolor=background_color)
    plt.close()
    
    # Format columns for LaTeX output
    def time_to_seconds(t):
        mins, secs = t.split(':')
        return int(mins) * 60 + int(secs)

    df_formatted = comparison_df.copy().drop(columns=['walking_time'])
    df_formatted['Coverage (%)'] = df_formatted['Coverage (%)'].apply(lambda x: f"{x:.2f}")
    df_formatted['Population Outside 10min'] = df_formatted['Population Outside 10min'].apply(lambda x: f"{int(round(x, 0))}")
     
    # Create LaTeX file with comparison table
    logging.info("Creating LaTeX table")
    latex_file = snakemake.output.latex_table  # file name remains same, but now it's a .tex file
    with open(latex_file, 'w') as f:
        # Write the LaTeX document header
        f.write('\\documentclass{article}\n')
        f.write('\\usepackage[utf8]{inputenc}\n')
        f.write('\\usepackage{booktabs}\n')
        f.write('\\usepackage{array}\n')
        f.write('\\begin{document}\n')
        f.write('\\section*{Method Comparison Results}\n\n')
        
        # Begin the table environment
        f.write('\\begin{table}[ht]\n\\centering\n')
        
        # Determine column alignment: using left alignment for all columns
        num_cols = len(df_formatted.columns.tolist())
        alignment = "l" * num_cols
        f.write(f'\\begin{{tabular}}{{{alignment}}}\n')
        f.write('\\toprule\n')
        
        # Write table header
        headers = df_formatted.columns.tolist()
        f.write(' & '.join(headers) + ' \\\\\n')
        f.write('\\midrule\n')
        
        # Write table rows
        for _, row in df_formatted.iterrows():
            f.write(' & '.join(str(val) for val in row.values) + ' \\\\\n')
        
        f.write('\\bottomrule\n')
        f.write('\\end{tabular}\n')
        f.write('\\caption{Method Comparison Results}\n')
        f.write('\\end{table}\n\n')
        
        # End of document
        f.write('\\end{document}\n')
    
    logging.info("Method analysis completed successfully")

if __name__ == "__main__":
    main()
