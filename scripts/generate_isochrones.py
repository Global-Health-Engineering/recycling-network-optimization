import geopandas as gpd
import requests
from snakemake.logging import logger
from shapely.ops import unary_union
from shapely.geometry import shape
import pandas as pd
import sys
import os
import openrouteservice
from openrouteservice.exceptions import ApiError

# Add path to import utility functions
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scripts.util import generate_isochrone, merge_isochrones_preserve_time, get_ors_client

# Obtain paths from snakemake
INPUT_FLATS = snakemake.input['flats']
INPUT_RCPS = snakemake.input['rcps']
OUTPUT_GPKG_all = snakemake.output['iso_all']
OUTPUT_GPKG_merged = snakemake.output['iso_merged']

# Get routing engine from params
ROUTING_ENGINE = snakemake.params.get('routing_engine', 'valhalla')
logger.info(f"Using {ROUTING_ENGINE} routing engine for isochrones")

# Time limits in seconds
TIME_LIMITS = [60, 120, 180, 240, 300, 360, 420, 480, 540, 600]  # 1 to 10 minutes

# Initialize ORS client if using ORS
ors_client = None
if ROUTING_ENGINE == "ors":
    ors_client = get_ors_client()
    logger.info("Initialized OpenRouteService client for isochrones")

# Load data
flats = gpd.read_file(INPUT_FLATS).to_crs("EPSG:4326")
rcps = gpd.read_file(INPUT_RCPS).to_crs("EPSG:4326")

# Generate isochrones for each RCP
all_isochrones = []
for idx, rcp in rcps.iterrows():
    logger.info(f"Processing RCP {idx+1}/{len(rcps)}")
    point_coord = (rcp.geometry.x, rcp.geometry.y)
    
    for time_limit in TIME_LIMITS:
        logger.info(f"  Generating {time_limit/60} minute isochrone")
        
        if ROUTING_ENGINE == "valhalla":
            # Use existing Valhalla implementation
            isochrone_data = generate_isochrone(point_coord, time_limit)
        else:
            # Use ORS client
            try:
                isochrone_data = ors_client.isochrones(
                    locations=[[point_coord[0], point_coord[1]]],
                    profile='foot-walking',
                    range=[time_limit],  # time in seconds
                    units='m',
                    location_type='start',
                )
            except ApiError as e:
                logger.error(f"ORS API error: {e}")
                continue
            except Exception as e:
                logger.error(f"Error generating isochrone: {e}")
                continue
        
        if isochrone_data is None:
            logger.warning(f"Failed to generate isochrone for RCP {idx+1} with time {time_limit/60} minutes")
            continue
            
        # Process the isochrone data
        if ROUTING_ENGINE == "valhalla":
            # Process Valhalla response
            # ...existing code for Valhalla...
            pass
        else:
            # Process ORS response
            for feature in isochrone_data['features']:
                feature_properties = feature['properties']
                range_value = feature_properties['value']
                
                # Create a polygon from the geometry
                isochrone_polygon = shape(feature['geometry'])
                
                # Add to the list of all isochrones
                all_isochrones.append({
                    'geometry': isochrone_polygon,
                    'time': range_value,
                    'rcp_id': rcp['poi_id'],
                    'rcp_address': rcp['adresse'] if 'adresse' in rcp else ''
                })

# Create GeoDataFrame from all isochrones
if all_isochrones:
    isochrones_gdf = gpd.GeoDataFrame(all_isochrones, crs="EPSG:4326")
    
    # Save all individual isochrones
    isochrones_gdf.to_file(OUTPUT_GPKG_all, driver="GPKG")
    
    # Merge isochrones by time
    merged_isochrones = merge_isochrones_preserve_time(isochrones_gdf)
    merged_isochrones.to_file(OUTPUT_GPKG_merged, driver="GPKG")
    
    logger.info(f"Successfully generated {len(all_isochrones)} isochrones and {len(merged_isochrones)} merged isochrones")
else:
    logger.error("No valid isochrones were generated!")
