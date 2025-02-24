import openrouteservice
import os
from snakemake.logging import logger
import sys
import geopandas as gpd
import pandas as pd

# Constants
DUMP_COORDS = [8.512281878574365, 47.38447647508825]  # [longitude, latitude]
TRUCK_GARAGE_COORDS = [8.575500, 47.414889]  # [longitude, latitude]
DEPOT_COORDS = [DUMP_COORDS, TRUCK_GARAGE_COORDS]

# Read GeoDataFrames
potential_locations_gdf = gpd.read_file(snakemake.input.potential_locations)
if not isinstance(potential_locations_gdf, gpd.GeoDataFrame):
    potential_locations_gdf = gpd.GeoDataFrame(
        potential_locations_gdf, geometry=potential_locations_gdf.geometry
    )

demand_points_gdf = gpd.read_file(snakemake.input.demand_points)
if not isinstance(demand_points_gdf, gpd.GeoDataFrame):
    demand_points_gdf = gpd.GeoDataFrame(
        demand_points_gdf, geometry=demand_points_gdf.geometry
    )

# Make sure the CRS is EPSG:4326
potential_locations_gdf = potential_locations_gdf.to_crs("EPSG:4326")
demand_points_gdf = demand_points_gdf.to_crs("EPSG:4326")

def get_ors_client():
    """Initialize OpenRouteService client"""
    return openrouteservice.Client(base_url='http://localhost:8080/ors')

def calculate_distance_matrix(source_coords, depot_coords):
    """Calculate distance matrix between source coordinates and depot coordinates"""
    try:
        # Convert source coordinates from Points to list format
        source_coords_list = [[p.x, p.y] for p in source_coords]
        
        # depot_coords is already in the correct format
        depot_coords_list = depot_coords  # Already [longitude, latitude] pairs
        
        # Get distances from ORS
        distances = []
        for coord in source_coords_list:
            response = client.directions(
                coordinates=[coord, depot_coords_list[0]],  # To dump
                profile='driving-car',
                format='geojson'
            )
            dump_dist = response['features'][0]['properties']['segments'][0]['distance']
            
            response = client.directions(
                coordinates=[coord, depot_coords_list[1]],  # To garage
                profile='driving-car',
                format='geojson'
            )
            garage_dist = response['features'][0]['properties']['segments'][0]['distance']
            
            distances.append([dump_dist, garage_dist])

        # Create DataFrame
        df = pd.DataFrame(distances, columns=['Distance_to_Dump', 'Distance_to_Garage'])
        return df
    except Exception as e:
        logger.error(f"Error calculating distance matrix: {e}")
        return None

def calculate_walking_distance_matrix(potential_locations, demand_points):
    """Calculate walking duration matrix between potential locations and demand points.
    
    The returned DataFrame includes:
      - ID from potential_locations,
      - cluster_ID from demand_points,
      - Walking_Duration_Minutes for each pair.
    """
    try:
        data = []
        # Iterate over full GeoDataFrames to have access to attributes and geometry.
        for _, loc in potential_locations.iterrows():
            for _, dp in demand_points.iterrows():
                response = client.directions(
                    coordinates=[[loc.geometry.x, loc.geometry.y], [dp.geometry.x, dp.geometry.y]],
                    profile='foot-walking',
                    format='geojson'
                )
                duration = response['features'][0]['properties']['segments'][0]['duration']
                duration_minutes = duration / 60  # Convert seconds to minutes
                data.append({
                    'ID': loc['ID'],
                    'cluster_ID': dp['cluster_id'],
                    'Walking_Duration_Minutes': duration_minutes
                })
        walking_df = pd.DataFrame(data)
        return walking_df
    except Exception as e:
        logger.error(f"Error calculating walking duration matrix: {e}")
        return None

# Main execution
try:
    # Get source coordinates 
    source_coords = potential_locations_gdf.geometry.tolist()
    
    # Initialize client
    client = get_ors_client()
    
    # Calculate matrix
    matrix = calculate_distance_matrix(source_coords, DEPOT_COORDS)
    if matrix is not None:
        matrix.to_csv(snakemake.output.matrix_trucks, index=False)
        logger.info("Calculated and saved truck distance matrix.")
    else:
        logger.error("Failed to calculate truck distance matrix.")
    
    # Calculate and save walking distance matrix
    walking_matrix = calculate_walking_distance_matrix(potential_locations_gdf, demand_points_gdf)
    if walking_matrix is not None:
        walking_matrix.to_csv(snakemake.output.matrix_walking, index=False)
        logger.info("Calculated and saved walking distance matrix.")
    else:
        logger.error("Failed to calculate walking distance matrix.")
finally:
    logger.info("Script finished.")