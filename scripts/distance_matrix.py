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

# Datasets
potential_locations = snakemake.input.potential_locations
demand_points = snakemake.input.demand_points

# Read datasets
potential_locations_gdf = gpd.read_file(potential_locations)
demand_points_gdf = gpd.read_file(demand_points)
# Ensure demand_points is a GeoDataFrame
if isinstance(demand_points, gpd.GeoSeries):
    demand_points = demand_points.to_frame().reset_index(drop=True)
# Ensure CRS is EPSG:4326
if potential_locations_gdf.crs != "EPSG:4326":
    potential_locations_gdf = potential_locations_gdf.to_crs("EPSG:4326")
    logger.info("Converted potential_locations CRS to EPSG:4326.")

if demand_points_gdf.crs != "EPSG:4326":
    demand_points_gdf = demand_points_gdf.to_crs("EPSG:4326")
    logger.info("Converted demand_points CRS to EPSG:4326.")

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

def calculate_walking_distance_matrix(potential_locations_gdf, demand_points_gdf):
    """Calculate walking duration matrix between potential locations and demand points GeoDataFrames"""
    try:
        data = []
        for i, loc_row in potential_locations_gdf.iterrows():
            loc = loc_row.geometry.centroid  # Get centroid of location geometry
            for j, dp_row in demand_points_gdf.iterrows():
                dp = dp_row.geometry  # Get demand point geometry
                response = client.directions(
                    coordinates=[[loc.x, loc.y], [dp.x, dp.y]],
                    profile='foot-walking',
                    format='geojson'
                )
                duration = response['features'][0]['properties']['segments'][0]['duration']
                data.append({
                    'Location_ID': i,
                    'Object_ID': loc_row['object_id'],
                    'Demand_Point_ID': j,
                    'Walking_Duration': duration
                })
        walking_df = pd.DataFrame(data)
        return walking_df
    except Exception as e:
        logger.error(f"Error calculating walking duration matrix: {e}")
        return None

# Main execution
try:
    os.chdir("/home/silas/rcp_project/rcp_project")
    logger.info("Changed working directory.")
    
    # Read and ensure CRS matches EPSG:4326
    gdf = gpd.read_file(snakemake.input[0])
    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
        logger.info("Converted CRS to EPSG:4326.")
    
    source_coords = gdf.geometry.tolist()
    
    # Initialize client
    client = get_ors_client()
    
    # Calculate matrix
    matrix = calculate_distance_matrix(source_coords, DEPOT_COORDS)
    if matrix is not None:
        matrix.to_csv(snakemake.output.matrix_trucks, index=False)
        logger.info("Calculated and saved truck distance matrix.")
    else:
        logger.error("Failed to calculate distance matrix.")
    
    # Constants for testing - adjust subset size as needed
    TEST_SUBSET_SIZE = 5  # Number of locations to process during testing

    # Calculate and save walking distance matrix
    # Use head() to subset the data for testing
    walking_matrix = calculate_walking_distance_matrix(
        potential_locations_gdf.head(TEST_SUBSET_SIZE),
        demand_points_gdf.head(TEST_SUBSET_SIZE)
    )
    if walking_matrix is not None:
        walking_matrix.to_csv(snakemake.output.matrix_walking, index=False)
        logger.info(f"Calculated and saved walking distance matrix (using {TEST_SUBSET_SIZE} test locations).")
    else:
        logger.error("Failed to calculate walking distance matrix.")
finally:
    logger.info("Script finished.")