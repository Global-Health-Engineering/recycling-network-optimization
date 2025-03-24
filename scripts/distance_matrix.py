import requests
from snakemake.logging import logger
import geopandas as gpd
import pandas as pd

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

def calculate_walking_distance_matrix(potential_locations, demand_points, valhalla_url="http://localhost:8002/route"):
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
                # Prepare Valhalla request
                valhalla_params = {
                    "locations": [
                        {"lat": loc.geometry.y, "lon": loc.geometry.x},
                        {"lat": dp.geometry.y, "lon": dp.geometry.x}
                    ],
                    "costing": "pedestrian",
                    "directions_options": {
                        "units": "kilometers"
                    },
                    "format": "json"
                }
                
                # Make the request to Valhalla
                response = requests.post(valhalla_url, json=valhalla_params)
                
                if response.status_code == 200:
                    valhalla_data = response.json()
                    duration_seconds = valhalla_data["trip"]["legs"][0]["summary"]["time"]
                    duration_minutes = duration_seconds / 60  # Convert seconds to minutes
                    
                    data.append({
                        'ID': loc['ID'],
                        'cluster_ID': dp['cluster_id'],
                        'Walking_Duration_Minutes': duration_minutes
                    })
                else:
                    logger.warning(f"Error calculating walking duration: {response.status_code}")
                    logger.warning(response.text)
                    
        walking_df = pd.DataFrame(data)
        return walking_df
    except Exception as e:
        logger.error(f"Error calculating walking duration matrix: {e}")
        return None

# Main execution
try:
    # Setup Valhalla URL
    valhalla_url = "http://localhost:8002/route"
    
    # Calculate and save walking distance matrix
    walking_matrix = calculate_walking_distance_matrix(potential_locations_gdf, demand_points_gdf, valhalla_url)
    if walking_matrix is not None:
        walking_matrix.to_csv(snakemake.output.matrix_walking, index=False)
        logger.info("Calculated and saved walking distance matrix.")
    else:
        logger.error("Failed to calculate walking distance matrix.")
finally:
    logger.info("Script finished.")