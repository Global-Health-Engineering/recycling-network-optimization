import geopandas as gpd
import pandas as pd
import os
import openrouteservice
import util

def create_high_pop_unserved():
    """
    Create a GeoDataFrame of flats with high population that are not served within
    the isochrone threshold, including their walking duration to nearest RCP.
    """
    # Import required data
    flats_with_pop = gpd.read_file('/home/silas/projects/msc_thesis/data/derived_data/flats_population.gpkg')
    merged_isochrones_gdf = gpd.read_file('/home/silas/projects/msc_thesis/data/derived_data/isochrones_1-10min.gpkg')
    rcps = gpd.read_file('/home/silas/projects/msc_thesis/data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp')

    # Remove flats with population 0
    flats_with_pop = flats_with_pop[flats_with_pop['est_pop'] > 0]

    # Convert CRS to EPSG:4326
    flats_with_pop.to_crs(epsg=4326, inplace=True)
    rcps.to_crs(epsg=4326, inplace=True)

    # Initialize BallTree
    tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps)

    # Verify 'time' column exists
    if 'time' not in merged_isochrones_gdf.columns:
        raise KeyError("'time' column is missing in merged_isochrones_gdf")

    # Spatial join to retain all flats
    joined = gpd.sjoin(
        flats_with_pop, 
        merged_isochrones_gdf[['geometry', 'time']], 
        how='left', 
        predicate='within'
    )

    # Assign default high time value to unserved flats
    iso_threshold = 10
    joined['time'] = joined['time'].fillna(iso_threshold + 1)

    # Get the shortest time for each flat
    joined = joined.groupby('egid', as_index=False).agg({
        'est_pop': 'first',
        'geometry': 'first',
        'time': 'min'
    })

    # Identify unserved flats
    high_pop_unserved = joined[joined['time'] >= iso_threshold].copy()

    # Get ORS key from environment variable
    ors_key = os.getenv('ORS_API_KEY')
    client = openrouteservice.Client(key=ors_key)

    # Calculate duration to the nearest RCP
    high_pop_unserved[['nearest_rcp_id', 'duration_to_rcp_min']] = high_pop_unserved['geometry'].apply(
        lambda geom: util.find_nearest_rcp_duration(geom, tree, rcp_coords, rcp_ids, client)
    ).apply(pd.Series)

    # Update 'time' with the calculated duration
    joined.loc[joined['time'] >= iso_threshold, 'time'] = high_pop_unserved['duration_to_rcp_min'].values

    # Ensure GeoDataFrame consistency
    high_pop_unserved = gpd.GeoDataFrame(high_pop_unserved, geometry='geometry', crs="EPSG:4326")

    # Optional: Update 'time' column in 'high_pop_unserved'
    high_pop_unserved['time'] = high_pop_unserved['duration_to_rcp_min']
    high_pop_unserved.drop(columns=['nearest_rcp_id', 'duration_to_rcp_min'], inplace=True)

    # Export to file
    high_pop_unserved.to_file(
        '/home/silas/projects/msc_thesis/data/derived_data/high_pop_unserved_with_durations.gpkg', 
        driver='GPKG'
    )

    return high_pop_unserved

if __name__ == "__main__":
    high_pop_unserved = create_high_pop_unserved()
