import geopandas as gpd
import os
import pandas as pd
import openrouteservice
import logging
import time

# Configure logging
logging.basicConfig(
    filename='distance_calc.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

# Main script
start_time = time.perf_counter()

try:
    os.chdir("/home/silas/projects/msc_thesis")
    logging.info("Changed working directory.")

    # Set parameters
    n = 100
    buffer_distance = 600
    api_key= "5b3ce3597851110001cf624865e19fb4d0c2400e9aba8877785f6853"

    # Import datasets
    flats_zh = gpd.read_file('./data/raw_data/geodata_stadt_Zuerich/building_stats/data/ssz.gwr_stzh_wohnungen.shp')
    rcps = gpd.read_file('./data/raw_data/geodata_stadt_Zuerich/recycling_sammelstellen/data/stzh.poi_sammelstelle_view.shp')
    logging.info("Imported datasets.")

    # Filter data
    flats_zh_existing = flats_zh.query('wstatlang=="Bestehend"').drop_duplicates(subset=['egid'])
    flats_zh_existing['egid'] = flats_zh_existing['egid'].astype(int)
    logging.info("Filtered existing flats.")

    # Check CRS
    logging.info(f"Flats CRS: {flats_zh_existing.crs}")
    logging.info(f"RCPs CRS: {rcps.crs}")

    # Create subset and buffer
    flats_subset = flats_zh_existing.iloc[1:n].copy()
    flats_subset['buffer'] = flats_subset.geometry.buffer(buffer_distance)
    logging.info("Created buffer around flats.")

    # Find points in buffer
    results = pd.DataFrame(columns=['flat_id', 'rcp'])
    for idx, row in flats_subset.iterrows():
        points_in_buffer = rcps[rcps.geometry.within(row['buffer'])]
        if not points_in_buffer.empty:
            for point_idx in points_in_buffer.index:
                results.loc[len(results)] = {'flat_id': row['egid'], 'rcp': point_idx}
    logging.info("Mapped RCPs to flats within buffer.")

    # Initialize ORS client
    client = openrouteservice.Client(key=api_key)
    flats_subset = flats_subset.to_crs(epsg=4326)
    rcps = rcps.to_crs(epsg=4326)
    logging.info("Initialized OpenRouteService client and transformed CRS.")

    # Calculate routes
    results['distance'] = 0.0
    results['duration'] = 0.0

    for idx, row in results.iterrows():
            flat_coords = flats_subset.loc[flats_subset['egid'] == row['flat_id'], 'geometry'].values[0]
            rcp_coords = rcps.geometry[row['rcp']]
            coords = ([flat_coords.x, flat_coords.y], [rcp_coords.x, rcp_coords.y])
            route = client.d5b3ce3597851110001cf624865e19fb4d0c2400e9aba8877785f6853irections(coordinates=coords, profile='foot-walking', format='geojson')
            distance = route['features'][0]['properties']['segments'][0]['distance']
            duration = route['features'][0]['properties']['segments'][0]['duration']
            results.at[idx, 'distance'] = distance
            results.at[idx, 'duration'] = duration / 60

    logging.info("Calculated walking routes.")

    # Map closest RCP
    closest_rcp = results.loc[results.groupby('flat_id')['duration'].idxmin()]
    flats_subset_with_rcp = flats_subset.merge(closest_rcp[['flat_id', 'rcp', 'distance', 'duration']], 
                                              left_on='egid', right_on='flat_id', how='left')
    flats_subset_with_rcp.drop(columns=['buffer', 'flat_id'], inplace=True)
    flats_subset_with_rcp.to_file('./data/derived_data/flats_subset_with_rcp.shp')
    logging.info("Mapped closest RCPs and saved shapefile.")
    logging.info("Process completed.")
    elapsed_time = time.perf_counter() - start_time
    minutes, seconds = divmod(elapsed_time, 60)
    logging.info(f"Elapsed time: {int(minutes)} minutes and {int(seconds)} seconds.")

except Exception as e:
    logging.critical(f"An unexpected error occurred: {e}")
