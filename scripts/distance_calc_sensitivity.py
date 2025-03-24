import geopandas as gpd
import logging
from sklearn.neighbors import BallTree
import numpy as np
import scripts.util as util

# Set up logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    # Read input data
    logger.info("Reading input data")
    flats = gpd.read_file(snakemake.input.flats)
    rcps = gpd.read_file(snakemake.input.rcps)

    #convert to epsg 4326
    flats = flats.to_crs(epsg=4326)
    rcps = rcps.to_crs(epsg=4326)
    
    # Valhalla routing service URL
    valhalla_url = "http://localhost:8002/route"
    
    # Create a BallTree for efficient nearest neighbor search
    tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps, 'ID')
    
    # Calculate distances for each flat to nearest RCP
    logger.info(f"Calculating distances for {len(flats)} flats")
    buildings = flats.copy()
    buildings['nearest_rcp_id'], buildings['duration'] = zip(*buildings['geometry'].apply(
        lambda geom: util.find_nearest_rcp_duration(geom, tree, rcp_coords, rcp_ids, valhalla_url)
    ))
    
    # Save the output
    logger.info(f"Saving output to {snakemake.output.duration}")
    buildings.to_file(snakemake.output.duration, driver='GPKG')
    
    logger.info("Distance calculation completed successfully")
except Exception as e:
    logger.error(f"Error in distance calculation: {e}")
    raise