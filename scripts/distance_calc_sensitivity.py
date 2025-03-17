import geopandas as gpd
import logging
from sklearn.neighbors import BallTree
import numpy as np
import scripts.util as util
import openrouteservice as ors

# Set up logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    # Read input data
    logger.info("Reading input data")
    flats = gpd.read_file(snakemake.input.flats)
    rcps = gpd.read_file(snakemake.input.rcps)
    
    # Initialize ORS client
    client = ors.Client(base_url='http://localhost:8080/ors')
    
    # Create a BallTree for efficient nearest neighbor search
    tree, rcp_coords, rcp_ids = util.initialize_ball_tree(rcps, 'id')
    
    # Calculate distances for each flat to nearest RCP
    logger.info(f"Calculating distances for {len(flats)} flats")
    buildings = flats.copy()
    buildings['nearest_rcp_id'], buildings['duration'] = zip(*buildings['geometry'].apply(
        lambda geom: util.find_nearest_rcp_duration(geom, tree, rcp_coords, rcp_ids, client)
    ))
    
    # Save the output
    logger.info(f"Saving output to {snakemake.output.duration}")
    buildings.to_file(snakemake.output.duration, driver='GPKG')
    
    logger.info("Distance calculation completed successfully")
except Exception as e:
    logger.error(f"Error in distance calculation: {e}")
    raise