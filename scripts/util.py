def calculate_duration(origin, destination, client):
    """
    Calculate walking duration between origin and destination using ORS.

    Parameters:
    - origin: Tuple of (longitude, latitude)
    - destination: Tuple of (longitude, latitude)
    - client: OpenRouteService client instance

    Returns:
    - Duration in minutes if successful, else None
    """
    try:
        route = client.directions(
            coordinates=[origin, destination],
            profile='foot-walking',
            format='geojson'
        )
        duration_seconds = route['features'][0]['properties']['segments'][0]['duration']
        return duration_seconds / 60  # Convert to minutes
    except openrouteservice.exceptions.ApiError as e:
        logging.error(f"ORS API error for origin {origin} to destination {destination}: {e}")
    except Exception as e:
        logging.error(f"Unexpected error for origin {origin} to destination {destination}: {e}")
    return None

def find_nearest_rcp_duration(flat_geom):
    """
    Find the nearest RCP within a 5 km radius and calculate walking duration.

    Parameters:
    - flat_geom: Shapely geometry of the flat

    Returns:
    - Tuple of (rcp_id, duration_min)
    """
    flat_coord = (flat_geom.x, flat_geom.y)
    flat_rad = np.radians([flat_coord[::-1]])  # [lat, lon] in radians
    distance, index = tree.query(flat_rad, k=5)  # Get top 5 nearest

    for dist, idx in zip(distance[0], index[0]):
        actual_distance = dist * 6371000  # Earth radius in meters
        if actual_distance <= 5000:  # 5 km radius
            rcp_coord = rcp_coords[idx]
            duration = calculate_duration(flat_coord, rcp_coord, client)
            if duration is not None:
                return rcp_ids[idx], duration
    return None, None

