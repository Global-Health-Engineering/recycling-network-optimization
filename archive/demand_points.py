import numpy as np
import geopandas as gpd
import folium
from sklearn.cluster import KMeans
from shapely.geometry import Point

# n_clusters should be obtained from Snakemake parameters
n_clusters = int(snakemake.params.n_clusters)

# Load buildings population data from Snakemake input
buildings_pop = gpd.read_file(snakemake.input.flats)

# Load RCP data from Snakemake input
rcp_data = gpd.read_file(snakemake.input.rcps)

# Ensure correct CRS
buildings_pop = buildings_pop.to_crs("EPSG:4326")
rcp_data = rcp_data.to_crs("EPSG:4326")

print("Data loaded successfully. Shape:", buildings_pop.shape)

# Extract coordinates and population weights
coordinates = np.column_stack([buildings_pop.geometry.x, buildings_pop.geometry.y])
weights = buildings_pop['est_pop'].values

# Initialize KMeans with n clusters
kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)

# Fit KMeans with sample weights
kmeans.fit(coordinates, sample_weight=weights)

# Assign cluster labels to buildings_pop
buildings_pop['kmeans_cluster'] = kmeans.labels_

# Get cluster centers
cluster_centers = kmeans.cluster_centers_

# Create a GeoDataFrame for cluster centers
cluster_centers_gdf = gpd.GeoDataFrame(
    {'cluster_id': range(len(cluster_centers))},
    geometry=[Point(xy) for xy in cluster_centers],
    crs=buildings_pop.crs
)

# Compute population per cluster
cluster_pop = buildings_pop.groupby('kmeans_cluster')['est_pop'].sum().reset_index()
cluster_centers_gdf = cluster_centers_gdf.merge(cluster_pop, left_on='cluster_id', right_on='kmeans_cluster', how='left')
cluster_centers_gdf.rename(columns={'est_pop': 'total_est_pop'}, inplace=True)

# Create a Folium map centered around the mean coordinates
m_clusters = folium.Map(
    location=[buildings_pop.geometry.y.mean(), buildings_pop.geometry.x.mean()],
    zoom_start=12,
    control_scale=True,
    tiles='cartodbpositron'
)

# Add cluster centers to the map with population
for _, row in cluster_centers_gdf.iterrows():
    folium.CircleMarker(
        location=[row.geometry.y, row.geometry.x],
        radius=2,
        color='blue',
        fill=True,
        fill_color='blue',
        fill_opacity=0.6,
        popup=f"Cluster {row.cluster_id}<br>Population: {int(row.total_est_pop)}"
    ).add_to(m_clusters)

# Add a legend (optional)
legend_html = '''
<div style="position: fixed; 
            bottom: 50px; left: 50px; width: 150px; height: 60px; 
            border:2px solid grey; z-index:9999; font-size:14px;
            background-color:white;
            ">
    &nbsp;<b>Cluster Centers</b><br>
    &nbsp;<i style="background: blue; width: 10px; height: 10px; display: inline-block;"></i>&nbsp; Cluster Center
</div>
'''
m_clusters.get_root().html.add_child(folium.Element(legend_html))

# Save and display the map
m_clusters.save(snakemake.output.html_map)
cluster_centers_gdf.to_file(snakemake.output.gpkg, driver='GPKG')
