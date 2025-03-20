#!/usr/bin/env python
import numpy as np
import geopandas as gpd
import fiona
import pandas as pd
import rasterio
from rasterio.warp import reproject, calculate_default_transform
from rasterio.enums import Resampling
from rasterstats import zonal_stats
import folium
from shapely.wkt import loads, dumps
from pathlib import Path
import tempfile

DST_CRS = "EPSG:4326"

def compute_slope(elev_path: str) -> str:
    with rasterio.open(elev_path) as src:
        elevation = src.read(1)
        transform = src.transform
        dx = transform.a
        dy = -transform.e  # pixel height (absolute value)
        grad_y, grad_x = np.gradient(elevation, abs(dy), abs(dx))
        slope_array = np.arctan(np.sqrt(grad_x**2 + grad_y**2)) * 180 / np.pi

        slope_temp_path = tempfile.mktemp(suffix='.tif')
        profile = src.profile
        profile.update(dtype=rasterio.float32, count=1)
        with rasterio.open(slope_temp_path, "w", **profile) as dst:
            dst.write(slope_array.astype(rasterio.float32), 1)
    return slope_temp_path

def reproject_raster(src_path: str, dst_path: str, dst_crs: str = DST_CRS) -> str:
    with rasterio.open(src_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        with rasterio.open(dst_path, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest
                )
    return dst_path

def suitability_analysis(
    buffer_dist_residential: float,
    building_footprints: gpd.GeoDataFrame,
    buffer_dist_vbz: float,
    area_threshold: float,
    buffer_trees: float,
    max_slope: float,
    parking_lots: gpd.GeoDataFrame,
    slope_raster,   # expects a rasterio dataset with slope data in DST_CRS
    trees: gpd.GeoDataFrame,
    vbz_lines: gpd.GeoDataFrame,
    vbz_points: gpd.GeoDataFrame
):
    parking_filtered = parking_lots[
        ~parking_lots['parking'].isin(['underground', 'multi-storey'])
    ].copy()

    vbz_line_buffers = vbz_lines.buffer(buffer_dist_vbz)
    vbz_point_buffers = vbz_points.buffer(buffer_dist_vbz)
    tree_buffers = trees.buffer(buffer_trees)
    buffer_buildings = building_footprints.geometry.apply(lambda x: x.buffer(buffer_dist_residential))

    all_buffers = gpd.GeoSeries(
        list(vbz_line_buffers) + list(vbz_point_buffers) +
        list(tree_buffers) +
        list(buffer_buildings),
        crs=parking_lots.crs
    ).unary_union

    slope_stats = zonal_stats(
        parking_filtered,
        slope_raster.name,
        stats=['mean', 'max'],
        nodata=slope_raster.nodata
    )

    parking_filtered['slope_mean'] = [s['mean'] for s in slope_stats]
    parking_filtered['slope_max'] = [s['max'] for s in slope_stats]

    slope_filtered = parking_filtered[parking_filtered['slope_mean'] <= max_slope]
    final_areas = slope_filtered.geometry.difference(all_buffers)

    result = gpd.GeoDataFrame(geometry=final_areas, crs=parking_lots.crs)
    result['area'] = result.geometry.area
    suitable = result[result['area'] >= area_threshold]

    return suitable

if __name__ == "__main__":
    elev_path = str(snakemake.input.elevation_model)
    slope_temp_path = compute_slope(elev_path)

    # Reproject computed slope to EPSG:4326
    slope_reproj_path = tempfile.mktemp(suffix='.tif')
    reproject_raster(slope_temp_path, slope_reproj_path, DST_CRS)
    slope_raster = rasterio.open(slope_reproj_path)

    # Load datasets and reproject to EPSG:4326
    tree_dataset = gpd.read_file(str(snakemake.input.trees)).to_crs(DST_CRS)
    parking_lots = gpd.read_file(str(snakemake.input.parking)).to_crs(DST_CRS)
    rcps = gpd.read_file(str(snakemake.input.rcps)).to_crs(DST_CRS)
    building_footprints = gpd.read_file(str(snakemake.input.buildings)).to_crs(DST_CRS)

    building_footprints = building_footprints[building_footprints['art_txt'].str.contains('wohn', case=False, na=False)]
    building_footprints['geometry'] = building_footprints.geometry.apply(
        lambda geom: loads(dumps(geom, output_dimension=2))
    )

    layers = fiona.listlayers(str(snakemake.input.vbz))
    line_layers = []
    point_layers = []

    for layer in layers:
        gdf = gpd.read_file(str(snakemake.input.vbz), layer=layer)
        gdf = gdf.to_crs(DST_CRS)
        if gdf.geom_type.str.startswith("Line").all():
            line_layers.append(gdf)
        elif gdf.geom_type.str.startswith("Point").all():
            point_layers.append(gdf)

    vbz_lines = gpd.GeoDataFrame(
        pd.concat(line_layers, ignore_index=True), crs=line_layers[0].crs
    )
    vbz_points = gpd.GeoDataFrame(
        pd.concat(point_layers, ignore_index=True), crs=point_layers[0].crs
    )

    buffer_dist_vbz = snakemake.params.get("buffer_dist_vbz", 2)
    buffer_trees = snakemake.params.get("buffer_trees", 2)
    max_slope = snakemake.params.get("max_slope", 5)
    area_threshold = snakemake.params.get("area_threshold", 16)
    buffer_buildings = snakemake.params.get("buffer_buildings", 14)

    suitable_areas = suitability_analysis(
        building_footprints=building_footprints,
        buffer_dist_residential=buffer_buildings,
        area_threshold=area_threshold,
        buffer_dist_vbz=buffer_dist_vbz,
        buffer_trees=buffer_trees,
        max_slope=max_slope,
        parking_lots=parking_lots,
        slope_raster=slope_raster,
        trees=tree_dataset,
        vbz_lines=vbz_lines,
        vbz_points=vbz_points
    )

    existing_sites = rcps.copy()
    existing_sites['ID'] = ['e_{}'.format(i+1) for i in existing_sites.index]
    existing_sites['status'] = 'open'
    existing_sites = existing_sites[['geometry', 'ID', 'status']]

    buffer_union = existing_sites.geometry.buffer(125).unary_union

    filtered_potential_sites = suitable_areas[~suitable_areas.geometry.intersects(buffer_union)].copy()
    filtered_potential_sites['ID'] = ['p_{}'.format(i+1) for i in filtered_potential_sites.index]
    filtered_potential_sites['status'] = 'potential'
    filtered_potential_sites = filtered_potential_sites[['geometry', 'ID', 'status']]

    merged_sites = pd.concat([existing_sites, filtered_potential_sites], ignore_index=True)

    if not merged_sites['geometry'].geom_type.str.startswith('Point').all():
        merged_sites = gpd.GeoDataFrame(merged_sites, geometry=merged_sites.geometry.centroid, crs=merged_sites.crs)
    merged_sites['geometry'] = merged_sites.geometry.centroid

    merged_sites.to_file(str(snakemake.output.sites), driver='GPKG')

    m = folium.Map(location=[47.3769, 8.5417], zoom_start=12)
    folium.GeoJson(
        merged_sites,
        name='Suitable Areas',
        popup=folium.GeoJsonPopup(fields=['ID', 'status']),
    ).add_to(m)
    m.save(str(snakemake.output.map))
