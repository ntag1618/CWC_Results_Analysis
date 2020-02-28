# Import Raw Lidar layers - merge and clip to extent. resample and align with Sep17 DSM.
# Create Canopy height model from DSM and DTM.
# Identify canopy >2m and classify
# create polygon of woodland area.

import rasterio
from rasterio.merge import merge
from rasterio.plot import show
from rasterio import features
from rasterio.mask import mask
import geopandas as gpd
from matplotlib import pyplot as plt
from shapely.geometry import box
from shapely.geometry import shape
import os
from glob import glob
import json

raw_dtm_folder = os.path.abspath('C:/HG_Projects/CWC_Drone_work/DTM_Raw/Download_1462811/england-dtm-2m_3439972/st')
chm_out_path = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/CWC_CHM_Sep17.tif')
trees_vec_out = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/Rip_Vec_Sep17.gpkg')

home = os.path.abspath("C:/HG_Projects/CWC_Drone_work/Prec_Anal_Exports/Rasters_v3")
dsm_path = os.path.join(home, "dsm3.tif")

def main():

    asc_files = glob(os.path.join(raw_dtm_folder,'*.asc'))

    src_files_to_mosaic = []
    for asc in asc_files:
        src = rasterio.open(asc)
        src_files_to_mosaic.append(src)

    mosaic, out_trans = merge(src_files_to_mosaic)

    # show(mosaic, cmap='terrain')

    with rasterio.open(dsm_path) as src:
        dsm_arr = src.read(1)
        bounds = src.bounds
        bbox = box(bounds.left, bounds.bottom, bounds.right, bounds.top)
        geo = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=src.crs.data)
        coords = getFeatures(geo)
        dsm_meta = src.meta
        loc_crs = src.crs
        # show(src.read(1, masked=True), cmap='terrain', transform=src.transform)

    tmpfile1 = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/temp_ras.tif')
    with rasterio.open(tmpfile1, 'w+', count=1, dtype='float32', driver='GTiff', transform=out_trans,
                      height=mosaic.shape[-2], width=mosaic.shape[-1]) as dataset:
        dataset.crs = loc_crs
        dataset.write(mosaic)
        # show(dataset.read(1), cmap='viridis', transform=dataset.transform)

        out_img, out_transform = mask(dataset, shapes=coords, crop=True)

    out_img = out_img.reshape((out_img.shape[-2], out_img.shape[-1]))

    new_meta = dsm_meta.copy()
    new_meta.update(count=3)
    new_meta.update(dtype='float32')

    with rasterio.open(chm_out_path, 'w+', **new_meta) as clip_dtm:

        clip_dtm.write(out_img, 1) # writing the array like this automates the resampling process

        clip_arr = clip_dtm.read(1)
        canopy_arr = clip_arr.copy()

        canopy_arr[dsm_arr == -999] = -999
        canopy_arr[dsm_arr != -999] = dsm_arr[dsm_arr != -999] - canopy_arr[dsm_arr != -999]
        canopy_arr[canopy_arr > 80] = -999

        can_class1 = canopy_arr.copy()
        can_class1[canopy_arr > 2] = 1
        can_class1[canopy_arr < 2] = 0

        clip_dtm.write(canopy_arr, 1)
        clip_dtm.write(can_class1, 2)

        show(clip_dtm.read(2, masked=True), cmap='viridis', transform=clip_dtm.transform)

        masker = can_class1 == 1
        generator = features.shapes(can_class1, mask=masker, transform=clip_dtm.transform)

    geom_list = []
    for geom, value in generator:
        geom = shape(geom)
        geom_list.append(geom)

    canopy_gdf = gpd.GeoDataFrame(
            gpd.GeoSeries(geom_list),
            columns=['geometry'], crs=loc_crs.data)

    canopy_gdf = canopy_gdf.buffer(1, join_style=1).buffer(-1, join_style=1)

    canopy_gdf = canopy_gdf[canopy_gdf.area > 0.5]

    canopy_gdf.to_file(trees_vec_out, driver="GPKG")

    canopy_gdf.plot(color='green')
    plt.show()

    if os.path.isfile(tmpfile1):
        os.remove(tmpfile1)

    print("BREAK")

def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features'][0]['geometry']]



if __name__ == '__main__':
    main()
