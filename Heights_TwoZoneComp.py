# Compare the relative change of the different woodland zones.

import geopandas as gpd
import rasterio
from rasterio.mask import mask
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from scipy import stats
import pandas as pd
# from tqdm import tqdm

# suppress warnings for now...
import warnings
warnings.filterwarnings("ignore")

home = os.path.abspath("C:/HG_Projects/CWC_Drone_work/Prec_Anal_Exports/Rasters_v4")
shps_root = os.path.abspath("C:/HG_Projects/CWC_Drone_work/shp_files")
shrub_zones = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/Woodland_Zones20m.gpkg')

#
ts1_name = "Dec16"
ts1_path = os.path.join(home, "chm1.tif")
# ts1_fs = os.path.join(shps_root, 'FS_1617.shp')
ts1_fs = os.path.join(shps_root, 'Feed_signs_All.gpkg')
ts1_tup = (ts1_name, ts1_path, ts1_fs)


ts3_name = "Sep17"
ts3_path = os.path.join(home, "chm3.tif")
# ts3_fs = os.path.join(shps_root, 'FS_1617.shp')
ts3_fs =os.path.join(shps_root, 'Feed_signs_All.gpkg')
ts3_tup = (ts3_name, ts3_path, ts3_fs)

#
ts4_name = "Jan18"
ts4_path = os.path.join(home, "chm4.tif")
# ts4_fs = os.path.join(shps_root, 'FS_1718.shp')
ts4_fs = os.path.join(shps_root, 'Feed_signs_All.gpkg')
ts4_tup = (ts4_name, ts4_path, ts4_fs)
#
ts6_name = "Sep18"
ts6_path = os.path.join(home, "chm6.tif")
# ts6_fs = os.path.join(shps_root, 'FS_1718.shp')
ts6_fs = os.path.join(shps_root, 'Feed_signs_All.gpkg')
ts6_tup = (ts6_name, ts6_path, ts6_fs)

# #
# ts2_name = "Feb17"
# ts2_path = os.path.join(home, "chm2.tif")
# # ts2_fs = os.path.join(shps_root, 'FS_1617.shp')
# ts2_fs = os.path.join(shps_root, 'Feed_signs_All.gpkg')
# ts2_tup = (ts2_name, ts2_path, ts2_fs)
#
#
# #
# ts5_name = "Mar18"
# ts5_path = os.path.join(home, "chm5.tif")
# # ts5_fs = os.path.join(shps_root, 'FS_1718.shp')
# ts5_fs =os.path.join(shps_root, 'Feed_signs_All.gpkg')
# ts5_tup = (ts5_name, ts5_path, ts5_fs)

# all_feeding = os.path.abspath("C:/HG_Projects/CWC_Drone_work/shp_files/CWC_FS_clip1.shp")

CWC_CanChange_df = os.path.abspath("C:/HG_Projects/CWC_Drone_work/CWC_Results_Analysis/data/CWC_can_heights_df.csv")

def main():
    print("running Zone-Compare pipline")

    df_list = []
    name_list = []
    for ts in [ts1_tup,  ts3_tup, ts4_tup,  ts6_tup]: # ts5_tup,ts2_tup,
        change_ras = ts[1]
        name = ts[0]
        fs = ts[2]

        out_df = compare_zones(shrub_zones, change_ras, fs, name)

        df_list.append(out_df)
        name_list.append(name)

    tupzip = zip(name_list, df_list)

    for i in tupzip:
        df = i[1]
        name = i[0]
        df['time_step'] = name

    dfconcat = pd.concat(df_list)

    dfconcat.to_csv(CWC_CanChange_df, na_rep='NA', index=False)

def compare_zones(zones, diff_ras, feed_signs, name):
    z_gdf = gpd.read_file(zones)
    z_gdf['id'] = z_gdf.index
    z_gdf['area'] = z_gdf.area

    fs_gdf = gpd.read_file(feed_signs)
    fs_gdf_buff = fs_gdf.copy()
    fs_gdf_buff.geometry = fs_gdf.geometry.buffer(10)

    beav_zone = get_union(z_gdf, fs_gdf_buff)
    beav_zone['signs_YN'] = 1
    beav_zone = beav_zone.dissolve(by='signs_YN')
    beav_zone = beav_zone.reset_index()

    no_beav_zone = gpd.overlay(z_gdf, beav_zone, how='symmetric_difference')
    no_beav_zone['signs_YN'] = 0
    no_beav_zone = no_beav_zone.dissolve(by='signs_YN')
    no_beav_zone = no_beav_zone.reset_index()

    zones_comb = pd.concat([beav_zone, no_beav_zone])

    ax = zones_comb.plot(column='signs_YN', colormap='Dark2', edgecolor='None')
    fs_gdf.plot(color='red', ax=ax, markersize=8, alpha=0.3)
    plt.title(name + ' beaver and non beaver zones')
    plt.show()


    beav_zone_df = mask_ras_get_df(beav_zone, diff_ras)
    beav_zone_df['signs_YN'] = 1
    beav_zone_df['signs_YNf'] = 'Foraging Observed'

    no_beav_zone_df = mask_ras_get_df(no_beav_zone, diff_ras)
    no_beav_zone_df['signs_YN'] = 0
    no_beav_zone_df['signs_YNf'] = 'No Foraging'

    cwc_df = pd.concat([beav_zone_df, no_beav_zone_df])


    cwc_df['signs_YNf'] = cwc_df['signs_YNf'].astype('category')
    cwc_df['signs_YNf'] = cwc_df['signs_YNf'].cat.reorder_categories(['No Foraging', 'Foraging Observed'])

    sns.set(style="whitegrid")

    ax = sns.boxplot(x="signs_YNf", y="canopy_height", data=cwc_df, showfliers=True)
    # ax = sns.swarmplot(x="signs_YNf", y="canopy_height", data=cwc_df, color=".25")
    plt.title(name + ': canopy height')
    plt.ylabel('Canopy height (m)')
    plt.xlabel('')
    plt.show()

    np.random.seed(12345678)
    no_beav_group_o = no_beav_zone_df.canopy_height.dropna()
    yes_beav_group_o = beav_zone_df.canopy_height.dropna()

    t, p = stats.ttest_ind(no_beav_group_o, yes_beav_group_o, equal_var=False)
    omean1 = np.mean(no_beav_group_o)
    ostd1 = np.std(no_beav_group_o)
    ocount1 = len(no_beav_group_o)
    omean2 = np.mean(yes_beav_group_o)
    ostd2 = np.std(yes_beav_group_o)
    ocount2 = len(yes_beav_group_o)

    ocohens_d = effect_size_cohensD(omean1, ostd1, ocount1, omean2, ostd2, ocount2)
    print("-------- {0}: T test result ------------".format(name + 'Overall change'))
    print('                t value: {:.6f}'.format(t))
    print('                p value: {:.6f}'.format(p))
    print('                d value: {:.6f}'.format(ocohens_d))
    print("--------------------------------------------------------------")

    # -------------------------------------------------------------------- #
    z_df = z_gdf.drop(columns=['geometry'])

    return cwc_df

def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features'][0]['geometry']]


def effect_size_cohensD(mean1, std1, count1, mean2, std2, count2): #  https://gist.github.com/sriisking/10716f107ed30f4911639d695e3fbe49
    cohens_d = (mean1 - mean2) / np.sqrt(
        ((count1 - 1) * std1 ** 2 + (count2 - 1) * std2 ** 2) / (count1 + count2 - 2))
    return cohens_d

def get_union(to_select, union_shp):
    """Alternative union function"""
    union = gpd.GeoDataFrame(
        gpd.GeoSeries([union_shp.unary_union]),
        columns=['geometry'],
        crs=union_shp.crs)

    clip_gdf = gpd.sjoin(to_select, union, op='intersects')
    clip_gdf = clip_gdf.drop(columns=['index_right'])

    return clip_gdf

def mask_ras_get_df(gdf, ras):

        geom = getFeatures(gdf)  # returns geometries for AOIs

        with rasterio.open(ras) as src:
            out_image, out_transform = rasterio.mask.mask(src, geom, crop=True)

            res = src.res[0] * src.res[1]

        out_image_1d = out_image[0].flatten()
        out_image_1d = out_image_1d[out_image_1d != -999]

        # print('mean: {0}'.format(np.nanmean(out_image_1d)))
        # print('min: {0}'.format(np.nanmin(out_image_1d)))
        # print('max: {0}'.format(np.nanmax(out_image_1d)))
        # print('stdev: {0}'.format(np.nanstd(out_image_1d)))
        # print('sum: {0}'.format(np.nansum(out_image_1d)))
        # print('sum of canopy vol gain: {0}'.format(np.nansum(out_image_1d[out_image_1d > 0])* res))
        # print('sum of canopy volloss: {0}'.format(np.nansum(out_image_1d[out_image_1d < 0])* res))
        #
        # print('canopy vol gain/m^2: {0}'.format(np.nansum(out_image_1d[out_image_1d > 0]) / len(
        #     out_image_1d[out_image_1d > 0]) * res))
        #
        # print('canopy vol loss/m^2: {0}'.format(np.nansum(out_image_1d[out_image_1d < 0]) / len(
        #     out_image_1d[out_image_1d < 0]) * res))
        #
        # print('zone area (m^2): {0}'.format(gdf.area))

        out_df = pd.DataFrame(out_image_1d, columns=['canopy_height'])
        return out_df

if __name__ == '__main__':
    main()
