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

# suppress warnings for now...
import warnings
warnings.filterwarnings("ignore")

home = os.path.abspath("C:/HG_Projects/CWC_Drone_work/Prec_Anal_Exports/Rasters_v4")
shps_root = os.path.abspath("C:/HG_Projects/CWC_Drone_work/shp_files")
shrub_zones = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/Woodland_Zones20m.gpkg')

# LATE SUMMER SEP-SEP
cr1_name = ("Sep17 - Sep18")
cr1_path = os.path.join(home, "DoD_Sep17_Sep18.tif")
cr1_fs = os.path.join(shps_root, 'FS_1618.shp')
cr1_tup = (cr1_name, cr1_path, cr1_fs)


# WINTERS
cr2_name = "Dec16 - Feb17"
cr2_path = os.path.join(home, "DOD_Dec16_Feb17.tif")
cr2_fs = os.path.join(shps_root, 'FS_1617.shp')
cr2_tup = (cr2_name, cr2_path, cr2_fs)

cr3_name = "Feb17 - Jan18"
cr3_path = os.path.join(home, "DOD_Feb17_Jan18.tif")
cr3_fs = os.path.join(shps_root, 'FS_1718.shp')
cr3_tup = (cr3_name, cr3_path, cr3_fs)

cr4_name = "Jan18 - Mar18"
cr4_path = os.path.join(home, "DOD_Jan18_March18.tif")
cr4_fs = os.path.join(shps_root, 'FS_1718.shp')
cr4_tup = (cr4_name, cr4_path, cr4_fs)

# WINTER START-END
cr5_name = "Dec16 - Mar18"
cr5_path = os.path.join(home, "DoD_Dec16_Mar18.tif")
cr5_fs = os.path.join(shps_root, 'FS_1618.shp')
cr5_tup = (cr5_name, cr5_path, cr5_fs)


# all_feeding = os.path.abspath("C:/HG_Projects/CWC_Drone_work/shp_files/CWC_FS_clip1.shp")


def main():
    print("running Zone-Compare pipline")

    for ts in [cr1_tup, cr2_tup, cr3_tup, cr4_tup, cr5_tup]:
        change_ras = ts[1]
        name = ts[0]
        fs = ts[2]
        # geodf = ras_to_points(change_ras)

        # gdf_dist = dist_to_nearest_point(geodf, fs)

        # gdf = plotting(gdf_dist, name)

        compare_zones(shrub_zones, change_ras, fs, name)

def compare_zones(zones, diff_ras, feed_signs, name):
    z_gdf = gpd.read_file(zones)
    z_gdf['id'] = z_gdf.index
    z_gdf['area'] = z_gdf.area
    z_gdf['mean'] = np.nan
    z_gdf['min'] = np.nan
    z_gdf['max'] = np.nan
    z_gdf['stdev'] = np.nan
    z_gdf['sumall'] = np.nan
    z_gdf['sumpos'] = np.nan
    z_gdf['sumneg'] = np.nan
    z_gdf['n_signs'] = np.nan
    z_gdf['signs_YN'] = np.nan
    z_gdf['signs_YNf'] = np.nan

    fs_gdf = gpd.read_file(feed_signs)
    fs_gdf_buff = fs_gdf.copy()
    fs_gdf_buff.geometry = fs_gdf.geometry.buffer(10)

    for idx, row in z_gdf.iterrows():
        tempgdf = gpd.GeoDataFrame(
            gpd.GeoSeries(z_gdf.loc[idx, 'geometry']),
            columns=['geometry'])
        tempgdf.crs = z_gdf.crs

        geom = getFeatures(tempgdf)  # returns buffered geometries for AOIs

        with rasterio.open(diff_ras) as src:
            out_image, out_transform = rasterio.mask.mask(src, geom, crop=True)
            res = src.res[0] * src.res[1]
            out_meta = src.meta

        out_image_1d = out_image.flatten()
        out_image_1d = out_image_1d[out_image_1d != -999]
        if len(out_image_1d) < 1:
            z_gdf.loc[idx, 'mean'] = np.nan
            z_gdf.loc[idx, 'min'] = np.nan
            z_gdf.loc[idx, 'max'] = np.nan
            z_gdf.loc[idx, 'stdev'] = np.nan
            z_gdf.loc[idx, 'sumall'] = np.nan
            z_gdf.loc[idx, 'sumpos'] = np.nan
            z_gdf.loc[idx, 'sumneg'] = np.nan
            z_gdf.loc[idx, 'percpos'] = np.nan
            z_gdf.loc[idx, 'percneg'] = np.nan
        else:
            z_gdf.loc[idx, 'mean'] = np.nanmean(out_image_1d)
            z_gdf.loc[idx, 'min'] = np.nanmin(out_image_1d)
            z_gdf.loc[idx, 'max'] = np.nanmax(out_image_1d)
            z_gdf.loc[idx, 'stdev'] = np.nanstd(out_image_1d)
            z_gdf.loc[idx, 'sumall'] = np.nansum(out_image_1d)
            z_gdf.loc[idx, 'sumpos'] = np.nansum(out_image_1d[out_image_1d > 1])
            z_gdf.loc[idx, 'sumneg'] = np.nansum(out_image_1d[out_image_1d < -1])
            z_gdf.loc[idx, 'percpos'] = np.nansum(out_image_1d[out_image_1d > 0]) / len(
                out_image_1d[out_image_1d > 0]) * res
            z_gdf.loc[idx, 'percneg'] = np.nansum(out_image_1d[out_image_1d < 0]) / len(
                out_image_1d[out_image_1d < 0]) * res

        union_fs = gpd.overlay(tempgdf, fs_gdf_buff, how='intersection')
        z_gdf.loc[idx, 'n_signs'] = len(union_fs)
        if len(union_fs) == 0:
            z_gdf.loc[idx, 'signs_YN'] = 0
            z_gdf.loc[idx, 'signs_YNf'] = 'No Foraging'
        else:
            z_gdf.loc[idx, 'signs_YN'] = 1
            z_gdf.loc[idx, 'signs_YNf'] = 'Foraging Observed'

    z_gdf['signs_YNf'] = z_gdf['signs_YNf'].astype('category')
    z_gdf['signs_YNf'] = z_gdf['signs_YNf'].cat.reorder_categories(['No Foraging', 'Foraging Observed'])

    sns.set(style="whitegrid")

    ax = z_gdf.plot(column='signs_YNf', colormap='Dark2', edgecolor='None')
    fs_gdf.plot(color='red', ax=ax, markersize=8, alpha=0.3)
    plt.title(name + ' beaver and non beaver zones')
    plt.show()

    ax = sns.boxplot(x="signs_YNf", y="percneg", data=z_gdf, showfliers=False, )
    ax = sns.swarmplot(x="signs_YNf", y="percneg", data=z_gdf, color=".25")
    plt.title(name + ': canopy change in areas of loss')
    plt.ylabel('Canopy height loss/m$^2$')
    plt.xlabel('')
    plt.show()

    ax = sns.boxplot(x="signs_YNf", y="percpos", data=z_gdf, showfliers=False, )
    ax = sns.swarmplot(x="signs_YNf", y="percpos", data=z_gdf, color=".25")
    plt.title(name + ': canopy change in areas of growth')
    plt.ylabel('Canopy height gain/m$^2$')
    plt.xlabel('')
    plt.show()

    ax = sns.boxplot(x="signs_YNf", y="mean", data=z_gdf, showfliers=False, )
    ax = sns.swarmplot(x="signs_YNf", y="mean", data=z_gdf, color=".25")
    plt.title(name + ': overall canopy change per m$^2$')
    plt.ylabel('Canopy height cahnge/m$^2$')
    plt.xlabel('')
    plt.show()

    kw = dict(column='percneg', k=100, edgecolor='black', colormap='coolwarm_r')
    ax = z_gdf.plot(scheme='quantiles', **kw)
    fs_gdf.plot(color='black', ax=ax, markersize=8, alpha=0.5)
    plt.title(name + ' feeding locations with canopy loss/m$^2$')
    plt.show()

    kw = dict(column='percpos', k=100, edgecolor='black', colormap='coolwarm')
    ax = z_gdf.plot(scheme='quantiles', **kw)
    fs_gdf.plot(color='black', ax=ax, markersize=8, alpha=0.5)
    plt.title(name + ' feeding locations with canopy gain/m$^2$')
    plt.show()

    np.random.seed(12345678)
    no_beav_group_n = z_gdf['percneg'][z_gdf['signs_YN'] == 0].dropna()
    yes_beav_group_n = z_gdf['percneg'][z_gdf['signs_YN'] == 1].dropna()

    t, p = stats.ttest_ind(no_beav_group_n, yes_beav_group_n, equal_var=False)

    mean1 = np.mean(no_beav_group_n)
    std1 = np.std(no_beav_group_n)
    count1 = len(no_beav_group_n)
    mean2 = np.mean(yes_beav_group_n)
    std2 = np.std(yes_beav_group_n)
    count2 = len(yes_beav_group_n)

    cohens_d = effect_size_cohensD(mean1, std1, count1, mean2, std2, count2)

    print("-------- {0}: T test result ------------".format(name + '_p_neg'))
    print('                t value: {:.6f}'.format(t))
    print('                p value: {:.6f}'.format(p))
    print('                d value: {:.6f}'.format(cohens_d))
    print("--------------------------------------------------------------")

    no_beav_group_p = z_gdf['percpos'][z_gdf['signs_YN'] == 0].dropna()
    yes_beav_group_p = z_gdf['percpos'][z_gdf['signs_YN'] == 1].dropna()

    t, p = stats.ttest_ind(no_beav_group_p, yes_beav_group_p, equal_var=False)

    pmean1 = np.mean(no_beav_group_p)
    pstd1 = np.std(no_beav_group_p)
    pcount1 = len(no_beav_group_p)
    pmean2 = np.mean(yes_beav_group_p)
    pstd2 = np.std(yes_beav_group_p)
    pcount2 = len(yes_beav_group_p)

    pcohens_d = effect_size_cohensD(pmean1, pstd1, pcount1, pmean2, pstd2, pcount2)
    print("-------- {0}: T test result ------------".format(name + '_p_pos'))
    print('                t value: {:.6f}'.format(t))
    print('                p value: {:.6f}'.format(p))
    print('                d value: {:.6f}'.format(pcohens_d))
    print("--------------------------------------------------------------")

    no_beav_group_o = z_gdf['mean'][z_gdf['signs_YN'] == 0].dropna()
    yes_beav_group_o = z_gdf['mean'][z_gdf['signs_YN'] == 1].dropna()

    t, p = stats.ttest_ind(no_beav_group_o, yes_beav_group_o, equal_var=False)

    omean1 = np.mean(no_beav_group_o)
    ostd1 = np.std(no_beav_group_o)
    ocount1 = len(no_beav_group_o)
    omean2 = np.mean(yes_beav_group_o)
    ostd2 = np.std(yes_beav_group_o)
    ocount2 = len(yes_beav_group_o)

    ocohens_d = effect_size_cohensD(omean1, ostd1, ocount1, omean2, ostd2, ocount2)
    print("-------- {0}: T test result ------------".format(name + 'Overall cahnge'))
    print('                t value: {:.6f}'.format(t))
    print('                p value: {:.6f}'.format(p))
    print('                d value: {:.6f}'.format(ocohens_d))
    print("--------------------------------------------------------------")

    # -------------------------------------------------------------------- #


def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    return [json.loads(gdf.to_json())['features'][0]['geometry']]


def effect_size_cohensD(mean1, std1, count1, mean2, std2, count2): #  https://gist.github.com/sriisking/10716f107ed30f4911639d695e3fbe49
    cohens_d = (mean1 - mean2) / np.sqrt(
        ((count1 - 1) * std1 ** 2 + (count2 - 1) * std2 ** 2) / (count1 + count2 - 2))
    return cohens_d


if __name__ == '__main__':
    main()
