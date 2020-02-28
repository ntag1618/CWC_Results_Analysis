# -------------------------------------------------------------------------------------------- #
# ------------------- CANOPY STRUCTURE CHANGE IN RESPONSE TO BEAVER FORAGING------------------ #
# ---- A workflow to generate statistics and visualisations for SFM point cloud data --------- #
# ---- from a beaver-impacted riparian area to understand changes in canopy structure -------- #
# -------------------------------------------------------------------------------------------- #

# Written by:

# Version:

# Notes:

# -------------------------------------------------------------------------------------------- #
import geopandas as gpd
import pandas as pd
import rasterio
import os
import numpy as np
from shapely.geometry import Point
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
import pathlib

# suppress warnings for now...
import warnings
warnings.filterwarnings("ignore")

home = os.path.abspath("C:/HG_Projects/CWC_Drone_work/Prec_Anal_Exports/Rasters_v3")
shps_root = os.path.abspath("C:/HG_Projects/CWC_Drone_work/shp_files")
shrub_zones = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/Woodland_Zones.gpkg')

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


def main():
    print("running distances to points pipline")

    for ts in [cr1_tup, cr2_tup, cr3_tup, cr4_tup, cr5_tup]:
        change_ras = ts[1]
        name = ts[0]
        fs = ts[2]
        geodf = ras_to_points(change_ras)

        gdf_dist = dist_to_nearest_point(geodf, fs)

        gdf = plotting(gdf_dist, name)


def ras_to_points(ras_path):
    print("converting raster to points.")
    with rasterio.open(ras_path) as src:

        print("retrieving raster attributes")
        arr = src.read(1)
        height, width = np.shape(arr)
        x_size, y_size = src.res
        left = src.bounds[0]
        top = src.bounds[3]
        bottom = src.bounds[1]
        right = src.bounds[2]
        nodata = src.nodata

        arr_flip = np.flip(arr, 0) # not sure why this is needed but solves the upaide down plotting issue.

        print("get raster cell geometries where data is not NULL")
        geom_list = []
        values = []
        for row in tqdm(range(height)):
            for col in range(width):
                if nodata != arr_flip[row, col]:
                    x = col * x_size + left + (x_size/2)
                    y = row * y_size + bottom + (y_size/2)

                    geom = Point(x, y)

                    geom_list.append(geom)
                    values.append(arr_flip[row, col])
                    # coords.append((x, y))

        print("convert geometries into geo data frame")
        gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(geom_list))
        gdf['canopy_change'] = values

    # gdf.plot(column='canopy_change')
    # plt.show()
    # print(len(gdf.index))
    return gdf


def dist_to_nearest_point(grid_points, feed_signs):
    print("calculating distance to closest feeding sign")

    findpoints = gpd.read_file(feed_signs)

    findpoints = findpoints[['geometry']]

    def ckdnearest(gdA, gdB):
        nA = np.array(list(zip(gdA.geometry.x, gdA.geometry.y)))
        nB = np.array(list(zip(gdB.geometry.x, gdB.geometry.y)))
        btree = cKDTree(nB)
        dist, idx = btree.query(nA, k=1)
        gdf = pd.concat(
            [gdA.reset_index(drop=True), gdB.loc[idx, gdB.columns != 'geometry'].reset_index(drop=True),
             pd.Series(dist, name='dist')], axis=1)
        return gdf

    out_gdf = ckdnearest(grid_points, findpoints)

    return out_gdf


def plotting(gdf, name):
    plot_folder = os.path.join(pathlib.Path(__file__).parent.absolute(), 'plots')
    if os.path.isdir(plot_folder):
        pass
    else:
        os.mkdir(plot_folder)
    # ----------------------------Multi Category density facet plot------------------------------------ #
    lossm5 = gdf['dist'][gdf['canopy_change'] < -5].copy()
    lossm2m5 = gdf['dist'][(gdf['canopy_change'] > -5) & (gdf['canopy_change'] < -2)].copy()
    lossm2z = gdf['dist'][(gdf['canopy_change'] > -2) & (gdf['canopy_change'] < 0)].copy()
    zero = gdf['dist'][(gdf['canopy_change'] == 0)].copy()
    gainz2 = gdf['dist'][(gdf['canopy_change'] > 0) & (gdf['canopy_change'] < 2)].copy()
    gain2_5 = gdf['dist'][(gdf['canopy_change'] > 2) & (gdf['canopy_change'] < 5)].copy()
    gain5 = gdf['dist'][gdf['canopy_change'] > 5].copy()

    # plot
    fig, axes = plt.subplots(1, 7, figsize=(12, 7), sharey=True)
    sns.distplot(lossm5, color="dodgerblue", ax=axes[0], axlabel='<-5 m', bins=10)
    sns.distplot(lossm2m5, color="deeppink", ax=axes[1], axlabel='-5:-2 m', bins=10)
    sns.distplot(lossm2z, color="gold", ax=axes[2], axlabel='-2:0 m', bins=10)
    sns.distplot(zero, color="grey", ax=axes[3], axlabel='-2:0 m', bins=10)
    sns.distplot(gainz2, color="green", ax=axes[4], axlabel='0:2 m', bins=10)
    sns.distplot(gain2_5, color="darkorange", ax=axes[5], axlabel='2:5 m', bins=10)
    sns.distplot(gain5, color="indigo", ax=axes[6], axlabel='> 5 m', bins=10)

    axes[3].set_title('{0}: Distribution of distance to nearest feeding sign across change categories'.format(name))
    axes[0].set_ylabel('Density')
    for i in range(7):
        axes[i].set_xlim(-10, 150)

    plt.show()

    dens_dist_plot = os.path.join(plot_folder, "{}_dist_density_facet.png".format(name))
    fig.savefig(fname=dens_dist_plot, dpi=300, format='png')

    # ----------------------------scatter plot dist vs canopy vol change---------------------------------- #

    gdf['direc'] = 'None'
    gdf.loc[gdf['canopy_change'] > 0, 'direc'] = 'gain'
    gdf.loc[gdf['canopy_change'] < 0, 'direc'] = 'loss'

    change_df = gdf[(gdf['direc'] == 'gain') | (gdf['direc'] == 'loss')].copy()
    change_df['volume change'] = change_df['canopy_change'] * 0.5

    fig, axes = plt.subplots(figsize=(10, 6))

    newPal = dict(gain='royalblue', loss='orangered')
    sns.scatterplot(x="dist", y="volume change", hue="direc",
                         data=change_df, alpha=0.3, linewidth=0, palette=newPal, s=1.5)

    plt.xlabel("Distance from nearest feeding sign in this period (m)")
    plt.ylabel("Canopy volume change per cell (m$^3$)")
    axes.set_axisbelow(False)
    plt.grid(linestyle='--')
    axes.get_legend().remove()
    plt.title(name)

    plt.ylim(-3, 3)
    plt.show()

    dist_change_plot = os.path.join(plot_folder, "{}_dist_change.png".format(name))
    fig.savefig(fname=dist_change_plot, dpi=300, format='png')

    # --------------------------stem plot - summed canopy volume for binned distance-------------------- #
    group_df = gdf.copy()
    tot_range = np.linspace(group_df['dist'].min(), group_df['dist'].max(), 250)

    group_df['direc'] = 'None'
    group_df.loc[group_df['canopy_change'] > 0, 'direc'] = 'gain'
    group_df.loc[group_df['canopy_change'] < 0, 'direc'] = 'loss'

    gains_df = group_df[group_df['canopy_change'] > 0].copy()
    loss_df = group_df[group_df['canopy_change'] < 0].copy()

    def grouping_func(df):
        df = df.drop(columns=['geometry'])
        df = df.sort_values(by='dist')
        df['groups'] = pd.cut(df.dist, tot_range)
        df['categories'] = (df.groups != df.groups.shift()).cumsum()

        df_a = df.drop(columns=['dist'])
        df_b = df.drop(columns=['canopy_change'])

        df_a = df_a.groupby(by='categories').sum()
        df_b = df_b.groupby(by='categories').mean()

        df = df_a.copy()
        df['dist'] = df_b.dist

        df['volume change'] = df['canopy_change'] * 0.5

        return df

    gains_df = grouping_func(gains_df)
    gains_df['direc'] = 'gain'
    loss_df = grouping_func(loss_df)
    loss_df['direc'] = 'loss'

    fig, axes = plt.subplots(1, 1, figsize=(10, 6), sharey=True)

    plt.stem(gains_df['dist'], gains_df['volume change'], 'royalblue',
             markerfmt='royalblue', use_line_collection=True)
    plt.stem(loss_df['dist'], loss_df['volume change'], 'orangered',
             markerfmt='orangered', use_line_collection=True)

    plt.xlabel("Distance from nearest feeding sign in this period (m)")
    plt.ylabel("Canopy volume change (m$^3$)")
    axes.set_axisbelow(False)
    plt.grid(linestyle='--')
    plt.title(name)
    # axes.set_xticklabels(0, 150)

    plt.show()

    binned_dist_change = os.path.join(plot_folder, "{}_binned_dist_change.png".format(name))
    fig.savefig(fname=binned_dist_change, dpi=300, format='png')

    return gdf

if __name__ == '__main__':
    main()


