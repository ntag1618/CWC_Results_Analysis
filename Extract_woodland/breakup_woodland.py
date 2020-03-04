# create vector grid to the extent of the DSM - 30-40m?
# delete any grid cells that don't intersect woodland vector.
# Use this modified grid as

import geopandas as gpd
from shapely.geometry import Polygon
import numpy as np
import os
from matplotlib import pyplot as plt
zones_out = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/Woodland_Zones50m.gpkg')

trees_vec = os.path.abspath('C:/HG_Projects/CWC_Drone_work/CHM/Rip_Vec_Sep17.gpkg')

def main():
    trees_gdf = gpd.read_file(trees_vec)

    xmin,ymin,xmax,ymax = trees_gdf.total_bounds

    length = 50
    wide = 50

    cols = list(range(int(np.floor(xmin)), int(np.ceil(xmax)), wide))
    rows = list(range(int(np.floor(ymin)), int(np.ceil(ymax)), length))
    rows.reverse()

    polygons = []
    for x in cols:
        for y in rows:
            polygons.append( Polygon([(x,y), (x+wide, y), (x+wide, y-length), (x, y-length)]) )

    grid = gpd.GeoDataFrame({'geometry':polygons})

    zones = gpd.overlay(grid, trees_gdf, how='intersection')
    zones.crs = trees_gdf.crs
    print(zones.crs)
    zones.plot(edgecolor='black')
    plt.show()
    zones.to_file(zones_out, driver='GPKG')

if __name__ == '__main__':
    main()