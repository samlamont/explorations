import os
import time

import pandas as pd
import numpy as np
import geopandas as gpd
from geopandas.tools import sjoin

gpd.options.use_pygeos = True  # huge gains!

SCS_CURVE_JAPAN = {
    0: ("unclassified", 100),
    1: ("water", 100),
    2: ("urban and built-up", 98),
    3: ("rice paddy", 94),
    4: ("crops", 89),
    5: ("Grassland", 78),
    6: ("deciduous broadleaf forest", 79),
    7: ("deciduous needle-leaved forest", 79),
    8: ("Evergreen broadleaf forest", 79),
    9: ("Evergreen needle-leaved forest", 79),
    10: ("bare land", 89),
    11: ("snow and ice", 100),
    253: ("other", 100),
    255: ("no data", 100),
}

NUM_SCHISM_GRIDS = 11
PATH_TO_CELL_DIMS = "/sam/urban/mesh_level/cell_dims_v1.16.1"
PATH_TO_CITY_POLYS = "/sam/gis/cities_107_with_grid_num.geojson"
PATH_T0_MESH_LANDCOVER = "/sam/urban/mesh_level/grid_landcover"
PATH_TO_BOUNDARY = "/sam/gis/polbnda_jpn_diss_4326.geojson"
OUTPUT_DIR = "/sam/urban/mesh_level/grid_maps"


def get_high_res_mask(
    df_cells: pd.DataFrame,
    gdf_geos: gpd.GeoDataFrame,
    grid_num: int,
    path_to_japan_boundary: str,
) -> pd.DataFrame:
    """
    This function produces the final mapping table for JP-wide SCHISM grid cells to use in the
    urban model

    Args:
        df_cells : Cell dimensions dataframe for a SCHISM grid with initial abstraction added
        gdf_geos : Geodataframe of 107 city polygons with ID and grid attributes
        grid_num : The SCHISM grid number (1-11)
        path_to_japan_boundary : Directory path to JP land mask polygon

    Returns: Dataframe of the final mapping table

    """

    gdf_this_grid = gdf_geos[gdf_geos.subregion == f"sub{grid_num}"].copy()
    gdf_jp_bnds = gpd.read_file(path_to_japan_boundary)
    gdf_pts = gpd.GeoDataFrame(
        df_cells,
        geometry=gpd.points_from_xy(df_cells.centroid_lon, df_cells.centroid_lat),
        crs="EPSG:4326",
    )
    # First create spatial index then join
    gdf_pts.sindex
    gdf_sj = sjoin(
        gdf_pts,
        gdf_this_grid[["geometry", "ID"]],
        how="left",
        predicate="within",  # faster than "intersects"?
    )
    gdf_sj["grid_number"] = grid_num
    gdf_sj.drop(columns="index_right", inplace=True)

    # Limit points to those only on JP mainland (not out in the sea)
    gdf_sj = sjoin(
        gdf_sj, gdf_jp_bnds[["geometry", "all"]], how="left", predicate="within"
    )
    gdf_sj = gdf_sj[gdf_sj["all"] == 1]

    gdf_sj.fillna(value=0, inplace=True)
    gdf_sj.drop(columns=["all", "index_right"], inplace=True)
    gdf_sj["resolution"] = 0
    gdf_sj.loc[gdf_sj.ID > 0, "resolution"] = 1
    gdf_sj["resolution"] = gdf_sj.resolution.astype(int)
    gdf_sj["ID"] = gdf_sj.ID.astype(int)
    gdf_sj["grid_number"] = gdf_sj.grid_number.astype(int)

    df_map = pd.DataFrame(gdf_sj.drop(columns=["geometry"]))
    df_map.rename(
        columns={
            "centroid_lat": "y",
            "centroid_lon": "x",
            "ID": "city_id",
            "cell_id": "schism_cell_id",
        },
        inplace=True,
    )

    return df_map


def precalculate_initial_abstraction(
    grid_num: int, path_to_mesh_landcover: str
) -> np.array:
    """
    This function pre-calculates initial abstraction for use in the urban model
    based on land cover and SCS_CURVE_JAPAN

    Args:
        grid_num : The SCHISM grid number
        path_to_mesh_landcover : Directory path to the JP land cover raster

    Returns: Initial abstraction as a numpy array

    """
    # Get the pre-calculated land cover array
    arr_lulc = np.load(f"{path_to_mesh_landcover}/grid_{grid_num}_landcover.npy")[:, 1]

    # Reclassify land cover to the CN
    for lc_class in SCS_CURVE_JAPAN:
        arr_lulc[arr_lulc == lc_class] = SCS_CURVE_JAPAN[lc_class][1]

    # Then calculate initial abstraction  (arr_lulc now holds CN values)
    max_sm = 25.4 * (1000.0 / arr_lulc - 10)  # metric conversion

    return max_sm


def build_mapping_table(
    path_to_cell_dims: str,
    path_to_city_polys: str,
    path_to_mesh_landcover: str,
    path_to_japan_boundary: str,
    output_dir: str,
):
    """
    The main function to build the urban-schism mapping table for running the pluvial
    model at the mesh-level. It first calculates initial abstraction based on summaries
    of land cover to each mesh cell (pre-calculated using a separate script)
    then adds city and grid IDs, and limits the points to those that intersect land

    Args:
        path_to_cell_dims : Path to the cell dimensions files
        path_to_city_polys : Path to city polygons layer
        path_to_mesh_landcover : Path to pre-calculated land cover summarized to mesh cells
        path_to_japan_boundary : Path to JP land boundary
        output_dir : User-specified output directory

    Returns: None. Saves individual grid tables and concatenated table to disk in output_dir

    """

    gdf_geos = gpd.read_file(path_to_city_polys)

    lst_output_paths = []
    for grid_num in range(1, NUM_SCHISM_GRIDS + 1):

        # Pre-calculate initial abstraction based on land cover
        max_sm = precalculate_initial_abstraction(grid_num, path_to_mesh_landcover)

        # Get a high-res/low-res mask, limit to JP land, and reformat columns
        df_cells = pd.read_csv(
            f"{path_to_cell_dims}/grid_{grid_num}_cell_dimensions.csv"
        )
        df_cells["initial_abstraction"] = max_sm
        df_map = get_high_res_mask(df_cells, gdf_geos, grid_num, path_to_japan_boundary)
        df_map.to_parquet(f"{output_dir}/urban_grid_map_{grid_num}.parquet")
        lst_output_paths.append(f"{output_dir}/urban_grid_map_{grid_num}.parquet")

        print(f"\tDone with grid {grid_num}, num cells: {len(df_map.index)}")

    # Concat all and format dtypes
    df_map_total = pd.concat(
        [pd.read_parquet(filepath) for filepath in lst_output_paths]
    )

    df_map_total["area"] = df_map_total.area.astype(np.float32)
    df_map_total["x"] = df_map_total.x.astype(np.float32)
    df_map_total["y"] = df_map_total.y.astype(np.float32)
    df_map_total["initial_abstraction"] = df_map_total.initial_abstraction.astype(
        np.float32
    )
    df_map_total["grid_number"] = df_map_total.grid_number.astype(np.uint)
    df_map_total["schism_cell_id"] = df_map_total.schism_cell_id.astype(np.uint)
    df_map_total["city_id"] = df_map_total.city_id.astype(np.uint)
    df_map_total["resolution"] = df_map_total.resolution.astype(np.uint8)

    df_map_total.to_parquet(f"{output_dir}/urban_grid_map.parquet")


if __name__ == "__main__":

    t0 = time.time()

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    build_mapping_table(
        PATH_TO_CELL_DIMS,
        PATH_TO_CITY_POLYS,
        PATH_T0_MESH_LANDCOVER,
        PATH_TO_BOUNDARY,
        OUTPUT_DIR,
    )
    print(f"Success! Total elapsed time: {(time.time() - t0) / 60:.2f} mins")
