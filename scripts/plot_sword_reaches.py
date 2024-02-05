import os
from pathlib import Path

import folium
from folium import IFrame
import geopandas as gpd
import matplotlib.pyplot as plt
import base64
import pandas as pd
import glob

plt.style.use("ggplot")


# def create_folium_map():
#     print("Getting reaches")
#     gdf = gpd.read_parquet("/mnt/data/af_sword_reaches_hb13_v14.parquet")

#     xmin, ymin, xmax, ymax = gdf.total_bounds

#     x = (xmax + xmin) / 2
#     y = (ymax + ymin) / 2

#     print("Building map...")

#     # Popup setup...
#     resolution = 80  # originally set at 75 (try lower for size??)
#     width = 6.5
#     height = 4.5

#     #    str_map_title='GagesII-Ref above minor flood stage'
#     str_map_title = "Streamflow Forecasts"

#     m = folium.Map([y, x], zoom_start=4, tiles="cartodbpositron")

#     folium.GeoJson(gdf.geometry).add_to(m)

#     m.save("/mnt/data/test_map.html")


def embed_plots():
    reaches_gdf = gpd.read_parquet(SWORD_REACHES)
    xmin, ymin, xmax, ymax = reaches_gdf.total_bounds
    x = (xmax + xmin) / 2
    y = (ymax + ymin) / 2

    print("Building map...")

    # Popup setup...
    resolution = 80  # originally set at 75 (try lower for size??)
    width = 6.5
    height = 4.5

    #    str_map_title='GagesII-Ref above minor flood stage'
    str_map_title = "Streamflow Forecasts"

    m = folium.Map([y, x], zoom_start=8, tiles=None)
    
    tile = folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Esri',
        name='Esri Satellite',
        overlay=False,
        control=True
       ).add_to(m)

    # reaches_gdf = gpd.read_parquet("/mnt/data/af_sword_reaches_hb13_v14.parquet")
    reaches_gdf["reach_id"] = reaches_gdf.reach_id.astype(str)
    reaches_gdf.set_index("reach_id", inplace=True)

    # Get list of filenames and paths in selected directory
    # Path(DATA_PATH, "plots"
    lst_filenames = glob.glob("/mnt/data/sword/sacramento_2024/plots/" + "*.png")
    filepaths = []
    for this_path in lst_filenames:
        # Get the plot filename and reach id
        str_file = os.path.basename(this_path)
        reach_id = os.path.splitext(str_file)[0]
        filepaths.append({"filepath": this_path, "reach_id": reach_id})
    df_files = pd.DataFrame(filepaths)
    df_files["reach_id"] = df_files.reach_id.astype(str)
    df_files.set_index("reach_id", inplace=True)

    cntr = 0
    fg_markers = folium.FeatureGroup(
        name="Reaches"
    )  # Add markers to this layer
    fg_popups = folium.FeatureGroup(name="Plots")
    for tpl in reaches_gdf.itertuples():
        geom = tpl.geometry
        reach_id = tpl.Index

        try:
            this_path = df_files.loc[reach_id].filepath

            print("yes!")

            encoded = base64.b64encode(open(this_path, "rb").read())
            html = '<img src="data:image/png;base64,{}">'.format
            iframe = IFrame(
                html(encoded.decode("UTF-8")),
                width=(width * resolution) + 75,
                height=(height * resolution) + 50,
            )

            popup = folium.Popup(iframe, max_width=2650)

            b = folium.GeoJson(geom, popup=popup)
            fg_popups.add_child(b)
            # fg_markers.add_child(b)
        except:
            print("skipping!")

        b2 = folium.GeoJson(geom)
        fg_markers.add_child(b2)

        # b.add_child(popup)

        cntr += 1

        # if cntr > 100:
        #     break

    m.add_child(fg_markers)
    m.add_child(fg_popups)
    m.add_child(folium.LayerControl())
    m.save(Path(DATA_PATH, "sacramento_map.html"))


def build_one_ensemble_plot(reach_id, q_arr, dt_rng):
    fig, ax = plt.subplots()

    # if not any(q_arr):
    #     return
    if q_arr[q_arr > 0].size == 0:
        return

    ax.plot(
        dt_rng,
        q_arr.T,
        linewidth=1.0,
        # color="blue",
        linestyle="-",
        label=reach_id,
    )
    
    # Ensemble mean
    ax.plot(
        dt_rng,
        q_arr.T.mean(axis=1),
        linewidth=1.0,
        color="black",
        linestyle="--",
    )    

    # Formatting and annotations...
    ax.set_title(
        f"Reach {reach_id}",
        fontsize=10,
    )

    fig.autofmt_xdate()

    ax.set_ylabel("Q (m^3/s)")

    # Save plots as png...
    # figname = f"/mnt/data/plots/{tpl.Index}.png" 
    figname = Path(DATA_PATH, "plots", f"{reach_id}.png")
    fig.savefig(figname, dpi=90, bbox_inches="tight")

    fig.clf()
    plt.close(fig)
    plt.close("all")
    
    
def build_one_plot(i, tpl, vals, dt_rng):
    fig, ax = plt.subplots()

    if not any(vals):
        return
    if vals[vals > 0].size == 0:
        return

    ax.plot(
        dt_rng,
        vals,
        linewidth=3.0,
        color="blue",
        linestyle="-",
        label=tpl.Index,
    )

    # Formatting and annotations...
    ax.set_title(
        f"Reach {tpl.Index}",
        fontsize=10,
    )

    fig.autofmt_xdate()

    ax.set_ylabel("Q")

    # Save plots as png...
    # figname = f"/mnt/data/plots/{tpl.Index}.png" 
    figname = Path(DATA_PATH, "plots", f"{tpl.Index}.png")
    fig.savefig(figname, dpi=90, bbox_inches="tight")

    fig.clf()
    plt.close(fig)
    plt.close("all")


def create_plots():
    reaches_gdf = gpd.read_parquet(SWORD_REACHES_DATA)
    reaches_gdf.set_index("reach_id", inplace=True)
    start = pd.to_datetime("2024-02-03 00:00")  # 2024020300
    end = pd.to_datetime("2024-02-12 18:00")
    dt_rng = pd.date_range(start, end, freq="6H")
    dt_rng_str = dt_rng.strftime("%Y%m%d%H")
    # q_arr = reaches_gdf[dt_rng_str].values
    
    gps = reaches_gdf.groupby("reach_id")
    
    for reach_id, df in gps:
        q_arr = df[dt_rng_str].values
        build_one_ensemble_plot(reach_id, q_arr, dt_rng)
        # break

    # for i, tpl in enumerate(reaches_gdf.itertuples()):
    #     vals = q_arr[i, :]
    #     build_one_plot(i, tpl, vals, dt_rng)
        # break


if __name__ == "__main__":
    
    DATA_PATH = "/mnt/data/sword/sacramento_2024/"
    
    # POLYGONS = "/mnt/data/congo_discharge_polygons_repaired.parquet"
    # SWORD_REACHES = "/mnt/data/af_sword_reaches_hb13_v14.parquet"
    
    POLYGONS = Path(DATA_PATH, "sacramentovalley6hdata", "RRM_202402.shp")
    SWORD_REACHES = Path(DATA_PATH, "sacramentovalley6hdata", "reaches_clipped.parquet")   
    SWORD_NODES = Path(DATA_PATH, "sacramentovalley6hdata", "nodes_clipped.parquet")  
    
    SWORD_REACHES_DATA = Path(DATA_PATH, "sacramentovalley6hdata", "sword_reaches_with_discharge.parquet")
    
    # # Open files
    # polys_gdf = gpd.read_file(POLYGONS)
    # nodes_gdf = gpd.read_parquet(SWORD_NODES)
    # reaches_gdf = gpd.read_parquet(SWORD_REACHES)

    # # Join
    # polys_gdf["CATID"] = polys_gdf.CATID.astype(str)

    # # Find most-occurring node in each polygon
    # gdf_overlay = gpd.overlay(nodes_gdf, polys_gdf[["geometry", "CATID"]])
    # gp = gdf_overlay.groupby(by="reach_id")
    # reference = []
    # for reach_id, gdf in gp:
    #     most_commmon_catid = gdf.CATID.value_counts().idxmax()
    #     reference.append({"reach_id": reach_id, "catid": most_commmon_catid})
    # df_lookup = pd.DataFrame(reference)

    # # Join CATID to reaches
    # reaches_gdf["reach_id"] = reaches_gdf.reach_id.astype(str)
    # reaches_gdf.set_index("reach_id", inplace=True)

    # df_lookup["reach_id"] = df_lookup.reach_id.astype("str")
    # df_lookup.set_index("reach_id", inplace=True)

    # gdf_catid_lines = reaches_gdf.join(df_lookup)
    # gdf_catid_lines.reset_index(inplace=True)
    # gdf_catid_lines.dropna(subset=["catid"], inplace=True)

    # gdf_catid_lines["catid"] = gdf_catid_lines.catid.astype(str)
    # gdf_catid_lines.set_index("catid", inplace=True)

    # polys_gdf.set_index("CATID", inplace=True)

    # # Join streamflow data to reaches
    # reaches_gdf_data = gdf_catid_lines.join(polys_gdf, rsuffix="_r")
    # reaches_gdf_data.drop(columns=["geometry_r"], inplace=True)
    # # reaches_gdf_data.loc["121520"]
    # reaches_gdf_data.to_parquet(SWORD_REACHES_DATA)

    # Create individual popup plots and save as png for each reach
    # create_plots()

    # Build the final map
    embed_plots()

    pass
