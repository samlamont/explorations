import os

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
    gdf = gpd.read_parquet("/mnt/data/af_sword_reaches_hb13_v14.parquet")
    xmin, ymin, xmax, ymax = gdf.total_bounds
    x = (xmax + xmin) / 2
    y = (ymax + ymin) / 2

    print("Building map...")

    # Popup setup...
    resolution = 80  # originally set at 75 (try lower for size??)
    width = 6.5
    height = 4.5

    #    str_map_title='GagesII-Ref above minor flood stage'
    str_map_title = "Streamflow Forecasts"

    m = folium.Map([y, x], zoom_start=4, tiles="cartodbpositron")

    gdf_lines = gpd.read_parquet("/mnt/data/af_sword_reaches_hb13_v14.parquet")
    gdf_lines["reach_id"] = gdf_lines.reach_id.astype(str)
    gdf_lines.set_index("reach_id", inplace=True)

    # Get list of filenames and paths in selected directory
    lst_filenames = glob.glob("/mnt/data/plots/" + "*.png")
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
    for tpl in gdf_lines.itertuples():
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
    m.save("/mnt/data/test_map.html")


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
    figname = f"/mnt/data/plots/{tpl.Index}.png"
    fig.savefig(figname, dpi=90, bbox_inches="tight")

    fig.clf()
    plt.close(fig)
    plt.close("all")


def create_plots():
    gdf_lines = gpd.read_parquet(
        "/mnt/data/sword_reaches_with_discharge.parquet"
    )
    gdf_lines.set_index("reach_id", inplace=True)
    start = pd.to_datetime("2021-01-01")
    end = pd.to_datetime("2022-11-30")
    dt_rng = pd.date_range(start, end)
    dt_rng_str = dt_rng.strftime("%Y-%m-%d")
    q_arr = gdf_lines[dt_rng_str].values

    for i, tpl in enumerate(gdf_lines.itertuples()):
        vals = q_arr[i, :]
        build_one_plot(i, tpl, vals, dt_rng)
        # break


if __name__ == "__main__":
    POLYGONS = "/mnt/data/congo_discharge_polygons_repaired.parquet"
    SWORD_REACHES = "/mnt/data/af_sword_reaches_hb13_v14.parquet"

    # Join
    gdf_polys = gpd.read_parquet(
        "/mnt/data/congo_discharge_polygons_repaired.parquet"
    )
    gdf_polys["CATID"] = gdf_polys.CATID.astype(str)

    # # Find most-occurring node in each polygon
    # gdf_nodes = gpd.read_parquet("/mnt/data/af_sword_nodes_hb13_v14.parquet")
    # gdf_overlay = gpd.overlay(gdf_nodes, gdf_polys[["geometry", "CATID"]])
    # gp = gdf_overlay.groupby(by="reach_id")
    # reference = []
    # for reach_id, gdf in gp:
    #     most_commmon_catid = gdf.CATID.value_counts().idxmax()
    #     reference.append({"reach_id": reach_id, "catid": most_commmon_catid})
    # df_lookup = pd.DataFrame(reference)

    # # Join CATID to reaches
    # gdf_lines = gpd.read_parquet("/mnt/data/af_sword_reaches_hb13_v14.parquet")
    # gdf_lines["reach_id"] = gdf_lines.reach_id.astype(str)
    # gdf_lines.set_index("reach_id", inplace=True)

    # df_lookup["reach_id"] = df_lookup.reach_id.astype("str")
    # df_lookup.set_index("reach_id", inplace=True)

    # gdf_catid_lines = gdf_lines.join(df_lookup)
    # gdf_catid_lines.reset_index(inplace=True)
    # gdf_catid_lines.dropna(subset=["catid"], inplace=True)

    # gdf_catid_lines["catid"] = gdf_catid_lines.catid.astype(str)
    # gdf_catid_lines.set_index("catid", inplace=True)

    # gdf_polys.set_index("CATID", inplace=True)

    # # Join streamflow data to reaches
    # gdf_lines_data = gdf_catid_lines.join(gdf_polys, rsuffix="_r")
    # gdf_lines_data.drop(columns=["geometry_r"], inplace=True)
    # # gdf_lines_data.loc["121520"]
    # gdf_lines_data.to_parquet("/mnt/data/sword_reaches_with_discharge.parquet")

    # Create individual popup plots and save as png for each reach
    # create_plots()

    # Build the final map
    embed_plots()

    pass
