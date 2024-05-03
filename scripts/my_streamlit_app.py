# my_app.py
import folium
import folium.vector_layers
from streamlit_folium import st_folium
import streamlit as st
import geopandas as gpd
from shapely.geometry import Point
from folium import IFrame
import matplotlib.pyplot as plt
# import base64
import pandas as pd

# Popup setup...
resolution = 60  # originally set at 75 (try lower for size??)
width = 4.5
height = 2.5
# (-2.1967272417616583, 22.6318359375)
CENTER_START = [-1.65, 33.3]
ZOOM_START = 5 

# def plotDot(point):
#     '''input: series that contains a numeric named latitude and a numeric named longitude
#     this function creates a CircleMarker and adds it to your this_map'''
#     folium.CircleMarker(location=[point.lat, point.lon],
#                         radius=2,
#                         weight=5).add_to(m)

def build_plot(lat, lng, gdf, filepath):
    pt = Point(lng, lat)
    mask = gdf.geometry.intersects(pt)

    subset = gdf[mask]
    cell_id = subset.CELLID.values[0]

    if len(subset) == 0:
        return
    
    sel = [("CELLID", "=", cell_id)] # , engine="pyarrow"
    gdf_cellid = gpd.read_parquet(filepath, filters=sel)

    dt_rng = pd.to_datetime(gdf_cellid.columns.values[2:32])
    vals_arr = gdf_cellid.iloc[:, 2: 32].values.T
    cols = gdf_cellid.ens.values

    df = pd.DataFrame(vals_arr, columns=cols, index=dt_rng)

    st.write(cell_id)
    st.line_chart(df)

    return df

@st.cache_data
def load_geodataframe(filepath: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_parquet(filepath)
    gdf = gdf.drop_duplicates("geometry")
    return gdf

@st.cache_data
def create_geojson_feature(_gdf: gpd.GeoDataFrame) -> folium.GeoJson:
    geom = gdf.drop_duplicates("geometry").geometry
    geo_j = geom.to_json()
    geo_j = folium.GeoJson(data=geo_j)
    return geo_j

m = folium.Map(location=CENTER_START, zoom_start=5)

# if "center" not in st.session_state:
#     st.session_state["center"] = CENTER_START
# if "zoom" not in st.session_state:
#     st.session_state["zoom"] = ZOOM_START
# if "markers" not in st.session_state:
#     st.session_state["markers"] = []

# if "counter" not in st.session_state:
#     st.session_state.counter = 0

# df = pd.read_parquet("/Users/sam/chandana/shapeApr282024/cell_centers.parquet")
# df.apply(plotDot, axis=1)



# st.session_state["markers"] = [
#     folium.CircleMarker(
#         location=[point.lat, point.lon],
#         radius=2,
#         weight=5
#     ) for point in df.itertuples()
# ]

# for marker in st.session_state["markers"]:
#     fg.add_child(marker)

filename = "S585-R99p"
filepath = f"../data/{filename}.parquet"

gdf = load_geodataframe(filepath)

geo_j = create_geojson_feature(gdf)

fg = folium.FeatureGroup(name="Markers")
fg.add_child(geo_j)

# Render folium map in streamlit and listen to the screen
st_data = st_folium(
    m,
    center=CENTER_START,
    zoom=ZOOM_START,
    feature_group_to_add=fg,
    key="new",    
    height=650,
    width=1200
)



data = None
if st_data.get("last_clicked"):
    data = build_plot(st_data["last_clicked"]["lat"], st_data["last_clicked"]["lng"], gdf, filepath)

if data is not None:
    st.write(data) # Writes to the app
    st.write((st_data["last_clicked"]["lat"], st_data["last_clicked"]["lng"]))

# if __name__ == "__main__":
