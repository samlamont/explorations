import streamlit as st
import streamlit.components.v1 as components
import httpx
from folium import Map, TileLayer

st.set_page_config(
    page_title="flood mapping demo",
    page_icon=":world_map:️",
    # layout="wide",
)

st.title("Flood Maps")


# The following taken from: https://developmentseed.org/titiler/examples/notebooks/Working_with_CloudOptimizedGeoTIFF_simple/
titiler_endpoint = "https://titiler.xyz"  # Developmentseed Demo endpoint. Please be kind.
url = "https://opendata.digitalglobe.com/events/mauritius-oil-spill/post-event/2020-08-12/105001001F1B5B00/105001001F1B5B00.tif"

# Fetch File Metadata to get min/max rescaling values (because the file is stored as float32)
r = httpx.get(
    f"{titiler_endpoint}/cog/info",
    params={"url": url}).json()

bounds = r["bounds"]

# Fetch File Metadata to get min/max rescaling values (because the file is stored as float32)
r = httpx.get(
    f"{titiler_endpoint}/cog/statistics",
    params={"url": url}).json()

r = httpx.get(
    f"{titiler_endpoint}/cog/tilejson.json",
    params={"url": url}).json()

m = Map(
    location=((bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2),
    zoom_start=13
)

aod_layer = TileLayer(
    tiles=r["tiles"][0],
    opacity=1,
    attr="DigitalGlobe OpenData"
)
aod_layer.add_to(m)

# Add map to streamlit
height = 500
width = 700
components.html(m._repr_html_(), height=height + 10, width=width)

"""
Does this add text? Why yes, yes it does. And you can use `regular markdown` type *stuff*
"""