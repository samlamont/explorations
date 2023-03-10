import streamlit as st
import streamlit.components.v1 as components
import httpx
from folium import Map, TileLayer

st.set_page_config(
    page_title="flood mapping demo",
    page_icon=":world_map:️",
    # layout="wide",
)

st.title("SAR - Flood Mapping Demo")

"""
`Notes` We can potentially use this as a demo page to display flood mapping outputs on an interactive webpage. 
The data (cloud-optimized geotiffs, COGs) can be kept in s3 buckets for cheap and displayed here with a fairly simple
app that is easy to customize.

`What other functionality do you think is needed?`
:thinking_face:

The data shown in the map below is just a sample file from Digital Globe hosted at:
 https://opendata.digitalglobe.com/events/mauritius-oil-spill/post-event/2020-08-12/105001001F1B5B00/105001001F1B5B00.tif
"""

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
    zoom_start=13,
    tiles="cartodbdark_matter"  # "cartodbpositron", "stamentoner"
)

aod_layer = TileLayer(
    tiles=r["tiles"][0],
    opacity=1,
    attr="DigitalGlobe OpenData"
)
aod_layer.add_to(m)

# Add map to streamlit
height = 700
width = 900
components.html(m._repr_html_(), height=height + 10, width=width)

"""
See other StreamLit examples here: https://streamlit.io/gallery
"""

st.info('This is a purely informational message')

with st.form("my_form"):
   st.write("Inside the form")
   slider_val = st.slider("Form slider")
   checkbox_val = st.checkbox("Form checkbox")

   # Every form must have a submit button.
   submitted = st.form_submit_button("Submit")
   if submitted:
       st.write("slider", slider_val, "checkbox", checkbox_val)

st.write("Outside the form")
