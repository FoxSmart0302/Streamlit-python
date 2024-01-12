import streamlit as st
import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_agg import RendererAgg
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from ipyleaflet import Map, Marker, TileLayer, AwesomeIcon, SearchControl, FullScreenControl
from streamlit_js_eval import streamlit_js_eval
from PIL import Image
from io import BytesIO
from osgeo import gdal
from streamlit_drawable_canvas import st_canvas
import folium
from folium import plugins
from streamlit_folium import st_folium
from geopy.geocoders import Nominatim
from geopy import geocoders
from geopy.point import Point
from datetime import datetime
current_time = datetime.now()

import gmaps
# *********vectorize start*******
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd

from scipy.ndimage import convolve
from osgeo import gdal, gdalconst
from rasterio.features import shapes
from shapely.geometry import shape, LineString, MultiLineString, Polygon
from shapely.ops import linemerge
from skimage.morphology import skeletonize
# *********vectorize end*********

print("------------------pre-loading-------------------")
# Initialization
if 'position_list' not in st.session_state:
    st.session_state['position_list'] = []

if 'display_position_list' not in st.session_state:
    st.session_state['display_position_list'] = []

if 'map_position_list' not in st.session_state:
    st.session_state['map_position_list'] = []

if 'upload_file_name' not in st.session_state:
    st.session_state['upload_file_name'] = ''

if 'show_map_image' not in st.session_state:
    st.session_state['show_map_image'] = None

if 'select_key_option' not in st.session_state:
    st.session_state['select_key_option'] = 'VALLADO'

if 'select_type_option' not in st.session_state:
    st.session_state['select_type_option'] = 'Polygon'

if 'step' not in st.session_state:
    st.session_state['step'] = 0

if 'category_item' not in st.session_state:
    st.session_state['category_item'] = []

if 'file_current_time' not in st.session_state:
    st.session_state['file_current_time'] = []

st.sidebar.title("Georeferenced Image")
st.markdown(
    """
    <style>
    .stephome {
        display: flex;
        align-items: center;
        margin-bottom: 10px;
    }
    .stepletter {
        color: #813ed7;
        margin-left: 5px;
        font-size: 20px;
        font-weight: bold;
    }
    .circle {
        display: flex;
        justify-content: center;
        align-items: center;
        width: 30px;
        height: 30px;
        border-radius: 50%;
        background-color: #813ed7;
        color: #fff;
        font-weight: bold;
        margin-right: 5px;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# Display numbered circles with letters in the sidebar
numbers = range(1, 6)
letters = ['A', 'B', 'C', 'D', 'E']

def step_show(num, letter):
    circle_html = f"<div class='stephome'><div class='circle'>{num}</div><div class='stepletter'>{letter}</div></div>"
    st.sidebar.markdown(f"{circle_html}", unsafe_allow_html=True)


def get_image_bounds(image_path):
    # Open the image using GDAL
    dataset = gdal.Open(image_path)
    if dataset is None:
        st.error("Could not open the image")
        return None
    # Extract geotransformation information to get image bounds
    geotransform = dataset.GetGeoTransform()
    min_x = geotransform[0]
    max_x = min_x + geotransform[1] * dataset.RasterXSize
    max_y = geotransform[3]
    min_y = max_y + geotransform[5] * dataset.RasterYSize
    # Define bounds dictionary
    bounds = {
        "min_x": min_x,
        "max_x": max_x,
        "min_y": min_y,
        "max_y": max_y
    }
    # Convert bounds to corner coordinates
    esquina_inf_izq = (bounds['min_y'], bounds['min_x'])
    esquina_sup_der = (bounds['max_y'], bounds['max_x'])
    return [esquina_inf_izq, esquina_sup_der]
    
def get_coordinates_from_pixel(pixel_x, pixel_y, image_width, image_height, bounds):
    lower_lat, lower_lon = bounds[0]
    upper_lat, upper_lon = bounds[1]

    normalized_x = pixel_x / image_width
    normalized_y = pixel_y / image_height

    latitude = (1 - normalized_y) * lower_lat + normalized_y * upper_lat
    longitude = (1 - normalized_x) * lower_lon + normalized_x * upper_lon

    return latitude, longitude

def get_position_from_image(col1, uploaded_file):
    position_list = []
    display_position_list = []
    image_height = 800
    image_width = 1200
    image = None
    if uploaded_file is not None and st.session_state['step'] == 0:
        if st.session_state['upload_file_name'] != uploaded_file.name:
            st.session_state['position_list'] = []
            st.session_state['map_position_list'] = []

        st.session_state['upload_file_name'] = uploaded_file.name
        # get bound of georeferenced image end
        with open("temp_map_file", "wb") as file:
            file.write(uploaded_file.getvalue())
        # Convert TIF to PNG and get image bounds
        image_png = tif_to_png("temp_map_file", 'upload/image.png')
        bounds = get_image_bounds("temp_map_file")

        image = Image.open(BytesIO(uploaded_file.read()))
        image_width = image.width
        image_height = image.height
        # Display the canvas for interaction
        with col1:
            canvas_result = st_canvas(
                fill_color="#FFF",  # Color for marking points
                stroke_width=0,  # Stroke width for marking points
                update_streamlit=True,
                drawing_mode="point",
                height=image.height,  # Set the canvas height to match the image height
                width=image.width,  # Set the canvas width to match the image width
                background_image=image,  # Convert the image to a byte string
                key="image_data",  # Key for storing the image data
                
            )
            
            
        # Custom CSS for marker image
        marker_css = """
        <style>
            .marker-container {
                position: relative;
                width: 30px;  /* Adjust the size of the marker image as needed */
                height: 40px;
            }
            .marker-image {
                position: absolute;
                width: 100%;
                height: 100%;
                background-image: url("https://cdn.jsdelivr.net/npm/leaflet@1.9.3/dist/images/marker-icon.png");
                background-repeat: no-repeat;
                background-position: center;
                pointer-events: none;
            }
            .marker-letter {
                position: absolute;
                top: -15px;
                left: 50%;
                transform: translate(-50%, -50%);
                color: black;
                font-weight: bold;
                font-size: 20px;
            }
        </style>
        """
        # Draw marker image instead of dots
        st.markdown(marker_css, unsafe_allow_html=True)
        if image is not None and canvas_result.json_data is not None:
            for data in canvas_result.json_data["objects"]:
                if data["type"] == "circle":
                    x = data["left"]
                    y = data["top"]
                    # latitude, longitude
                    # latitude, longitude = get_coordinates_from_pixel(x, y, image.width, image.height, bounds)
                    # position_list.append((round(x, 1), round(y, 1)))
                    letter = chr(ord("A") + len(position_list))

                    position_list.append((x, y))
                    display_position_list.append((round(x, 1), round(y, 1)))

                    st.session_state['position_list'] = position_list
                    st.session_state['display_position_list'] = display_position_list

                    if len(position_list) > 3:
                        position_list = position_list[:3]
                        display_position_list = display_position_list[:3]

                        st.session_state['position_list'] = position_list
                        st.session_state['display_position_list'] = display_position_list
                        break
                    
                     # Create the marker container with image and letter
                    marker_html = f"""
                    <div class="marker-container" style="left: {x - 15}px; top: {y - 40}px;">
                        <div class="marker-image"></div>
                        <div class="marker-letter">{letter}</div>
                    </div>
                    """

                    st.markdown(marker_html, unsafe_allow_html=True)
                    
        # Truncate the position_list to length 3 if it exceeds
        # position_list = position_list[:3]
    if uploaded_file is not None:
        num_points = st.sidebar.text_input("Three points of the image", value=st.session_state['display_position_list'])
    return image_width, image_height

def get_color_from_image(uploaded_file):
    image_height = 800
    image_width = 1200
    image = None
    if uploaded_file is not None:
        # get bound of georeferenced image end
        with open("temp_map_file", "wb") as file:
            file.write(uploaded_file.getvalue())
        # Convert TIF to PNG and get image bounds
        image_png = tif_to_png("temp_map_file", 'upload/image.png')
        bounds = get_image_bounds("temp_map_file")

        image = Image.open(BytesIO(uploaded_file.read()))
        image_width = image.width
        image_height = image.height
        # Display the canvas for interaction
        
        canvas_result = st_canvas(
            fill_color="#FFF",  # Color for marking points
            stroke_width=1,  # Stroke width for marking points
            update_streamlit=True,
            drawing_mode="point",
            height=image.height,  # Set the canvas height to match the image height
            width=image.width,  # Set the canvas width to match the image width
            background_image=image,  # Convert the image to a byte string
            key="image_data",  # Key for storing the image data
            
        )
            
            
        # Custom CSS for marker image
        marker_css = """
        <style>
            .marker-image {
                position: absolute;
                width: 30px;  /* Adjust the size of the marker image as needed */
                height: 40px;
                background-image: url("https://cdn.jsdelivr.net/npm/leaflet@1.9.3/dist/images/marker-icon.png");
                background-repeat: no-repeat;
                background-position: center;
                pointer-events: none;
            }
        </style>
        """
        # Draw marker image instead of dots

        st.markdown(marker_css, unsafe_allow_html=True)
        rgb = None
        if image is not None and canvas_result.json_data is not None:
            for data in canvas_result.json_data["objects"]:
                if data["type"] == "circle":
                    x = data["left"]
                    y = data["top"]
        
                    pixel = image.getpixel((x, y))
                    rgb = f"RGB: {pixel[:3]}"
                    # latitude, longitude
                    # latitude, longitude = get_coordinates_from_pixel(x, y, image.width, image.height, bounds)
                    # position_list.append((round(x, 1), round(y, 1)))
                    
                    st.markdown(
                        f'<div class="marker-image" style="left: {x - 11}px; top: {y - 39}px;"></div>',
                        unsafe_allow_html=True,
                    )

            after_dots_length = len(canvas_result.json_data["objects"])
            print("====after_dots_length:", after_dots_length)       
            temp_category = st.session_state['category_item']
            if st.session_state['step'] == 2 and rgb is not None and len(st.session_state['category_item']) != after_dots_length:
                json_key = st.session_state['select_key_option']
                json_type = st.session_state['select_type_option']
                print("====json_key:", json_key)
                print("====json_type:", json_type)
                print("====rgb:", rgb)
                json_object = {
                    'key' : json_key,
                    'type' : json_type,
                    'rgb' : rgb
                }
                temp_category.append(json_object)
                st.session_state['category_item'] = temp_category
       

        # Truncate the position_list to length 3 if it exceeds
        # position_list = position_list[:3]
  
def tif_to_png(input_path, output_path):
    # Convert TIF image to PNG format using PIL
    image_tif = Image.open(input_path)
    image_tif.save(output_path)
    return output_path

def create_folium_map(image_tif):
    # Create a folium map with a raster layer overlay
    m = folium.Map(location=[40, -4], zoom_start=8,
        tiles="https://tms-pnoa-ma.idee.es/1.0.0/pnoa-ma/{z}/{x}/{-y}.jpeg",
        attr="PNOA",)

    # Convert TIF to PNG and get image bounds
    image_png = tif_to_png(image_tif, 'image.png')
    print("=====image_png:", image_png)
    bounds = get_image_bounds(image_tif)
    print("=====image_bound", bounds)
    # Check if PNG file exists and add it to the folium map
    if not os.path.isfile(image_png):
        st.error(f"Could not find {image_png}")
    else:
        opacity = st.sidebar.slider("Opacity", 0.0, 1.0, 0.6, 0.1)  # Add opacity slider
        img = folium.raster_layers.ImageOverlay(
            name=image_png,
            image=image_png,
            bounds=bounds,
            opacity=opacity,
            interactive=True,
            cross_origin=False,
            zindex=1,
        )
        img.add_to(m)

        # Add 3 markers to the map
        # marker_locations = [[40.1, -4.2], [40.2, -4.3], [40.3, -4.4]]
        # marker_list = []
        # for marker_location in marker_locations:
        #     marker = folium.Marker(location=marker_location, popup=folium.Popup(str(marker_location)))
        #     marker.add_to(m)
        #     marker_list.append(marker)
        
        folium.LayerControl().add_to(m)

    return m

def get_pos(lat,lng):
    return lat,lng

def folium_map_coordinates(image_width, image_height):
    # Create a Folium map
    m = folium.Map(location=[40, -4], zoom_start=8)
     # Add search control to the map
    # plugins.Geocoder().add_to(m)
    
    # Add click event listener to the map
    m.add_child(folium.ClickForMarker(popup=None))  # Add a marker on click
    m.add_child(folium.LatLngPopup())  # Display lat/lon on click
    # Display the map in Streamlit
    map = st_folium(m, width=image_width, height=image_height)
    data = None
    if map.get("last_clicked"):
        data = get_pos(map["last_clicked"]["lat"], map["last_clicked"]["lng"])

    if data is not None:
        st.write(data) # Writes to the app
        tempData = st.session_state['map_position_list']
        if len(tempData) == 0:
            print("Empty_insert:", data)
            tempData.append(data)
        else:
            last_element = tempData[-1]
            # print("last_element:", last_element)
            # print("current_element:", data)
            if last_element != data:
                tempData.append(data)
                if len(tempData) > 3:
                    tempData = tempData[:3]
                st.session_state['map_position_list'] = tempData

def main_style():
    st.markdown(
        f"""
        <style>
        .main .block-container {{
            max-width: 100%;
        }}
        .main .element-container{{
            position: absolute;
        }}
        .css-1aumxhk {{
            width: 100% !important;
        }}
        .stButton button,
        .stDownloadButton button
         {{
            width: 100%;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )

def generate_unique_filename(filename):
    base_name, extension = os.path.splitext(filename)
    counter = 1
    while os.path.exists(filename):
        filename = f"{base_name}_{counter}{extension}"
        counter += 1
    return filename

def georeferenced_image_download(url):
    # File path of the georeferenced image
    georeferenced_image_path = url
    file_name = os.path.basename(url)
    print("image_upload_url:", url, file_name)

    # Check if the georeferenced image exists
    print("georeferenced_image_path:", georeferenced_image_path)

    if os.path.exists(georeferenced_image_path) and os.path.isfile(georeferenced_image_path):
        # Display a download button
        with open(georeferenced_image_path, "rb") as file:
            st.sidebar.download_button(
                label="Download Georeferenced Image",
                data=file,
                file_name=file_name,
            )
    else:
        st.markdown("Georeferenced image not found.")

def geojson_file_download(url):
    # File path of the georeferenced image
    georeferenced_image_path = url
    file_name = os.path.basename(url)
    print("image_upload_url:", url, file_name)

    # Check if the georeferenced image exists
    if os.path.exists(georeferenced_image_path) and os.path.isfile(georeferenced_image_path):
        # Display a download button
        with open(georeferenced_image_path, "rb") as file:
            st.sidebar.download_button(
                label="Download GEOJSON File",
                data=file,
                file_name=file_name,
            )
    else:
        st.markdown("Georeferenced image not found.")

def show_image_map(url):
    image_tif = url
    folium_map = create_folium_map(image_tif)  
    folium_map.save('map.html')
    map_html = open('map.html', 'r').read()
    st.components.v1.html(map_html, height=1100)

def image_map_with_opaticy(session_image_path):
    print("----------------------------session_image_show---------------------------------")
    if session_image_path is not None:
        show_image_map(session_image_path)
        if st.sidebar.button("Next"):
            st.session_state['step'] = 2
            st.session_state['category_item'] = []
            st.experimental_rerun()
    
def init_button(uploaded_file):
    if uploaded_file is not None:
        # Button 1: Init Map
        if st.sidebar.button("Init"):
            # st.experimental_rerun()
            streamlit_js_eval(js_expressions="parent.window.location.reload()")

def select_option_with_type():
    # Define options for the select box
    # rgb_uploaded_file = st.sidebar.file_uploader("Upload an image for dictionary", type=['png', 'tif'])
    # if rgb_uploaded_file is not None:
    #     st.image(rgb_uploaded_file, caption='Uploaded Image', width=None)
    type_options = ['Polygon', 'LineString']
    key_options = ['VALLADO', 'VIALES', 'PANELES', 'POWERBLOCK', 'MT', 'SET', 'PLATAFORMAS', 'TALUDES', 'VIROLAS', 'Otros']

    # Create a select box
    selected_key_option = st.sidebar.selectbox('Select an key', key_options)
    selected_type_option = st.sidebar.selectbox('Select an type', type_options)
    category_item_count = st.sidebar.text_input("Category items count", len(st.session_state['category_item']))
    st.session_state['select_key_option'] = selected_key_option
    st.session_state['select_type_option'] = selected_type_option
        # print("select_key_option", st.session_state['select_key_option'])
        # print("select_type_option", st.session_state['select_type_option'])

# **********vectorize start***********

# Function to process raster data based on specified category data
def process_raster_with_category_data(raster_path, category_data):
    with rasterio.open(raster_path, 'r+') as dataset:
        bandas_rgb = dataset.read()
        alto, ancho = bandas_rgb.shape[1:]

        # Iterate through each pixel in the raster
        for x in range(ancho):
            for y in range(alto):
                color_actual = list(bandas_rgb[:, y, x])
                coincidencia = False

                # Check for each category's color range
                for categoria, data in category_data.items():
                    if all(color_range[0] <= color_val <= color_range[1] for color_val, color_range in zip(color_actual, data['range'])):
                        bandas_rgb[:3, y, x] = data['color']
                        coincidencia = True
                        break

                # If no match is found, set the color to [0, 0, 0]
                if not coincidencia:
                    bandas_rgb[:3, y, x][:3] = [0, 0, 0]

        # Write the modified data back to the raster
        dataset.write(bandas_rgb)

    print("Successfully overwritten bands in the TIFF raster file.")

# Function to create a new raster file based on specified category data using a reference raster
def create_new_raster_from_category_data(reference_raster_path, output_raster_path, category_data):
    reference_raster = gdal.Open(reference_raster_path)
    reference_band = reference_raster.GetRasterBand(1)

    # Extract information from the reference raster
    geotransform = reference_raster.GetGeoTransform()
    projection = reference_raster.GetProjection()
    width = reference_band.XSize
    height = reference_band.YSize
    banda_roja = reference_raster.GetRasterBand(1).ReadAsArray()
    banda_verde = reference_raster.GetRasterBand(2).ReadAsArray()
    banda_azul = reference_raster.GetRasterBand(3).ReadAsArray()

    # Iterate through each category to create a new raster file
    for idx, (categoria, data) in enumerate(category_data.items(), start=1):
        output_path = output_raster_path + f'_{idx}.tif'
        driver = gdal.GetDriverByName('GTiff')
        output_raster = driver.Create(output_path, width, height, 1, gdalconst.GDT_Float32)

        output_raster.SetGeoTransform(geotransform)
        output_raster.SetProjection(projection)

        output_band = output_raster.GetRasterBand(1)

        data_array = np.zeros((height, width), dtype=np.uint8)

        range_list = data['range']
        mask = (
            (banda_roja >= range_list[0][0]) & (banda_roja <= range_list[0][1]) &
            (banda_verde >= range_list[1][0]) & (banda_verde <= range_list[1][1]) &
            (banda_azul >= range_list[2][0]) & (banda_azul <= range_list[2][1])
        )
        data_array[mask] = idx

        output_band.WriteArray(data_array)

    reference_raster = None
    output_raster = None

    print(f"New raster created successfully at: {output_raster_path}")

# Function to apply a nearest neighbor filter to smooth a raster file
def apply_nearest_neighbor_filter(raster_path, output_path, category_data):
    for idx, (categoria, data) in enumerate(category_data.items(), start=1):
        path = raster_path + f'_{idx}.tif'
        output_raster_path = output_path + f'_{idx}.tif'
            
        # Open the raster file
        ds = gdal.Open(path)
        banda = ds.GetRasterBand(1)
        array_raster = banda.ReadAsArray()

        # Define a 3x3 kernel for the convolution
        kernel = np.array([[1, 1, 1],
                           [1, 1, 1],
                           [1, 1, 1]])

        # Apply convolution and set threshold values
        resultado = convolve(array_raster, kernel, mode='constant', cval=0.0)
        resultado[resultado <= 0] = 0
        resultado[resultado > 0] = idx
        
        if data['type'] == 'LineString':
            # Apply skeletonization
            skeleton_result = skeletonize(resultado.astype(bool))
             # Write the result to a new raster file
            driver = gdal.GetDriverByName('GTiff')
            ds_resultado = driver.Create(output_raster_path, ds.RasterXSize, ds.RasterYSize, 1, banda.DataType)
            ds_resultado.GetRasterBand(1).WriteArray(skeleton_result)

            ds_resultado.SetGeoTransform(ds.GetGeoTransform())
            ds_resultado.SetProjection(ds.GetProjection())

            ds = None
            ds_resultado = None
            
        else:
            # Write the result to a new raster file
            driver = gdal.GetDriverByName('GTiff')
            ds_resultado = driver.Create(output_raster_path, ds.RasterXSize, ds.RasterYSize, 1, banda.DataType)
            ds_resultado.GetRasterBand(1).WriteArray(resultado)

            ds_resultado.SetGeoTransform(ds.GetGeoTransform())
            ds_resultado.SetProjection(ds.GetProjection())

            ds = None
            ds_resultado = None
        
def convert_raster_to_geojson(raster_path, output_geojson_path, category_data):
    geojsons = []
    for idx, (categoria, data) in enumerate(category_data.items(), start=1):
        path = raster_path + f'_{idx}.tif'
        output_path = output_geojson_path + f'_{idx}.geojson'

        # Open the raster file and read raster array
        with rasterio.open(path) as src:
            raster_array = src.read(1)

        poligonos = []
        for poligono, valor in shapes(raster_array, transform=src.transform):
            if valor in [1, 2, 3]:
                poligonos.append(shape(poligono))

        if data['type'] == 'Polygon':
            lineas = [LineString(poligono.exterior.coords) for poligono in poligonos]
            merged_lines = linemerge(lineas)

            # Check if there are lines to merge
            if merged_lines.is_empty:
                print("No lines to merge.")
            else:
                if isinstance(merged_lines, MultiLineString):
                    merged_lines = list(merged_lines.geoms)
                    merged_lines = MultiLineString(merged_lines)

            # Create polygons from merged lines
            merged_lines = list(merged_lines.geoms)
            geometries = [Polygon(line) for line in merged_lines]

        elif data['type'] == 'LineString':
            lineas = [LineString(poligono.exterior.coords) for poligono in poligonos]
            merged_lines = linemerge(lineas)

            # Check if there are lines to merge
            if merged_lines.is_empty:
                print("No lines to merge.")
            else:
                if isinstance(merged_lines, MultiLineString):
                    # Extract coordinates from each line in MultiLineString
                    merged_lines = list(merged_lines.geoms)
                    geometries = [LineString(line.coords) for line in merged_lines]
                else:
                    # Convert single LineString to a list of LineStrings
                    geometries = [LineString(merged_lines.coords)]

        else:
            print(f"Unsupported geometry type: {data['type']}")
            continue

        gdf_geometries = gpd.GeoDataFrame(geometry=geometries, crs=src.crs)
        gdf_geometries['OCUPACION'] = categoria
        gdf_geometries.to_crs(epsg=25830)
        gdf_geometries.to_file(output_path, driver="GeoJSON")
        geojsons.append(output_path)

    return geojsons

# Function to unify multiple GeoJSON files into a single GeoJSON file
def unificar_geojson(file_paths, output_path):
    gdframes = [gpd.read_file(file_path) for file_path in file_paths]

    gdf_unificado = gpd.GeoDataFrame(pd.concat(gdframes, ignore_index=True))
    gdf_unificado.to_crs(epsg=25830)
    gdf_unificado.to_file(output_path, driver='GeoJSON')

    print(f"Unified and saved at {output_path}")

# Function to convert GeoJSON to ESRI Shapefile
def geojson_to_shp(geojson_path, shp_path):
    gdf = gpd.read_file(geojson_path)
    gdf.to_file(shp_path, driver='ESRI Shapefile')
    print(f"Converted and saved at {shp_path}")

def vectorize_main(uploaded_file):
    if st.sidebar.button("Vectorizar"):
        data = st.session_state['category_item']
        # print("pre_category_data:", data)
        category_data = {}

        for item in data:
            key = item['key']
            color = [int(x) for x in item['rgb'][len('RGB: '):][1:-1].split(', ')]
            color_range = []
            for value in color:
                lower_bound = max(value - 10, 0)
                upper_bound = min(value + 10, 255)
                color_range.append((lower_bound, upper_bound))
            
            category_data[key] = {
                'color': color,
                'range': color_range,  # Add the desired range values
                'type': item['type']
            }

        
        # with open("uploaded_image", "wb") as file:
        #     file.write(uploaded_file.getvalue())
        # image_path = "uploaded_image"
        image_path = uploaded_file
        raster_path = 'nuevo/nuevo_raster'
        output_raster_path = 'raster/raster_suavizado'
        output_geojson_path = 'geojson/FV_TRAMITACIÓN'
        shp_path = 'shp/FV_TRAMITACIÓN.shp'
        geojson_path = 'geojson/FV_TRAMITACIÓN.geojson'
        
        # Category data specifying colors and ranges
        # The value of the key in this dictionary will be the value entered by the user.
        # The value of 'color' will be the RGB value of the pixel on which the user clicks, and 
        # the 'range' will be +- 10 for each of the RGB channels. In conclusion, 
        # the dictionary will be dynamically filled during each execution.
        # category_data = {
        #     'Valle_1': {
        #         'color': [47, 110, 246],  # Green
        #         'range': [(37, 57), (100, 120), (235, 255)],
        #         'type': 'Polygon'
        #     },
        #     'Valle_2': {
        #         'color': [230, 72, 43],  # Blue
        #         'range': [(220, 240), (72, 82), (43, 53)],
        #         'type': 'Polygon'
        #     },
        #     'Linea_subterranea_220kV': {
        #         'color': [234, 51, 192],  # Yellow
        #         'range': [(224, 244), (51, 61), (192, 202)],
        #         'type': 'LineString'
        #     }
        # }
        # print("new_category_data", category_data)
        # Execute the processing workflow
        process_raster_with_category_data(image_path, category_data)
        create_new_raster_from_category_data(image_path, raster_path, category_data)
        apply_nearest_neighbor_filter(raster_path, output_raster_path, category_data)
        geojsons = convert_raster_to_geojson(output_raster_path, output_geojson_path, category_data)
        unificar_geojson(geojsons, geojson_path)
        geojson_file_download('geojson/FV_TRAMITACIÓN.geojson')
        # geojson_to_shp(geojson_path, shp_path)
        print("--------------vectorize----------------")

# **********vectorize end*************

def delete_all_from_directory(directory_path):
    if directory_path:
    # Ensure that the directory path exists
        if os.path.exists(directory_path) and os.path.isdir(directory_path):
            # Get a list of all files within the directory
            file_list = os.listdir(directory_path)
            
            if file_list:
                # Loop through the files and delete them
                for file_name in file_list:
                    file_path = os.path.join(directory_path, file_name)
                    os.remove(file_path)
                
                print("All files deleted successfully.")
            else:
                print("The directory is empty.")
        else:
            print("Invalid directory path.")

def delete_file_from_directory(directory_path, file_to_exclude):
    # Check if the directory path and file to exclude are provided
    if directory_path and file_to_exclude:
        # Ensure that the directory path exists
        if os.path.exists(directory_path) and os.path.isdir(directory_path):
            # Get a list of all files within the directory
            file_list = os.listdir(directory_path)
            
            if file_list:
                # Loop through the files and delete them, except for the specified file
                for file_name in file_list:
                    if file_name != file_to_exclude:
                        file_path = os.path.join(directory_path, file_name)
                        os.remove(file_path)
                
                print("All files except the specified file have been deleted.")
            else:
                print("The directory is empty.")
        else:
            print("Invalid directory path.")

def main():
    
    col1, col2 = st.columns(2)
    main_style()
    map_position_list = []
    
    # show_image_map(f"new/test.tif")
    # if st.sidebar.button("Continue"):
    #     folium_map_coordinates(1200, 800, map_position_list)
    step_show(1, "Location list settings")
    uploaded_file = st.sidebar.file_uploader("Upload an Image", type=['png', 'tif'])
    
    
    image_width, image_height = get_position_from_image(col1, uploaded_file)
    if uploaded_file is not None:
            with col2:
                if st.session_state['step'] == 0:
                    folium_map_coordinates(image_width, image_height - 15)

        
            map_position_list = st.session_state['map_position_list']
                    
            display_position_list = []
            for x, y in map_position_list:
                display_position_list.append((round(x, 1), round(y, 1)))
                
            map_points = st.sidebar.text_input("Three points of the map", display_position_list)
            if len(st.session_state['map_position_list']) == 3:
                # Create two columns within the sidebar layout
                button_col1, button_col2 = st.sidebar.columns(2)

            if (len(st.session_state['map_position_list']) == 3 and len(st.session_state['position_list']) == 3) or st.session_state['step'] == 1 or st.session_state['step'] == 2:
                step_show(2, "Create Image")
                if st.sidebar.button("Georeference"):
                    with open("uploaded_image", "wb") as file:
                        file.write(uploaded_file.getvalue())
                    print("-------------------Georeference--------------------", st.session_state['position_list'])
                    position_list = st.session_state['position_list']
                    current_time_str = current_time.strftime("%Y-%m-%d_%H-%M-%S")
                    print("Formatted_current_time:", current_time_str)
                    st.session_state['file_current_time'] = current_time_str
                    print("position_list:", position_list[1])
                    print("map_position_list:", st.session_state['map_position_list'])
                    pixel_x1, pixel_y1 = position_list[0]
                    map_y1, map_x1 = st.session_state['map_position_list'][0]
                    pixel_x2, pixel_y2 = position_list[1]
                    map_y2, map_x2 = st.session_state['map_position_list'][1]
                    pixel_x3, pixel_y3 = position_list[2]
                    map_y3, map_x3 = st.session_state['map_position_list'][2]
                    print(pixel_x1, pixel_y1, map_y3, map_x3)
                    
                    # Perform georeferencing using gdal_translate and gdalwarp
                    os.system(f"gdal_translate -of GTiff -a_srs EPSG:4326 -gcp {pixel_x1} {pixel_y1} {map_x1} {map_y1} -gcp {pixel_x2} {pixel_y2} {map_x2} {map_y2} -gcp {pixel_x3} {pixel_y3} {map_x3} {map_y3} uploaded_image new/output_georeferenced.tif")
                    os.system(f"gdalwarp -of GTiff -t_srs EPSG:4326 new/output_georeferenced.tif new/georeferenced_{current_time_str}.tif")
                    st.sidebar.text("Successfully created!")
                    
                    #file_download_function
                    georeferenced_image_download(f"new/georeferenced_{current_time_str}.tif")
                    st.session_state['show_map_image'] = f"new/georeferenced_{current_time_str}.tif" 
                    st.session_state['step'] = 1
    # print("show_map_image_session_path", session_image_path)
    
    if st.session_state['step'] == 1:
        image_map_with_opaticy(st.session_state['show_map_image'])

    if st.session_state['step'] == 2:
        get_color_from_image(uploaded_file)
        step_show(3, "Create Dictionary")
        select_option_with_type()
        print("category_item", st.session_state['category_item'])
        if st.session_state['category_item'] is not None:
            vectorize_main(st.session_state['show_map_image'])
    delete_all_from_directory('nuevo')
    delete_all_from_directory('raster')
    delete_all_from_directory('shp')
    delete_file_from_directory('geojson', 'FV_TRAMITACIÓN.geojson')
    init_button(uploaded_file)

if __name__ == "__main__":
    main()