# Import necessary libraries
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
                    merged_lines = MultiLineString(merged_lines)
                    merged_lines = list(merged_lines.geoms)

            # Create polygons from merged lines
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

# Main function to execute the entire workflow
def main():
    
    image_path = 'boe_3_georeferenced.tif'
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
    category_data = {
        'Valle_1': {
            'color': [47, 110, 246],  # Green
            'range': [(37, 57), (100, 120), (235, 255)],
            'type': 'Polygon'
        },
        'Valle_2': {
            'color': [230, 72, 43],  # Blue
            'range': [(220, 240), (72, 82), (43, 53)],
            'type': 'Polygon'
        },
        'Linea_subterranea_220kV': {
            'color': [234, 51, 192],  # Yellow
            'range': [(224, 244), (51, 61), (192, 202)],
            'type': 'LineString'
        }
    }

    # Execute the processing workflow
    process_raster_with_category_data(image_path, category_data)
    create_new_raster_from_category_data(image_path, raster_path, category_data)
    apply_nearest_neighbor_filter(raster_path, output_raster_path, category_data)
    geojsons = convert_raster_to_geojson(output_raster_path, output_geojson_path, category_data)
    unificar_geojson(geojsons, geojson_path)
    #geojson_to_shp(geojson_path, shp_path)

# Execute the main function if this script is run
if __name__ == "__main__":
    main()