# streamlit run 3dprint.py
# Georeference-BOE

## Load File
- Offer the option to load a PDF or an image.

## PDF Processing (if applicable)

- If a PDF is loaded, execute the image extraction script.
- Proceed with the same flow as loading an image.

## Image Processing

- The image is initially not georeferenced.
- Prompt the user to click on a point.
- Record the pixel coordinates (`pixel_x1`, `pixel_y1`).
- Open a map using Leaflet, centered over Spain, with the option to search for locations and go fullscreen.
- The map should be a base map of the PNOA (URL in the `geo_image_map.py` script), and optionally, add another base map like Google Maps.
- The user locates the same point on the map, places a marker, and the coordinates are recorded (`map_y1`, `map_x1`).
- Repeat this process three times, considering the possibility of user errors and providing an option to restart.

## Georeferencing

- Execute the georeferencing process using the obtained coordinates (last part of the `geo_image_map.py` script).
- Provide an option to download the georeferenced image in TIFF format.

## Image Overlay on Map

- Visualize the georeferenced image on a map.
- Use Leaflet, PNOA, Google Maps, and allow toggling visibility and adjusting transparency of the image.

## Data Collection

- Display the image to the user again.
- Allow the user to fill a dictionary at the end of the `vectorize.py` script.
- User chooses the keys (types) from options like VALLADO, VIALES, PANELES, POWERBLOCK, MT, SET, PLATAFORMAS, TALUDES, VIROLAS, OTROS.
- User clicks on the image to assign RGB values to each selected key.

## Vectorization

- After the user selects all desired types, execute the `vectorize.py` script.
- Generate an output folder containing necessary files for a shapefile and a 'FV_TRAMITACION.geojson' file.
- Provide an option to download in either .shp or .geojson format.

## Abort and Error Handling

- Consider the possibility of the user aborting the process at any point.
- Handle possible errors and provide user-friendly messages.

