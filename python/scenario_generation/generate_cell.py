import h3
import folium
import csv
import webbrowser

cluster_config_mode = 1  # Options: 1, 2, or 3 corresponding to different regions, here we fix at region 1

# Set the H3 resolution
resolution = 4

# Configure the center cell index, wrap-around indices, and output file path based on cluster_config_mode
if cluster_config_mode == 1:
    center_lat = 46
    center_lng = 127
    # centercell_index = '8414963ffffffff'
    # wrap_around = ['84312b3ffffffff', '84312c1ffffffff', '8414b3bffffffff','84149cdffffffff', '8431647ffffffff', '843171dffffffff']
    centercell_index = '8414959ffffffff'
    wrap_around = ['8414a23ffffffff', '8414815ffffffff', '84149ebffffffff','8414921ffffffff', '84312d5ffffffff', '8414b47ffffffff']

    file_path = '../../data/cell_clusters/output_cell_cluster1.csv'

elif cluster_config_mode == 2:
    center_lat = -38
    center_lng = 145

    centercell_index = '84bf48bffffffff'
    wrap_around = ['84bf593ffffffff','84bf431ffffffff','84bf6a9ffffffff','84d8b45ffffffff','84d8b37ffffffff','84d892dffffffff']
    file_path = '../../data/cell_clusters/output_cell_cluster2.csv'

elif cluster_config_mode == 3:

    center_lat = -33
    center_lng = -72

    centercell_index = '84b2c43ffffffff'

    wrap_around = ['84b2e1bffffffff', '84b2ebdffffffff', '84b2ce5ffffffff', '84b2d19ffffffff', '84b2f07ffffffff',
                   '84b2e61ffffffff']
    file_path = '../../data/cell_clusters/output_cell_cluster3.csv'
# Create a Folium map object
m = folium.Map(location=[center_lat, center_lng], zoom_start=8, tiles='cartodbpositron')
polygon_group = folium.FeatureGroup(name="Polygons")

# Initialize the list to store hexagon data
hexagons = []

# Function: Retrieve hexagon center and vertex coordinates
def get_hexagon_data(h3_index, region):
    center = h3.cell_to_latlng(h3_index)
    boundary = h3.cell_to_boundary(h3_index)
    return {
        'h3_index': h3_index,
        'center': center,
        'vertices': boundary,
        'region': region,
        'index': None  # Initially set index to None; will be assigned after sorting
    }

# Generate hexagons in the central region and assign indices based on grid_disk ordering
center_hexagons = [get_hexagon_data(h3_index, 'center') for h3_index in h3.grid_disk(centercell_index, 3)]
for index, hexagon in enumerate(center_hexagons,start=1):
    hexagon['index'] = index

# Function: Calculate distance between two coordinates
def distance(coord1, coord2):
    return h3.great_circle_distance(coord1, coord2, unit='km')


# Calculate latitude and longitude deltas relative to the central cell (index=0)
center_deltas = [
    (
        hex['center'][0] - center_hexagons[0]['center'][0],
        hex['center'][1] - center_hexagons[0]['center'][1]
    )
    for hex in center_hexagons
]

# Match indices for cells in each wrap-around region
for wrap_center in wrap_around:
    wrap_hexagons = [get_hexagon_data(h3_index, f'wrap_{wrap_around.index(wrap_center) + 1}') for h3_index in h3.grid_disk(wrap_center, 3)]
    wrap_center_coord = h3.cell_to_latlng(wrap_center)

    # Match cells
    for delta, center_hex in zip(center_deltas, center_hexagons):
        target_lat = wrap_center_coord[0] + delta[0]
        target_lng = wrap_center_coord[1] + delta[1]

        # Find the nearest hexagon in the wrap-around region to the target coordinates
        best_match = min(wrap_hexagons, key=lambda wrap_hex: distance((target_lat, target_lng), wrap_hex['center']))
        best_match['index'] = center_hex['index']
        hexagons.append(best_match)
        wrap_hexagons.remove(best_match)

# Add central region hexagons
hexagons.extend(center_hexagons)

# Add polygons to the map
for hexagon in hexagons:
    color = 'blue' if hexagon['region'] == 'center' else 'red'
    popup_info = f'Region: {hexagon["region"]}, Index: {hexagon["index"]}, H3: {hexagon["h3_index"]}'
    polygon = folium.Polygon(
        locations=hexagon['vertices'],
        color=color,
        weight=1,
        fill=True,
        fill_opacity=0.1 if hexagon['region'] == 'center' else 0.3,
        popup=popup_info
    )
    polygon_group.add_child(polygon)

# Add FeatureGroup to map
polygon_group.add_to(m)

# Add layer control
folium.LayerControl().add_to(m)

# Save the map to an HTML file
m.save('map.html')
print("Map has been saved to map.html")
webbrowser.open('map.html')


# Save cell data to a CSV file
with open(file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    # Write header
    writer.writerow(['region', 'index', 'h3_index', 'center_lat', 'center_lon',
                     'vertex1_lat', 'vertex1_lon', 'vertex2_lat', 'vertex2_lon',
                     'vertex3_lat', 'vertex3_lon', 'vertex4_lat', 'vertex4_lon',
                     'vertex5_lat', 'vertex5_lon', 'vertex6_lat', 'vertex6_lon'])

    for hexagon in hexagons:
        row = [hexagon['region'], hexagon['index'], hexagon['h3_index'], hexagon['center'][0], hexagon['center'][1]]
        for vertex in hexagon['vertices']:
            row.extend([vertex[0], vertex[1]])
        writer.writerow(row)

print(f"Data saved to {file_path}")