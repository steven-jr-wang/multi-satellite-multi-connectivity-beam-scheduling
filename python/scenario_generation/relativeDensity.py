import fiona
import geopandas as gpd


import pandas as pd
from shapely.geometry import Polygon
import folium
import webbrowser

# ---------------------------
# 1. Read cell CSV data and construct polygons (center cells only)
# ---------------------------
csv_file = '../../data/cell_clusters/output_cell_cluster1.csv'
df = pd.read_csv(csv_file)


def create_polygon(row):
    # Vertices in the CSV: vertex1_lat, vertex1_lon, ..., vertex6_lat, vertex6_lon
    vertices = [
        (row['vertex1_lon'], row['vertex1_lat']),
        (row['vertex2_lon'], row['vertex2_lat']),
        (row['vertex3_lon'], row['vertex3_lat']),
        (row['vertex4_lon'], row['vertex4_lat']),
        (row['vertex5_lon'], row['vertex5_lat']),
        (row['vertex6_lon'], row['vertex6_lat'])
    ]
    # Ensure the polygon is closed by repeating the first vertex at the end if necessary
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    return Polygon(vertices)


df['geometry'] = df.apply(create_polygon, axis=1)
gdf_cells = gpd.GeoDataFrame(df, geometry='geometry', crs="EPSG:4326")

# Filter to retain only the center-region cells
gdf_cells = gdf_cells[gdf_cells['region'] == 'center']

# ---------------------------
# 2. Read the Kontur GPKG population data and compute population distribution
# ---------------------------
gpkg_file = '../../data/kontur_population_CN_20220630.gpkg'
# read the layer "population"
gdf_population = gpd.read_file(gpkg_file, layer="population")
print("CRS of population data (GPKG):", gdf_population.crs)
print("Initial CRS of cell data:", gdf_cells.crs)

# Reproject cell data to match GPKG CRS if their coordinates are different
if gdf_population.crs != gdf_cells.crs:
    gdf_cells = gdf_cells.to_crs(gdf_population.crs)

# ---------------------------
# Compute the total population of each center cell by spatial intersection
# ---------------------------
pop_sum_list = []
for idx, cell in gdf_cells.iterrows():
    # Find all Kontur population units that intersect with the current cell
    intersecting = gdf_population[gdf_population.intersects(cell.geometry)]
    pop_sum = 0
    for i, pop_row in intersecting.iterrows():
        inter_geom = cell.geometry.intersection(pop_row.geometry)
        if not inter_geom.is_empty:
            # Compute the overlapping area ratio and proportionally allocate the population
            fraction = inter_geom.area / pop_row.geometry.area
            pop_sum += pop_row['population'] * fraction
    pop_sum_list.append(pop_sum)

gdf_cells['pop_sum'] = pop_sum_list

# Compute total population and normalized relative density for center cells
total_population = gdf_cells['pop_sum'].sum()
gdf_cells['relative_density'] = gdf_cells['pop_sum'] / total_population

# ---------------------------
# 4. Visualize center cells with interactive map using Folium
# ---------------------------
# Reproject back to EPSG:4326 (WGS84 lat/lon) for proper Folium display
if gdf_cells.crs != "EPSG:4326":
    gdf_cells = gdf_cells.to_crs("EPSG:4326")

# Estimate map center by averaging centroids of all center cells
centroids = gdf_cells.geometry.centroid
center_lat = centroids.y.mean()
center_lng = centroids.x.mean()

# Color mapping function: interpolate from light red to deep red based on normalized density
m = folium.Map(location=[center_lat, center_lng], zoom_start=8, tiles='cartodbpositron')


# Color mapping function: interpolate from light red to deep red based on normalized density
def get_color(value, min_val, max_val):
    # Normalize value t to [0,1]
    if max_val > min_val:
        t = (value - min_val) / (max_val - min_val)
    else:
        t = 0
    # Color transitions smoothly from (1, 0.8, 0.8) to (1, 0, 0)
    r = 1
    g = 0.8 * (1 - t)
    b = 0.8 * (1 - t)
    # Convert RGB values to 0–255 integer range
    r_int = int(r * 255)
    g_int = int(g * 255)
    b_int = int(b * 255)
    return '#{:02x}{:02x}{:02x}'.format(r_int, g_int, b_int)


min_density = gdf_cells['relative_density'].min()
max_density = gdf_cells['relative_density'].max()

# Add each center cell to the Folium map
for idx, row in gdf_cells.iterrows():
    geom = row.geometry
    # For Polygon use outer ring; for MultiPolygon take the first polygon
    if geom.geom_type == 'Polygon':
        coords = [[pt[1], pt[0]] for pt in list(geom.exterior.coords)]
    else:
        coords = [[pt[1], pt[0]] for pt in list(geom.geoms[0].exterior.coords)]

    fill_color = get_color(row['relative_density'], min_density, max_density)
    popup_text = f"Relative Density: {row['relative_density']:.4f}"

    folium.Polygon(
        locations=coords,
        color='black',
        weight=1,
        fill=True,
        fill_color=fill_color,
        fill_opacity=0.7,
        popup=popup_text
    ).add_to(m)

# Visualize the rendered map to an HTML file and open it in the browser
map_file = '../../figures/cell_density_map.html'
m.save(map_file)
webbrowser.open(map_file)
# Save each center cell's index and relative density to CSV
output_csv = '../../data/relative_density.csv'
gdf_cells[['index', 'relative_density']].to_csv(output_csv, index=False)
print(f"Relative density data saved to {output_csv}")
