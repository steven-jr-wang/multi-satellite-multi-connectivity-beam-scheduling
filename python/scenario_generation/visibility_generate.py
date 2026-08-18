import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
import astropy.units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy.coordinates import SkyCoord, ITRS, CartesianRepresentation
import cartopy.crs as ccrs
import os
import csv
import h3
import folium
import csv
import webbrowser
import pandas as pd

# ----------------------
# Global parameter settings
# ----------------------

constellation_config_mode = 2 # Choose the satellite constellation configuration mode 840 (DTC Phase-1) or 2016 (DTC Phase-2), 1 or 2
cluster_config_mode = 1  # Select cell cluster 1 within Heilongjiang Province
# ------------------------------------------------------------------
# Simulation parameters
# ------------------------------------------------------------------
P = 28  # number of orbital planes
if constellation_config_mode == 1:

    S = 30  # satellites per plane
    output_dir = r"../../data/snapshots_Phase1/"
else:

    S = 72 # satellites per plane
    output_dir = r"../../data/snapshots_Phase2/"
T = P * S  # total number of satellites

os.makedirs(output_dir, exist_ok=True)

# Inclination and altitude for the two orbital shells
num_planes_incl1 = 14
num_planes_incl2 = 14
incl1 = 53.0 * u.deg
incl2 = 43.0 * u.deg
h1 = 530 * u.km
h2 = 525 * u.km

# Earth radius approximation
R_earth = 6371 * u.km

# Snapshots sampling
start_time = Time("2024-01-01 00:00:00", scale='utc')
end_time = Time("2024-01-01 00:01:00", scale='utc')
num_steps = 10
time_steps = start_time + (end_time - start_time) * np.linspace(0, 1, num_steps)

# Walker Delta constellation configuration
# Compute the RAAN spacing between adjacent orbital planes for each sub-constellation（or formally different orbital shells for Starlink constellation） separately.
delta_RAAN_deg_sub1 = 360.0 / num_planes_incl1
delta_RAAN_deg_sub2 = 360.0 / num_planes_incl2

# Generate orbital parameters for each satellite
# Semi-major axis: a = R_earth + altitude
# Use h1 for the first 14 planes, and h2 for the latter 14 planes
def get_orbit_parameters(p, s, epoch, inc, altitude, delta_RAAN_deg, phase_offset, num_planes, raan_offset=0.0):
    # Apply offset on the orbital plane (RAAN): p + raan_offset (e.g., sub-constellation 2 passes raan_offset=0.5)
    RAAN = (p + raan_offset) * delta_RAAN_deg * u.deg
    a = R_earth + altitude
    sat_id = s + 1
    phase_deg = s * (360.0 / S) + (360.0 / (num_planes * S)) * p + phase_offset
    phase_deg = (phase_deg + 180) % 360 - 180
    nu = phase_deg * u.deg
    # Assume e=0, argp=0
    ecc = 0 * u.one
    argp = 0 * u.deg
    orbit = Orbit.from_classical(
        Earth,
        a,
        ecc,
        inc,
        RAAN,
        argp,
        nu,
        epoch=epoch
    )
    return orbit


def eci_to_geodetic(pos_vec, obstime):
    # Convert ECI coordinates to latitude, longitude, and altitude
    cart_rep = CartesianRepresentation(pos_vec)
    itrs = SkyCoord(cart_rep, frame=ITRS(obstime=obstime))
    lat = itrs.earth_location.lat
    lon = itrs.earth_location.lon
    height = itrs.earth_location.height
    return lat.to(u.deg).value, lon.to(u.deg).value, height.to(u.km).value


# ----------------------
# Visualization preparation
# ----------------------
fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=ccrs.Robinson())
ax.set_global()
ax.coastlines()

colors = plt.cm.viridis(np.linspace(0, 1, num_steps))




satellite_data=[]
# Before the loop starts, set a phase offset for each sub-constellation
phase_offset_sub1 = 0
phase_offset_sub2 = 180/S

# Generate satellites independently for sub-constellation 1 and 2
for step, t in enumerate(time_steps, start=1):
    lats_all = []
    lons_all = []

    # ---------------------
    # Generate satellites for sub-constellation 1
    # ---------------------
    for p in range(num_planes_incl1):
        for s in range(S):
            # Calculate satellite position
            # For sub-constellation 1: pass num_planes=num_planes_incl1, default raan_offset=0.0
            orbit_t0 = get_orbit_parameters(p, s, time_steps[0], incl1, h1,
                                            delta_RAAN_deg_sub1, phase_offset_sub1,
                                            num_planes=num_planes_incl1)
            dt = (t - time_steps[0]).to(u.second)
            orbit = orbit_t0.propagate(dt)
            pos_eci = orbit.r
            lat, lon, h = eci_to_geodetic(pos_vec=pos_eci, obstime=t)
            lats_all.append(lat)
            lons_all.append(lon)

            # Determine sub-constellation ID, orbital plane, and satellite number
            sub_constellation_id = 1
            plane_id = p + 1
            sat_id = s + 1

            # Use formatted strings to generate IDs, ensuring two-digit plane and satellite numbers
            satellite_number = f"{sub_constellation_id}{plane_id:02d}{sat_id:02d}"

            # Store this satellite's data at this time step in the list
            satellite_data.append({
                "time": step,
                "satellite_id": satellite_number,
                "latitude_deg": lat,
                "longitude_deg": lon,
                "height_km": h
            })

    # ---------------------
    # Generate satellites for sub-constellation 2
    # ---------------------
    for p in range(num_planes_incl2):
        for s in range(S):
            # Calculate satellite position
            # For sub-constellation 2: pass num_planes=num_planes_incl2, add raan_offset=0.5
            orbit_t0 = get_orbit_parameters(p, s, time_steps[0], incl2, h2,
                                            delta_RAAN_deg_sub2, phase_offset_sub2,
                                            num_planes=num_planes_incl2, raan_offset=0.5)
            dt = (t - time_steps[0]).to(u.second)
            orbit = orbit_t0.propagate(dt)
            pos_eci = orbit.r
            lat, lon, h = eci_to_geodetic(pos_vec=pos_eci, obstime=t)
            lats_all.append(lat)
            lons_all.append(lon)

            # Determine sub-constellation ID, orbital plane, and satellite number
            sub_constellation_id = 2
            plane_id = p + 1
            sat_id = s + 1

            # Use formatted strings to generate IDs, ensuring two-digit plane and satellite numbers
            satellite_number = f"{sub_constellation_id}{plane_id:02d}{sat_id:02d}"

            # Store this satellite's data at this time step(snapshot) in the list
            satellite_data.append({
                "time": step,
                "satellite_id": satellite_number,
                "latitude_deg": lat,
                "longitude_deg": lon,
                "height_km": h
            })

#     # ---------------------
#     # Plot satellite positions for visualization and verification
#     # ---------------------
#     fig = plt.figure(figsize=(10, 5))
#     ax = plt.axes(projection=ccrs.Robinson())
#     ax.set_global()
#     ax.coastlines()
#
#     ax.scatter(lons_all, lats_all, transform=ccrs.PlateCarree(), s=20, color=colors[step - 1])
#
#     plt.title(f"Starlink D2C Constellation Snapshot {step}\nTime Step: {step}")
#     plt.savefig(os.path.join(output_dir, f'snapshot_{step}.png'), dpi=300, bbox_inches='tight')
#     plt.close(fig)
#
#     print(f"Snapshot {step} saved to {os.path.join(output_dir, f'snapshot_{step}.png')}")
#
# print(f"All snapshots saved to directory: {output_dir}")

# Export satellite data to CSV file
output_csv_path = os.path.join(output_dir, f'satellite_positions_cluster{cluster_config_mode}_config{constellation_config_mode}.csv')


with open(output_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ["time", "satellite_id", "latitude_deg", "longitude_deg", "height_km"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for row in satellite_data:
        writer.writerow(row)

print(f"Satellite IDs and 3D coordinate information saved to: {output_csv_path}")

# - satellite_data already contains position data for all satellites at all snapshots
# - center_hexagons list contains 37 cells, each with ['index'], ['center'], and computed ECEF coordinates


epsilon_min = 30.0  # elevation angle
Re = 6371.0


def deg2rad(deg):
    return deg * np.pi / 180


def phi_threshold(Re, h_sat, epsilon_min):
    # This formula establishes the relationship between the elevation angle (epsilon) and the satellite's Earth-centered cone angle (φ_th).
    # The geometric definition is at .../.../python/Geometric_Illustration.pdf
    return np.degrees(np.arcsin((Re / (Re + h_sat)) * np.cos(deg2rad(epsilon_min))))


def latlonh_to_ecef(lat_deg, lon_deg, h_km):
    # Convert geodetic coordinates to ECEF coordinates
    lat_rad = deg2rad(lat_deg)
    lon_rad = deg2rad(lon_deg)
    r = Re + h_km
    x = r * np.cos(lat_rad) * np.cos(lon_rad)
    y = r * np.cos(lat_rad) * np.sin(lon_rad)
    z = r * np.sin(lat_rad)
    return np.array([x, y, z])


# Set file path according to cluster_config_mode
if cluster_config_mode == 1:
    file_path = r'../../data/cell_clusters/output_cell_cluster1.csv'
elif cluster_config_mode == 2:
    file_path = r'../../data/cell_clusters/output_cell_cluster2.csv'
else:
    file_path = r'../../data/cell_clusters/output_cell_cluster3.csv'

center_hexagons = []
with open(file_path, newline='') as f:
    reader = csv.DictReader(f)
    # Only read records where region is 'center', and convert required fields to appropriate data types
    for row in reader:
        if row['region'] == 'center':
            hexagon = {}
            hexagon['index'] = int(row['index'])
            hexagon['center'] = [float(row['center_lat']), float(row['center_lon'])]
            hexagon['h3_index'] = row['h3_index']
            center_hexagons.append(hexagon)


# Add ECEF coordinates to center_hexagons (if not previously added)
for hexagon in center_hexagons:
    lat_c, lon_c = hexagon['center']
    hexagon['ecef'] = latlonh_to_ecef(lat_c, lon_c, 0.0)

from collections import defaultdict

time_groups = defaultdict(list)
for row in satellite_data:
    time_groups[row["time"]].append(row)

# Create visibility relationship output list
visibility_data = []

# Iterate over each snapshot
for current_time, sats in time_groups.items():
    for sat_info in sats:
        sat_id = sat_info["satellite_id"]
        lat_s = float(sat_info["latitude_deg"])
        lon_s = float(sat_info["longitude_deg"])
        h_s = float(sat_info["height_km"])

        # Calculate satellite's phi_threshold(threshold for visibility)
        curr_phi_th = phi_threshold(Re, h_s, epsilon_min)

        # Calculate satellite's ECEF coordinates
        S = latlonh_to_ecef(lat_s, lon_s, h_s)

        # Calculate Nadir direction vector (pointing toward Earth's center)
        N = -S  # Nadir direction is from satellite toward Earth's center

        visible_cells = []
        for hexagon in center_hexagons:
            C = hexagon['ecef']
            V = C - S  # Vector from satellite to cell center

            # Calculate angle theta
            dot_product = np.dot(N, V)
            norm_N = np.linalg.norm(N)
            norm_V = np.linalg.norm(V)
            # Prevent floating-point errors causing out-of-range values
            cos_theta = np.clip(dot_product / (norm_N * norm_V), -1.0, 1.0)
            angle = np.degrees(np.arccos(cos_theta))
            # Calculate cosine of angle between O->C and O->S
            dot_CS = np.dot(C, S)
            norm_C = np.linalg.norm(C)
            norm_S = np.linalg.norm(S)
            cos_alpha = dot_CS / (norm_C * norm_S)

            # Horizon condition: cos(horizon) = Re / |S|
            horizon_cos = Re / norm_S

            # Check if satellite is above the cell horizon (no Earth blockage)
            if cos_alpha >= horizon_cos:
                # Satellite is above cell horizon, further check phi_th condition
                if angle < curr_phi_th:
                    visible_cells.append(hexagon['index'])

        if len(visible_cells) > 0:
            for cell_id in visible_cells:
                # Convert visible cell list to comma-separated string
                cells_str = ",".join(map(str, visible_cells))
                visibility_data.append({
                    "time": current_time,
                    "satellite_id": sat_id,
                    "sat_lat_deg": lat_s,
                    "sat_lon_deg": lon_s,
                    "sat_height_km": h_s,
                    "cell_id": cell_id
                })

# Write visibility relationships to a new CSV file
visibility_csv_path = os.path.join(
    output_dir,
    f'visibility_relations_cluster{cluster_config_mode}_config{constellation_config_mode}.csv'
)
with open(visibility_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
    fieldnames = ["time", "satellite_id", "sat_lat_deg", "sat_lon_deg", "sat_height_km", "cell_id"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in visibility_data:
        writer.writerow(row)

print(f"Visibility relationship data saved to: {visibility_csv_path}")


center_lats = [cell['center'][0] for cell in center_hexagons]
center_lons = [cell['center'][1] for cell in center_hexagons]

# Read satellite position data
df_sat = pd.read_csv(output_csv_path)

# Read visibility relationship data
df_visibility = pd.read_csv(visibility_csv_path)

# Use Pandas to remove duplicate satellite-time pairs to get unique satellite positions per snapshot,if needed
df_visible_sat = df_visibility[['time', 'satellite_id']].drop_duplicates()

# Merge satellite position data to get coordinates of visible satellites
df_visible_sat = pd.merge(df_visible_sat, df_sat, on=['time', 'satellite_id'])

if cluster_config_mode == 1:
    center_lat = 46
    center_lng = 127

elif cluster_config_mode == 2:
    center_lat = -12
    center_lng = -77

elif cluster_config_mode == 3:
    center_lat = -33
    center_lng = -72
#       ------------- Visualization and Verification(Not Necessary) ---------------
# Define center coordinates for initial map location
map_center_lat = center_lat
map_center_lng = center_lng

# Color mapping, assign different colors for each snapshot
snapshot_colors = plt.cm.viridis(np.linspace(0, 1, num_steps))

# Flag to indicate whether legend label has been added
legend_added = False

# Iterate over each time step to generate corresponding map snapshots
for i, t in enumerate(time_steps, start=1):
    time_iso = t.iso  # Get ISO format time string

    # Get all visible satellite records at current time step
    sats_at_t = df_visible_sat[df_visible_sat['time'] == i]

    # Get nadir coordinates of these satellites
    lats_visible = sats_at_t['latitude_deg'].values
    lons_visible = sats_at_t['longitude_deg'].values

    # Create new map
    fig = plt.figure(figsize=(20, 10))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.coastlines()

    # Plot visible satellites' nadir points
    ax.scatter(
        lons_visible,
        lats_visible,
        transform=ccrs.PlateCarree(),
        s=10,
        color=snapshot_colors[i-1],
        label='Visible Satellites' if not legend_added else ""
    )

    # Plot center coordinates of 37 cells
    ax.scatter(
        center_lons,
        center_lats,
        transform=ccrs.PlateCarree(),
        s=1,  # Size adjustable as needed
        color='red',
        label='Center Cells' if not legend_added else ""
    )


    # Add title
    plt.title(f"Visibility Snapshot {i}\nTime: {time_iso}")

    # Add legend (only on first plot)
    if not legend_added:
        plt.legend(loc='lower left')
        legend_added = True

    # Save map as PNG image
    snapshot_path = os.path.join(output_dir, f'snapshot_visible_sat{i}.png')
    plt.savefig(snapshot_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

    print(f"Snapshot {i} saved to {snapshot_path}")
