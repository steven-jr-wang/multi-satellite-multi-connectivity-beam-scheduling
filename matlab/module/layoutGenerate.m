function [h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users]  = layoutGenerate(selected_demand, csvFile, csvFile2, csvFile3,nbrOfUserPerCell,nbrOfRB, relDensityParams)

data = readtable(csvFile);
data2 = readtable(csvFile2);
% Select one visibility snapshot and remap satellite identifiers.
visibility_data = readtable(csvFile3);

unique_times = unique(visibility_data.time);

selected_time = unique_times(randi(length(unique_times)));

visibility_at_t = visibility_data(visibility_data.time == selected_time, :);

original_sats = unique(visibility_at_t.satellite_id);
num_sats = length(original_sats);

[~, loc] = ismember(visibility_at_t.satellite_id, original_sats);

visibility_at_t.SatelliteIndex = loc;

% Load center-cell geometry and generate the traffic density profile.
center_data = data(strcmp(data.region, 'center'), :);
num_cells = height(center_data);

if exist('relDensityParams','var') && relDensityParams.useGenerated

    switch relDensityParams.distType
        case 0
            relative_density = ones(1, num_cells) / num_cells;
        case 1
            values = gamrnd(relDensityParams.alpha, 1, 1, num_cells);
            relative_density = values / sum(values);
        case 2
            values = lognrnd(relDensityParams.mu, relDensityParams.sigma, 1, num_cells);
            relative_density = values / sum(values);
        case 3
            U = rand(1, num_cells);
            values = relDensityParams.xm * (1 - U).^(-1/relDensityParams.alpha);
            relative_density = values / sum(values);
    end
else

    relative_density = data2.relative_density;
end

total_users = height(center_data) * nbrOfUserPerCell;

total_density = sum(relative_density);
mean_users_per_cell = floor((relative_density / total_density) * total_users);

mu_mean_users = mean(mean_users_per_cell);

var_mean_users = var(mean_users_per_cell);

norm_var_users = var_mean_users / (mu_mean_users^2);

users_per_cell = poissrnd(mean_users_per_cell);
maxRBperCell = floor((nbrOfRB + 1)/2);
users_per_cell = min(users_per_cell, maxRBperCell);
users_per_cell(users_per_cell < 0) = 0;

% Sample users and their traffic demands within each center cell.
Full_Table = [];

user_index = 1;
for i = 1:height(center_data)
    num_users = users_per_cell(i);
    if num_users > 0

        vertex_lats = [center_data.vertex1_lat(i), center_data.vertex2_lat(i), center_data.vertex3_lat(i), ...
            center_data.vertex4_lat(i), center_data.vertex5_lat(i), center_data.vertex6_lat(i)];
        vertex_lons = [center_data.vertex1_lon(i), center_data.vertex2_lon(i), center_data.vertex3_lon(i), ...
            center_data.vertex4_lon(i), center_data.vertex5_lon(i), center_data.vertex6_lon(i)];

        for j = 1:num_users

            [lat, lon] = generate_random_point_in_polygon(vertex_lats, vertex_lons);

            demand_user_select = selected_demand;
            Full_Table = [Full_Table; user_index, i, lat, lon, demand_user_select];
            user_index = user_index + 1;
        end
    end
end

cell_demand = zeros(height(center_data), 1);

for i = 1:height(center_data)
    cell_demand(i) = sum(Full_Table(Full_Table(:, 2) == i, 5));
end

center_demand_table = array2table([(1:height(center_data))', cell_demand], 'VariableNames', {'Index', 'DemandMbps'});

region_list = ["center", "wrap_1", "wrap_2", "wrap_3", "wrap_4", "wrap_5", "wrap_6"];
demand_table = [];

for region = region_list
    temp_table = center_demand_table;
    temp_table.Region = repmat(region, height(center_demand_table), 1);
    temp_table = [temp_table(:, 3), temp_table(:, 1:2)];
    demand_table = [demand_table; temp_table];
end

Full_Table=array2table(Full_Table, 'VariableNames', {'UserIndex', 'CellIndex', 'Latitude', 'Longitude', 'DemandMbps'});

repUE_lat = center_data.center_lat;
repUE_lon = center_data.center_lon;

unique_sats_idx = unique(visibility_at_t.SatelliteIndex);

% Build satellite position and visibility tables for the selected snapshot.
satellite_positions = table();
satellite_positions.SatelliteIndex = unique_sats_idx;
satellite_positions.Latitude = zeros(num_sats,1);
satellite_positions.Longitude = zeros(num_sats,1);
satellite_positions.Altitude = zeros(num_sats,1);

for k = 1:num_sats

    sat_rows = visibility_at_t(visibility_at_t.SatelliteIndex == k, :);

    satellite_positions.Latitude(k) = sat_rows.sat_lat_deg(1);
    satellite_positions.Longitude(k) = sat_rows.sat_lon_deg(1);
    satellite_positions.Altitude(k) = sat_rows.sat_height_km(1);
end

satellite_fov_table = table();
satellite_fov_table.SatelliteIndex = visibility_at_t.SatelliteIndex;
satellite_fov_table.CellIndex = visibility_at_t.cell_id;

cell_to_sat_table = table([], [], [], [], [], 'VariableNames', {'CellIndex', 'SatelliteIndex', 'Latitude', 'Longitude', 'Altitude'});

for i = 1:height(center_data)

    rel_sats = satellite_fov_table(satellite_fov_table.CellIndex == i, :);

    for j = 1:height(rel_sats)
        sat_idx = rel_sats.SatelliteIndex(j);

        sat_pos = satellite_positions(satellite_positions.SatelliteIndex == sat_idx, :);
        new_row = {i, sat_idx, sat_pos.Latitude, sat_pos.Longitude, sat_pos.Altitude};
        cell_to_sat_table = [cell_to_sat_table; new_row];
    end
end

cell_to_sat_table = unique(cell_to_sat_table, 'rows');

nbrOfRB_disabled = users_per_cell - 1;
nbrOfRB_disabled(nbrOfRB_disabled < 0) = 0;

% Configure bandwidth, thermal noise, and transmit power.
B=nbrOfRB*180*1e3;
w=B*(1-nbrOfRB_disabled./nbrOfRB)/1e6;
kT = -174;
sigma2_dBm = kT + 10*log10(B);
sigma2 = 10^((sigma2_dBm - 30)/10);

sigma2FreqFrag_dBm = kT + 10*log10(w * 1e6);
sigma2FreqFrag = 10.^((sigma2FreqFrag_dBm - 30)/10);

EIRP_density = 34;
EIRP = EIRP_density + 10*log10(B/1e6);
Gt = 30;
Pt_dB = EIRP - Gt;
Ptx = 10^(Pt_dB/10);

wgs84 = wgs84Ellipsoid('kilometer');

Demand = zeros(num_cells, 1);
for i = 1:num_cells
    Demand(i) = demand_table.DemandMbps(i);
end

% Build visible-satellite sets and cross-cell channel matrices.
epsilon_cell = cell(num_cells, 1);
h_vector = cell(num_cells,num_cells);
h_vector_sqr = cell(num_cells,num_cells);

for m = 1:num_cells
    epsilon_cell{m} = [];

    visible_sats = cell_to_sat_table(cell_to_sat_table.CellIndex == m, :);
    for idx = 1:height(visible_sats)
        sat_idx=visible_sats.SatelliteIndex(idx);
        epsilon_cell{m} = [epsilon_cell{m}, sat_idx];
    end
end

for m = 1:num_cells
    for j = 1:num_cells
        h_vector{m, j} = [];
        h_vector_sqr{m, j} = [];

        epsilon_j = epsilon_cell{j};

        h_mj_vec = zeros(length(epsilon_j), 1);
        h_mj_vec_sqr = zeros(length(epsilon_j), 1);

        for idx = 1:length(epsilon_j)
            sat_idx = epsilon_j(idx);
            sat_lat = satellite_positions.Latitude(sat_idx);
            sat_lon = satellite_positions.Longitude(sat_idx);
            sat_alt = satellite_positions.Altitude(sat_idx);
            users_in_m = Full_Table(Full_Table.CellIndex == m, :);
            num_users_in_m = size(users_in_m, 1);
            avg_h = 0;
            if num_users_in_m > 0
                for user_idx = 1:height(users_in_m)
                    user_lat=users_in_m.Latitude(user_idx);
                    user_lon=users_in_m.Longitude(user_idx);

                    cells_j = data(data.index == j, :);
                    min_dist = inf;
                    nearest_cell_lat = 0;
                    nearest_cell_lon = 0;

                    for c = 1:height(cells_j)
                        cell_lat = cells_j.center_lat(c);
                        cell_lon = cells_j.center_lon(c);
                        cell_alt = 0;

                        [cellX, cellY, cellZ] = geodetic2ecef(wgs84, cell_lat, cell_lon, cell_alt, "degrees");
                        [ueX, ueY, ueZ] = geodetic2ecef(wgs84, user_lat, user_lon, 0, "degrees");

                        dist = norm([ueX - cellX, ueY - cellY, ueZ - cellZ]);

                        if dist < min_dist
                            min_dist = dist;
                            nearest_cell_lat = cell_lat;
                            nearest_cell_lon = cell_lon;
                        end
                    end

                    h = Channel_Generate(sat_lat, sat_lon, sat_alt,user_lat,user_lon,0,nearest_cell_lat,nearest_cell_lon,0);
                    avg_h = avg_h+h/num_users_in_m;
                end

            else

                user_lat = repUE_lat(m);
                user_lon = repUE_lon(m);

                cells_j = data(data.index == j, :);
                min_dist = inf;
                nearest_cell_lat = 0;
                nearest_cell_lon = 0;

                for c = 1:height(cells_j)
                    cell_lat = cells_j.center_lat(c);
                    cell_lon = cells_j.center_lon(c);
                    cell_alt = 0;

                    [cellX, cellY, cellZ] = geodetic2ecef(wgs84, cell_lat, cell_lon, cell_alt, "degrees");
                    [ueX, ueY, ueZ] = geodetic2ecef(wgs84, user_lat, user_lon, 0, "degrees");
                    dist = norm([ueX - cellX, ueY - cellY, ueZ - cellZ]);

                    if dist < min_dist
                        min_dist = dist;
                        nearest_cell_lat = cell_lat;
                        nearest_cell_lon = cell_lon;
                    end
                end

                avg_h = Channel_Generate(sat_lat, sat_lon, sat_alt, ...
                    user_lat, user_lon, 0, ...
                    nearest_cell_lat, nearest_cell_lon, 0);
            end

            h_mj_vec(idx) = avg_h;

            h_mj_vec_sqr(idx) = avg_h^2;
        end

        h_vector{m, j} = h_mj_vec;

        h_vector_sqr{m, j} = h_mj_vec_sqr;
    end
end

end

