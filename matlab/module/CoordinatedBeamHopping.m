function [Beam_pattern, C, R_over_D,J_fairness, mean_R_all] = CoordinatedBeamHopping(h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, sigma2, Ptx, ~, Demand, wgs84, ...
    nbrOfRB)
B=nbrOfRB*180*1e3;
w = B * ones(num_cells, 1)/1e6;

Nvr = 100;

[~, sortedIdx] = sort(Demand, 'descend');
sortedCells = sortedIdx;

Load = zeros(num_sats, 1);

% Assign each cell to the least-loaded visible satellite.
exchange_matrix = zeros(num_cells, num_sats);

for m = 1:num_cells
    cell_m = sortedCells(m);
    if Demand(cell_m) <= 0
        continue;
    end

    if isempty(epsilon_cell{cell_m})
        warning(['Cell ', num2str(cell_m), ' has no visible satellite and cannot be assigned.']);
        continue;
    end

    [~, minIdx] = min(Load(epsilon_cell{cell_m}));
    n_star = epsilon_cell{cell_m}(minIdx);

    exchange_matrix(cell_m, n_star) = 1;

    Load(n_star) = Load(n_star) + Demand(cell_m);
end

% Build the six-nearest-cell adjacency matrix.
global_adj_matrix = zeros(num_cells, num_cells);

for m = 1:num_cells

    m_lat = center_data.center_lat(m);
    m_lon = center_data.center_lon(m);
    m_alt = 0;

    distances = zeros(num_cells, 1);
    [mX, mY, mZ] = geodetic2ecef(wgs84, m_lat, m_lon, m_alt, "degrees");
    for j = 1:height(data)
        if m == data.index(j)
            distances(j) = inf;
        else
            [jX, jY, jZ] = geodetic2ecef(wgs84, data.center_lat(j), data.center_lon(j), 0, "degrees");
            distances(j) = sqrt((mX - jX)^2 + (mY - jY)^2 + (mZ - jZ)^2);
        end
    end

    [~, sorted_indices] = sort(distances, 'ascend');
    nearest_six = sorted_indices(1:6);

    nearest_six_indices = data.index(nearest_six);

    global_adj_matrix(m, nearest_six_indices) = 1;
end

addpath('YOUR_GUROBI_MATLAB_PATH');
gurobi_setup

% Generate feasible beam patterns for each satellite.
phi_n_cell = cell(num_sats, 1);
V_Pn_cell = cell(num_sats, 1);
xi_n_values = zeros(num_sats, 1);
E=cell(num_sats,1);
num_E=zeros(num_sats,1);
for n = 1:num_sats

    E{n} = find(exchange_matrix(:, n) == 1);
    num_E(n) = length(E{n});
    if num_E(n) == 0
        phi_n_cell{n} = [];
        V_Pn_cell{n} = {};
        xi_n_values(n) = 0;
        continue;
    end

    Adj = zeros(num_E(n), num_E(n));
    for i = 1:num_E(n)
        for j = 1:num_E(n)
            if i ~= j

                cell_i = E{n}(i);
                cell_j = E{n}(j);

                Adj(i,j) = global_adj_matrix(cell_i, cell_j);
            end
        end
    end

    Q = 1;
    if num_E(n) > 4 * Q
        Q = floor(num_E(n)/ 4);
    end

    V_Pn = {};
    combinations = nchoosek(1:num_E(n), Q);
    numComb = size(combinations, 1);
    for k = 1:numComb
        z_k = zeros(1, num_E(n));
        z_k(combinations(k, :)) = 1;

        if (z_k * Adj * z_k') == 0
            V_Pn{end+1} = z_k;
        end
    end

    if isempty(V_Pn)
        phi_n_cell{n} = [];
        V_Pn_cell{n} = {};
        xi_n_values(n) = 0;
        continue;
    end

    rate_vectors = cell(length(V_Pn), 1);
    for k = 1:length(V_Pn)
        z_k = V_Pn{k};
        r_est = zeros(num_E(n), 1);
        for j = 1:num_E(n)
            if z_k(j) == 1

                cell_j = E{n}(j);

                satellite_list = epsilon_cell{cell_j};
                k_sat = find(satellite_list == n, 1);

                if ~isempty(k_sat)
                    h_nj = h_vector{cell_j, cell_j}(k_sat);
                else
                    warning(['Satellite ', num2str(n), ' is not visible to cell ', num2str(cell_j), '.']);
                    h_nj = 0;
                end

                gamma = (Ptx * abs(h_nj)^2) / sigma2;
                r_est(j) = B/(1e6) * log2(1 + gamma);

            else

                r_est(j) = 0;
            end
        end
        rate_vectors{k} = r_est;
    end

    % Optimize the number of uses of each single-satellite pattern.
    model = struct();
    model.modelsense = 'max';

    num_phi = length(V_Pn);
    num_vars = num_phi + 1;

    phi_varnames = strcat('phi_', arrayfun(@num2str, 1:num_phi, 'UniformOutput', false), '_n_', num2str(n));

    xi_varname = ['xi_n_', num2str(n)];

    model.varnames = [phi_varnames, {xi_varname}];

    model.vtype = [repmat('I', num_phi, 1); 'C'];

    model.lb = [zeros(num_phi, 1); 0];

    model.ub = [ones(num_phi, 1) * Nvr; inf];

    model.obj = [zeros(num_phi, 1); 1];

    Aeq = [ones(1, num_phi), 0];
    beq = Nvr;

    num_cells_n = num_E(n);
    A_ineq = zeros(num_cells_n, num_vars);
    b_ineq = zeros(num_cells_n, 1);

    for j = 1:num_cells_n
        Demand_j = Demand(E{n}(j));

        for k = 1:num_phi
            r_kj = rate_vectors{k}(j);
            A_ineq(j, k) = (r_kj) / Nvr;
        end
        A_ineq(j, end) = -Demand_j;
        b_ineq(j) = 0;
    end

    model.A = sparse([Aeq; A_ineq]);
    model.rhs = [beq; b_ineq];
    model.sense = [repmat('=', size(Aeq, 1), 1); repmat('>', size(A_ineq, 1), 1)];
    if size(model.A, 1) ~= length(model.rhs)
        error('The constraint matrix and right-hand side have different row counts.');
    end

    params = struct();
    params.OutputFlag = 0;
    params.IterationLimit = 100000;

    try
        result = gurobi(model, params);
    catch ME
        warning(['Gurobi failed while designing beam patterns for satellite ', num2str(n), ': ', ME.message]);
        result = [];
    end

    if ~isfield(result, 'x') || isempty(result.x)
        warning(['Gurobi found no solution for satellite ', num2str(n), ' because its cells have no demand.']);
        phi_n = zeros(num_phi, 1);
        xi_n = 0;
    else
        phi_n = result.x(1:num_phi);
        xi_n = result.x(end);
    end

    phi_n_cell{n} = phi_n;
    V_Pn_cell{n} = V_Pn;
    xi_n_values(n) = xi_n;
end

% Coordinate satellite patterns by maximizing inter-pattern separation.
Beam = cell(num_sats, Nvr);

Beam_pattern = cell(Nvr, 1);

remaining_phi_n = phi_n_cell;
V_Pn_all = V_Pn_cell;

[~, sortedSatIdx] = sort(Load, 'ascend');
sortedSats = sortedSatIdx;

Rate_matrix = zeros(num_cells, Nvr);
Rate_matrix_ideal = zeros(num_cells, Nvr);

C = cell(Nvr, 1);

for t = 1:Nvr

    b_all = [];

    for idx = 1:num_sats
        n = sortedSats(idx);

        available_patterns = find(remaining_phi_n{n} > 0);
        b_all=b_all(:);

        if isempty(available_patterns)

            selected_pattern = [];

            Beam{n, t} = selected_pattern;
            continue;
        end

        if idx == 1

            pattern_idx = available_patterns(randi(length(available_patterns)));
            selected_pattern = V_Pn_all{n}{pattern_idx};

            active_cells_n = E{n}(selected_pattern == 1);

            remaining_phi_n{n}(pattern_idx) = remaining_phi_n{n}(pattern_idx) - 1;

            if ~isempty(active_cells_n)
                Beam{n, t} = active_cells_n;
                b_all = [b_all; active_cells_n];
            else
                Beam{n, t} = [];
            end
        else

            best_min_dist = -inf;
            best_pattern_idx = -1;
            best_pattern = [];

            for pattern_idx = available_patterns'
                candidate_pattern = V_Pn_all{n}{pattern_idx};

                active_cells_n = E{n}(candidate_pattern == 1);
                active_cells_n=active_cells_n(:);

                if isempty(b_all)
                    min_distance = inf;
                else
                    min_distance = inf;

                    for cell_a = active_cells_n'

                        p_lat = center_data.center_lat(cell_a);
                        p_lon = center_data.center_lon(cell_a);

                        for cell_b = b_all'

                            q_cell_instances = data(data.index == cell_b, :);
                            if isempty(q_cell_instances)
                                continue;
                            end

                            distances = distance([p_lat, p_lon], [q_cell_instances.center_lat, q_cell_instances.center_lon]);

                            current_min_dist_candidate = min(distances);

                            if current_min_dist_candidate < min_distance
                                min_distance = current_min_dist_candidate;
                            end
                        end
                    end
                end

                if min_distance > best_min_dist
                    best_min_dist = min_distance;
                    best_pattern_idx = pattern_idx;
                    best_pattern = candidate_pattern;
                end
            end

            if best_pattern_idx ~= -1
                selected_pattern = best_pattern;

                active_cells_n = E{n}(selected_pattern == 1);

                remaining_phi_n{n}(best_pattern_idx) = remaining_phi_n{n}(best_pattern_idx) - 1;
                if ~isempty(active_cells_n)
                    Beam{n, t} = active_cells_n;
                    b_all = [b_all; active_cells_n];
                else
                    Beam{n, t} = [];
                end
            else

                selected_pattern = [];
            end
        end
    end

    if ~isempty(b_all)
        Beam_pattern{t} = sort(b_all);
    else
        Beam_pattern{t} = [];
    end

    active_cells = Beam_pattern{t};

    C{t} = zeros(num_cells, num_sats);

    for m = active_cells'

        n = find(exchange_matrix(m, :) == 1, 1);
        if ~isempty(n)
            C{t}(m, n) = 1;
        else
            warning(['Cell ', num2str(m), ' has no connected satellite in slot ', num2str(t), '.']);
        end
    end

    a_new = cell(num_cells, 1);

    for m = active_cells'
        epsilon_m = epsilon_cell{m};
        a_new_pattern = zeros(length(epsilon_m), 1);
        for idx = 1:length(epsilon_m)
            sat_idx = epsilon_m(idx);
            a_value = exchange_matrix(m, sat_idx);
            a_new_pattern(idx) = a_value;
        end
        a_new{m} = a_new_pattern;
    end

    for m = active_cells'
        a_m = a_new{m};
        numerator = dot(a_m.^2, h_vector_sqr{m, m});
        interference = 0;

        for j = active_cells'
            if j ~= m
                a_j_vec = a_new{j};
                h_j_vec = h_vector_sqr{m, j};
                if isempty(a_j_vec) || isempty(h_j_vec)
                    continue;
                end
                interference = interference + dot(a_j_vec.^2, h_j_vec);
            end
        end

        denominator = interference + sigma2 / Ptx;
        gamma_var = numerator / denominator;

        Rate_matrix_ideal(m, t) = w(m) * log2(1 + numerator / (sigma2 / Ptx));
        Rate_matrix(m, t) = w(m) * log2(1 + gamma_var);
    end

end

% Compute per-cell throughput, demand satisfaction, and fairness.
slotLoss = 6/7;
mean_R_all = mean(Rate_matrix, 2) * slotLoss;

selected_idx = find(Demand > 0);

R_over_D = zeros(length(selected_idx), 1);
for idx_roverd = 1:length(selected_idx)
    idx_selected_roverd = selected_idx(idx_roverd);
    R_over_D(idx_roverd) = mean_R_all(idx_selected_roverd) / Demand(idx_selected_roverd);
end

R_over_D_adjusted = min(R_over_D,ones(size(R_over_D)));
N_sel = length(selected_idx);
J_fairness = (sum(R_over_D_adjusted))^2 / (N_sel * sum(R_over_D_adjusted.^2));

    % Return the three-dimensional separation in kilometers.
    function d = distance(coord1, coord2)

        wgs84 = wgs84Ellipsoid('kilometer');
        [x1, y1, z1] = geodetic2ecef(wgs84, coord1(1), coord1(2), 0, "degrees");
        [x2, y2, z2] = geodetic2ecef(wgs84, coord2(1), coord2(2), 0, "degrees");
        d = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    end

end
