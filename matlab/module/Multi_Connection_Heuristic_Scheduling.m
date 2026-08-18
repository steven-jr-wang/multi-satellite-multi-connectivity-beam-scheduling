function [C_ultimate, R_over_D, J_fairness, mean_R_all] = Multi_Connection_Heuristic_Scheduling( ...
    h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, sigma2, ...
    sigma2FreqFrag, Ptx, w, Demand, wgs84)

% Derive each satellite's beam capacity from its visible cells.
exchange_matrix_visible = zeros(num_cells, num_sats);
for m = 1:num_cells
    visibleSats = epsilon_cell{m};
    if ~isempty(visibleSats)
        exchange_matrix_visible(m, visibleSats) = 1;
    end
end

Q = cell(num_sats, 1);
for q = 1:num_sats
    numVisible = sum(exchange_matrix_visible(:, q));
    Q{q} = max(1, floor(numVisible / 2));
end

Nvr = 100;

C_ultimate = cell(Nvr, 1);
Rate_matrix_ultimate = zeros(num_cells, Nvr);

active_cells = 1:num_cells;

% Greedily retain links that increase the network sum rate.
for t = 1:Nvr

    start_m = randi(num_cells);
    ordered_cells = [start_m:num_cells, 1:start_m-1];

    sat_load = zeros(num_sats, 1);
    selected_cells_per_satellite = cell(num_sats, 1);

    exchange_matrix = zeros(num_cells, num_sats);

    for i = 1:length(ordered_cells)
        m = ordered_cells(i);
        epsilon_m = epsilon_cell{m};

        [~, Rate_vec, sumRate_cur] = computeSumRate_givenExchange( ...
            exchange_matrix, active_cells, epsilon_cell, h_vector_sqr, sigma2FreqFrag, Ptx, w, num_cells);

        for k = 1 : length(epsilon_m)
            n = epsilon_m(k);

            if sat_load(n) >= Q{n}
                continue;
            end

            exchange_matrix(m, n) = 1;

            [a_tmp, Rate_tmp, sumRate_tmp] = computeSumRate_givenExchange( ...
                exchange_matrix, active_cells, epsilon_cell, h_vector_sqr, sigma2FreqFrag, Ptx, w, num_cells);

            if sumRate_tmp > sumRate_cur + 1e-12
                sumRate_cur = sumRate_tmp;
                a_new = a_tmp;
                Rate_vec = Rate_tmp;

                sat_load(n) = sat_load(n) + 1;
                selected_cells_per_satellite{n} = [selected_cells_per_satellite{n}; m];
            else

                exchange_matrix(m, n) = 0;
            end
        end
    end

    C_ultimate{t} = exchange_matrix;
    Rate_matrix_ultimate(:, t) = Rate_vec;

end

% Compute per-cell throughput, demand satisfaction, and fairness.
slotLoss = 6/7;
mean_R_all = mean(Rate_matrix_ultimate, 2) * slotLoss;

selected_idx = find(Demand > 0);
R_over_D = zeros(length(selected_idx), 1);
for idx_roverd = 1:length(selected_idx)
    m = selected_idx(idx_roverd);
    R_over_D(idx_roverd) = mean_R_all(m) / Demand(m);
end

R_over_D_adjusted = min(R_over_D, ones(size(R_over_D)));
N_sel = length(selected_idx);
J_fairness = (sum(R_over_D_adjusted))^2 / (N_sel * sum(R_over_D_adjusted.^2) + 1e-12);

end

% Evaluate rates for a candidate exchange matrix.
function [a_new, Rate_vec, sumRate] = computeSumRate_givenExchange( ...
    exchange_matrix, active_cells, epsilon_cell, h_vector_sqr, sigma2FreqFrag, Ptx, w, num_cells)

a_new = cell(num_cells, 1);
Rate_vec = zeros(num_cells, 1);

for m = active_cells(:).'
    epsilon_m = epsilon_cell{m};
    a_new_pattern = zeros(length(epsilon_m), 1);
    for idx = 1:length(epsilon_m)
        sat_idx = epsilon_m(idx);
        a_new_pattern(idx) = exchange_matrix(m, sat_idx);
    end
    a_new{m} = a_new_pattern;
end

for m = active_cells(:).'
    a_m = a_new{m};
    numerator_m = dot(a_m.^2, h_vector_sqr{m, m});

    interference_m = 0;
    for j = active_cells(:).'
        if j ~= m
            a_j_vec = a_new{j};
            h_j_vec = h_vector_sqr{m, j};
            if isempty(a_j_vec) || isempty(h_j_vec)
                continue;
            end
            interference_m = interference_m + dot(a_j_vec.^2, h_j_vec);
        end
    end

    denominator_m = interference_m + sigma2FreqFrag(m) / Ptx;
    gamma_m = numerator_m / max(denominator_m, 1e-15);

    Rate_vec(m) = w(m) * log2(1 + gamma_m);
end

sumRate = sum(Rate_vec(active_cells));

end
