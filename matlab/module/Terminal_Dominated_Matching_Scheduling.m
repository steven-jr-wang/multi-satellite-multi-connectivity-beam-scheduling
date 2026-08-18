function [C_ultimate, R_over_D, J_fairness, mean_R_all] = Terminal_Dominated_Matching_Scheduling( ...
    h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, sigma2, ...
    sigma2FreqFrag, Ptx, w, Demand, wgs84)

% Derive each satellite's proposal capacity from its visible cells.
exchange_matrix_visible = zeros(num_cells, num_sats);
for m = 1:num_cells
    visibleSats = epsilon_cell{m};
    if ~isempty(visibleSats)
        exchange_matrix_visible(m, visibleSats) = 1;
    end
end

Q = cell(num_sats, 1);
for n = 1:num_sats
    numVisible = sum(exchange_matrix_visible(:, n));
    Q{n} = max(1, floor(numVisible / 2));
end

Nvr = 100;
I_max = 2;

C_ultimate = cell(Nvr, 1);
Rate_matrix_ultimate = zeros(num_cells, Nvr);

active_cells = 1:num_cells;

% Match satellite proposals to each cell's strongest visible links.
for t = 1:Nvr

    proposal_list = cell(num_cells, 1);
    for m = 1:num_cells
        proposal_list{m} = [];
    end

    for n = 1:num_sats
        visible_cells_n = find(exchange_matrix_visible(:, n) == 1);

        [~, idx_sort] = sort(Demand(visible_cells_n), 'descend');
        sorted_cells = visible_cells_n(idx_sort);

        num_to_propose = min(Q{n}, length(sorted_cells));
        proposed_cells = sorted_cells(1:num_to_propose);

        for k = 1:num_to_propose
            m = proposed_cells(k);
            proposal_list{m} = [proposal_list{m}, n];
        end
    end

    exchange_matrix = zeros(num_cells, num_sats);

    start_m = randi(num_cells);
    ordered_cells = [start_m:num_cells, 1:start_m-1];

    for m = ordered_cells
        sats_req = proposal_list{m};

        if isempty(sats_req)
            continue;
        end

        if length(sats_req) <= I_max
            exchange_matrix(m, sats_req) = 1;
        else

            epsilon_m = epsilon_cell{m};
            h_mm = h_vector_sqr{m, m};

            [tf, loc] = ismember(sats_req, epsilon_m);

            loc = loc(tf);

            h_quality = h_mm(loc);
            [~, idx_best] = sort(h_quality, 'descend');

            chosen_loc = loc(idx_best(1:I_max));
            chosen_sats = epsilon_m(chosen_loc);

            exchange_matrix(m, chosen_sats) = 1;
        end
    end

    [~, Rate_vec, ~] = computeSumRate_givenExchange( ...
        exchange_matrix, active_cells, epsilon_cell, h_vector_sqr, sigma2FreqFrag, Ptx, w, num_cells);

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

% Evaluate rates for the accepted exchange matrix.
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
