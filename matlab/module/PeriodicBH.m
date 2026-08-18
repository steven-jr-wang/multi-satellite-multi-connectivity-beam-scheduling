function [Beam_pattern, C, R_over_D,J_fairness] = PeriodicBH(h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, sigma2, Ptx, w, Demand, wgs84, satellite_positions)

B = 30e6;

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
        warning(['Cell ', num2str(cell_m), ' has no visible satellite.']);
        continue;
    end
    [~, minIdx] = min(Load(epsilon_cell{cell_m}));
    n_star = epsilon_cell{cell_m}(minIdx);
    exchange_matrix(cell_m, n_star) = 1;

    Load(n_star) = Load(n_star) + Demand(cell_m);
end

% Build each satellite's assigned-cell set and periodic service order.
E     = cell(num_sats,1);
num_E = zeros(num_sats,1);
for n = 1:num_sats
    E{n}     = find(exchange_matrix(:,n)==1);
    num_E(n) = numel(E{n});
end
Q = max(1, floor(num_E/4));

pb_order = cell(num_sats,1);
pb_ptr   = ones(num_sats,1);
for n = 1:num_sats
    pb_order{n} = E{n}(:);
    if num_E(n)==0, pb_ptr(n)=1; end
end

Beam = cell(num_sats, Nvr);
Beam_pattern = cell(Nvr, 1);
Rate_matrix = zeros(num_cells, Nvr);
Rate_matrix_ideal = zeros(num_cells, Nvr); 
C = cell(Nvr, 1);

% Advance each satellite's persistent round-robin pointer every slot.
for t = 1:Nvr
    b_all = [];

    for n = 1:num_sats
        if num_E(n) == 0
            Beam{n,t} = [];
            continue;
        end
        q = min(Q(n), num_E(n));
        ord = pb_order{n};
        s  = pb_ptr(n);
        idxs = mod((s-1):(s+q-2), num_E(n)) + 1;
        Beam{n,t} = ord(idxs);
        pb_ptr(n) = mod(s+q-1, num_E(n)) + 1;

        b_all = [b_all; Beam{n,t}(:)];
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
        if ~isempty(n), C{t}(m, n) = 1; end
    end

    a_new = cell(num_cells, 1);
    for m = active_cells'
        epsilon_m = epsilon_cell{m};
        a_vec = zeros(length(epsilon_m), 1);
        for k = 1:length(epsilon_m)
            sat_idx = epsilon_m(k);
            a_vec(k) = exchange_matrix(m, sat_idx);
        end
        a_new{m} = a_vec;
    end

    for m = active_cells'
        a_m = a_new{m};
        numerator = dot(a_m.^2, h_vector_sqr{m, m});
        interference = 0;
        for j = active_cells'
            if j == m, continue; end
            a_j = a_new{j};
            h_mj = h_vector_sqr{m, j};
            if isempty(a_j) || isempty(h_mj), continue; end
            interference = interference + dot(a_j.^2, h_mj);
        end
        denom = interference + sigma2 / Ptx;
        gamma_var = numerator / denom;

        Rate_matrix_ideal(m, t) = w(m) * log2(1 + numerator / (sigma2 / Ptx));
        Rate_matrix(m, t)       = w(m) * log2(1 + gamma_var);
    end
end

% Compute demand satisfaction and Jain fairness.
mean_R_all = mean(Rate_matrix, 2);

selected_idx = find(Demand > 0);

R_over_D = zeros(length(selected_idx), 1);
for idx_roverd = 1:length(selected_idx)
    idx_selected_roverd = selected_idx(idx_roverd);
    R_over_D(idx_roverd) = mean_R_all(idx_selected_roverd) / Demand(idx_selected_roverd);
end

R_over_D_adjusted = min(R_over_D,ones(size(R_over_D)));

N_sel = length(selected_idx);
J_fairness = (sum(R_over_D_adjusted))^2 / (N_sel * sum(R_over_D_adjusted.^2));

end
