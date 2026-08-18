function [C_ultimate, R_over_D,J_fairness,mean_R_all]=proposedMethod(h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
    ~,sigma2FreqFrag, Ptx, w, Demand, wgs84)
  % Configure the per-satellite beam limits.
  exchange_matrix_visible = zeros(num_cells, num_sats);
  for m = 1:num_cells
      visibleSats=epsilon_cell{m};
      exchange_matrix_visible(m,visibleSats)=1;
  end
  Q=cell(num_sats,1);
  for q = 1:num_sats
      numVisible = sum(exchange_matrix_visible(:, q));
      Q{q} = max(1, floor(numVisible / 2));
  end



% Build first-, second-, and third-ring cell neighborhoods.
adjacentList = cell(num_cells, 1);
adjacentList2 = cell(num_cells, 1);
adjacentList3 = cell(num_cells, 1);
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
    neighbors = data.index(nearest_six);
    adjacentList{m} = neighbors;
    ring2_indices = sorted_indices(7:19);
    ring2_neighbors = data.index(ring2_indices);
    adjacentList2{m} = ring2_neighbors;
    ring3_indices = sorted_indices(20:37);
    ring3_neighbors = data.index(ring3_indices);
    adjacentList3{m} = ring3_neighbors;
end
disp('Cell neighborhood lists generated.');
% Partition cells into demand-balanced groups.
clusterAssignments = zeros(num_cells, 1);
[~, sortedIdx] = sort(Demand, 'descend');
group_size = ceil(num_cells / 3);
delta = ceil(group_size / 2);
start_ptr = 1;
end_ptr = num_cells;
group_num = 1;
assigned = false(num_cells, 1);
while start_ptr <= end_ptr
    current_group = [];
    
    front_candidates = sortedIdx(start_ptr:min(start_ptr + delta - 1, num_cells));
    front_candidates = front_candidates(~assigned(front_candidates));
    current_group = [current_group; front_candidates];
    assigned(front_candidates) = true;
    
    back_candidates = sortedIdx(max(end_ptr - delta + 1, start_ptr):end_ptr);
    back_candidates = back_candidates(~assigned(back_candidates));
    current_group = [current_group; back_candidates];
    assigned(back_candidates) = true;
    
    clusterAssignments(current_group) = group_num;
    group_num = group_num + 1;
    
    start_ptr = start_ptr + delta;
    end_ptr = end_ptr - delta;
end
if any(~assigned)
    remaining_cells = ~assigned;
    clusterAssignments(remaining_cells) = group_num - 1;
end
group_num = group_num - 1;
disp('Demand-based slot grouping completed.');
% Initialize cell-to-satellite association variables.
a = cell(num_cells, 1);
for m = 1:num_cells
    a{m} = zeros(length(epsilon_cell{m}), 1);
end
a_new = cell(num_cells, 1);
a_var = cell(num_cells, 1);
groups = cell(group_num, 1);
for c = 1:group_num
    groups{c} = find(clusterAssignments == c);
end
numerator=zeros(num_cells, 1);
denominator = zeros(num_cells, 1);
numerator_new=zeros(num_cells, 1);
gamma = zeros(num_cells, 1);
y=zeros(num_cells, 1);
d=cell(num_cells, 1);
J=cell(num_cells, 1);
gamma_var= zeros(num_cells, 1);
R_prev= zeros(num_cells, 1);
S_prev= zeros(num_cells, 1);
R_update= zeros(num_cells, 1);
S_update= zeros(num_cells, 1);
delta_S= zeros(num_cells, 1);
delta_S_prev= zeros(num_cells, 1);
alpha = zeros(num_cells, 1);
beta = zeros(num_cells, 1);
mu = zeros(num_cells, 1);
nu = zeros(num_cells, 1);
thresholds = [1;4];
eta_threshold = [0.99, 0.99];
for m = 1:num_cells
    d{m} = zeros(length(epsilon_cell{m}), 1);
    J{m} = zeros(length(epsilon_cell{m}), 1);
end
function patterns = recursivePatternGeneration(num_cells, active_cells, potential_active_cells, a_var, adjacentList, adjacentList2, Demand, epsilon_cell, h_vector_sqr, h_vector, w, gamma, y, sigma2FreqFrag, Ptx, max_depth,thresholds,depth, d, J,...
       eta_threshold)
    % Recursively enumerate feasible routing patterns.
    patterns = {};
    if (isempty(potential_active_cells)) || (max_depth == 0)

        sc = active_cells(1);
        Rate_comparison = w(sc) * log2(1 + gamma(sc));
        eta_value = Rate_comparison / Demand(sc);
        pattern.active_cells = active_cells;
        pattern.a_var = a_var;
        pattern.sign_cell = active_cells(1);
        pattern.eta_value = eta_value;
        patterns{end+1} = pattern;
        return;
    end
 
    for cell_idx = potential_active_cells
        backup_gamma = zeros(num_cells,1);
        backup_y     = zeros(num_cells,1);
        backup_d     = cell(num_cells,1);
        backup_J     = cell(num_cells,1);
        for cID = active_cells
            backup_gamma(cID) = gamma(cID);
            backup_y(cID)     = y(cID);
            backup_d{cID}     = d{cID};
            backup_J{cID}     = J{cID};
        end

        [a_new_cell, feasible] = assignBeams(cell_idx, a_var, Demand, epsilon_cell, h_vector_sqr, h_vector, w, y, sigma2FreqFrag, Ptx, active_cells, eta_threshold);
        
        if ~feasible
            for cID = [active_cells, cell_idx]
                gamma(cID) = backup_gamma(cID);
                y(cID)     = backup_y(cID);
                d{cID}     = backup_d{cID};
                J{cID}     = backup_J{cID};
            end
            continue;
        end
        beam_count = nnz(a_new_cell);
        if beam_count > thresholds(1)
            new_potential_active_cells = setdiff(potential_active_cells, adjacentList{cell_idx});
            if beam_count >= thresholds(2)
                new_potential_active_cells = setdiff(new_potential_active_cells, adjacentList2{cell_idx});
            end
            new_potential_active_cells = setdiff(new_potential_active_cells, cell_idx);
        else
            new_potential_active_cells = setdiff(potential_active_cells, cell_idx);
        end
        new_active_cells = [active_cells, cell_idx];
        new_a_var = a_var;
        new_a_var{cell_idx} = a_new_cell;
        for cID = new_active_cells
            numerator_cID = dot(new_a_var{cID}.^2, h_vector_sqr{cID, cID});
            interference_cID = 0;
            for j = new_active_cells
                if j ~= cID
                    a_j_vec = new_a_var{j};
                    h_j_vec = h_vector_sqr{cID, j};
                    if isempty(a_j_vec) || isempty(h_j_vec)
                        continue;
                    end
                    interference_cID = interference_cID + dot(a_j_vec.^2, h_j_vec);
                end
            end
            denominator_cID = interference_cID + sigma2FreqFrag(cID)/Ptx;
            gamma(cID) = numerator_cID / denominator_cID;
            y(cID) = sqrt(w(cID)*(1 + gamma(cID)) * numerator_cID) / ((denominator_cID + numerator_cID) * sqrt(Ptx));
            d{cID} = 2 * y(cID) * sqrt(w(cID)*(1 + gamma(cID))) * h_vector{cID, cID};
            J_prev_cID = zeros(length(epsilon_cell{cID}), 1);
            for j = new_active_cells
                J_prev_cID = J_prev_cID + y(j)^2 * h_vector_sqr{j, cID};
            end
            J{cID} = J_prev_cID + y(cID)^2 * h_vector_sqr{cID, cID};
        end
        sc = new_active_cells(1);
        Rate_comparison = w(sc) * log2(1 + gamma(sc));
        eta_value = Rate_comparison / Demand(sc);
        sub_patterns = recursivePatternGeneration(num_cells, new_active_cells, new_potential_active_cells, new_a_var, adjacentList, adjacentList2, Demand, epsilon_cell, h_vector_sqr, h_vector, w, gamma, y, sigma2FreqFrag, Ptx, ...
            max_depth - 1,thresholds, depth+1, d, J, eta_threshold);
        if isempty(sub_patterns)
            pattern_partial.active_cells = new_active_cells;
            pattern_partial.a_var       = new_a_var;
            pattern_partial.sign_cell   = new_active_cells(1);
            pattern_partial.eta_value   = eta_value;
            patterns{end+1} = pattern_partial;
        else
            patterns = [patterns, sub_patterns];
        end
        for cID = [active_cells, cell_idx]
            gamma(cID) = backup_gamma(cID);
            y(cID)     = backup_y(cID);
            d{cID}     = backup_d{cID};
            J{cID}     = backup_J{cID};
        end
    end
    return;
end
function [a_m_new, feasible, gamma_m_new, y_m_new] = assignBeams(m, a_var, Demand, epsilon_cell, h_vector_sqr, h_vector, w,  y, sigma2FreqFrag, Ptx, active_cells, eta_threshold)
    
    prev_a_var_nonzero = [];
    no_change_count = 0;
    if isempty(active_cells)
        sign_cell = m;
    else 
        sign_cell = active_cells(1);
    end
    max_no_change = 2;
    a_m_new = ones(length(epsilon_cell{m}), 1);
    feasible = true;
    if Demand(m) == 0
        feasible = false;
        return;
    end
    numerator = dot(((a_m_new).^2) , h_vector_sqr{m,m});
    interference = 0;
    for j = active_cells
        if j ~= m
            a_j_vec = a_var{j};
            h_j_vec = h_vector_sqr{m, j};
            if isempty(a_j_vec) || isempty(h_j_vec)
                continue;
            end
            interference = interference + dot((a_j_vec.^2) , h_j_vec);
        end
    end
    denominator = interference + sigma2FreqFrag(m)/Ptx;
    gamma_m_new = numerator / denominator;
    y_m_new = sqrt(w(m)*(1+gamma_m_new)*numerator)/((denominator+numerator)*sqrt(Ptx));
    d_temp = 2 * y_m_new * sqrt(w(m) * (1 + gamma_m_new)) * h_vector{m, m};
    J_prev = zeros(length(epsilon_cell{m}), 1);
    for j = active_cells
        J_prev = J_prev + y(j).^2 * h_vector_sqr{j, m};
    end
    J_temp = J_prev + y_m_new .^2 * h_vector_sqr{m, m};
    history_a_var_prev = cell(max_no_change, 1);
    history_numerator_prev = zeros(max_no_change, 1);
    for i = 1:max_no_change
        history_a_var_prev{i} = zeros(length(epsilon_cell{m}), 1);
        history_numerator_prev(i) = 0;
    end
    bisearch_iter = 0;
    numiter = 0;
    maxnumiter = 2;
    lambda_lower = 0;
    lambda_upper = max(d_temp);
    lambda_prev = lambda_lower;
    R_prev = w(m) * log2(1 + gamma_m_new);
    S_prev = R_prev / Demand(m);
    % Search for a feasible beam assignment for cell m.
    a_var_temp = zeros(length(epsilon_cell{m}), 1);
    reweightedl1_a_var = ones(length(epsilon_cell{m}), 1);
    while true
        bisearch_iter = bisearch_iter + 1;
        lambda = (lambda_upper + lambda_lower) / 2;
        for idx = 1:length(epsilon_cell{m})
            if d_temp(idx) < lambda
                a_var_temp(idx) = 0;
            elseif ((d_temp(idx) - lambda * reweightedl1_a_var(idx)) / J_temp(idx)) <= 1 && d_temp(idx) >= lambda
                a_var_temp(idx) = (d_temp(idx) - lambda * reweightedl1_a_var(idx)) / J_temp(idx);
            elseif ((d_temp(idx) - lambda * reweightedl1_a_var(idx)) / J_temp(idx)) >= 1 && d_temp(idx) >= lambda
                a_var_temp(idx) = 1;
            end
            reweightedl1_a_var(idx)= 1/((a_var_temp(idx))^2+1);
        end
        numerator_new = dot((a_var_temp).^2, h_vector_sqr{m, m});
        gamma_m_new = numerator_new / denominator;
        R_update = w(m) * log2(1 + gamma_m_new);
        S_update = R_update / Demand(m);
        delta_S = abs(R_update - Demand(m)) / Demand(m);
        delta_S_prev_current = abs(R_prev - Demand(m)) / Demand(m);
        alpha = S_prev < 1;
        beta = S_update < 1;
        mu = delta_S <= delta_S_prev_current;
        nu = lambda_prev == lambda_upper;
        if alpha && ~beta && ~nu
            lambda_upper = lambda;
        elseif alpha && ~beta && nu
            lambda_lower = lambda;
        elseif ~alpha && beta && ~nu
            lambda_upper = lambda;
        elseif ~alpha && beta && nu
            lambda_lower = lambda;
        elseif alpha && beta && mu && ~nu
            lambda_lower = lambda;
        elseif alpha && beta && mu && nu
            lambda_upper = lambda;
        elseif alpha && beta && ~mu && ~nu
            lambda_upper = lambda;
        elseif alpha && beta && ~mu && nu
            lambda_lower = lambda;
        elseif ~alpha && ~beta && mu && ~nu
            lambda_lower = lambda;
        elseif ~alpha && ~beta && mu && nu
            lambda_upper = lambda;
        elseif ~alpha && ~beta && ~mu && ~nu
            lambda_upper = lambda;
        elseif ~alpha && ~beta && ~mu && nu
            lambda_lower = lambda;
        end
        lambda_prev = lambda;

        current_nonzero = find(a_var_temp ~= 0);
        if isequal(current_nonzero, prev_a_var_nonzero) && nnz(a_var_temp) == nnz(prev_a_var_nonzero)
            no_change_count = no_change_count + 1;
        else
            no_change_count = 0;
            prev_a_var_nonzero = current_nonzero;
        end
        for i = max_no_change:-1:2
            history_a_var_prev{i} = history_a_var_prev{i-1};
            history_numerator_prev(i) = history_numerator_prev(i-1);
        end
        history_a_var_prev{1} = a_var_temp;
        history_numerator_prev(1) = numerator_new;
        if no_change_count >= max_no_change
            numiter = numiter + 1;
            a_var_temp = history_a_var_prev{max_no_change};
            numerator_new = history_numerator_prev(max_no_change);
            gamma_m_new = numerator_new / denominator;

            if numiter >= maxnumiter
                break;
            end
            y_m_new = sqrt(w(m)*(1 + (numerator_new / denominator)) * numerator_new) / ((denominator + numerator_new) * sqrt(Ptx));
            d_temp = 2 * y_m_new * sqrt(w(m) * (1 + (numerator_new / denominator))) * h_vector{m, m};
            J_temp = J_prev + y_m_new.^2 * h_vector_sqr{m, m};

            lambda_lower = 0;
            lambda_upper = max(d_temp);
            numerator_reset = dot(((a_m_new).^2) , h_vector_sqr{m,m});
            gamma_reset = numerator_reset / denominator;
            R_prev = w(m) * log2(1 + gamma_reset);
            S_prev = R_prev / Demand(m);
            no_change_count = 0;
        else
            R_prev = R_update;
            S_prev = S_update;
        end
        if bisearch_iter > 100
            warning('Reached too many iterations, forced break.');
            break;
        end
    end
    if isempty(active_cells)
        return;
    end
    for n = active_cells
        numerator_n = dot(a_var{n}.^2, h_vector_sqr{n, n});
        interference_n = 0;
        for j = active_cells
            if j ~= n
                a_j_vec = a_var{j};
                h_j_vec = h_vector_sqr{n, j};
                if isempty(a_j_vec) || isempty(h_j_vec)
                    continue;
                end
                interference_n = interference_n + dot(a_j_vec.^2, h_j_vec);
            end
        end
        denominator_n = interference_n + sigma2FreqFrag(n) / Ptx;
        gamma_n = numerator_n / denominator_n;
        Rate_comparison_n = w(n) * log2(1 + gamma_n);

        if n == sign_cell
            required_capacity = eta_threshold(1) * Demand(n);
        else
            required_capacity = eta_threshold(2) * Demand(n);
        end
        if Rate_comparison_n <= required_capacity
            feasible = false;
            return;
        end
    end

    a_m_new = a_var_temp;
end
% Generate and merge routing patterns from all groups.
all_patterns = {};
for c = 1:group_num
    cells_in_cluster = groups{c};
    cells_in_cluster = cells_in_cluster';
    active_cells_initial = [];
    potential_active_cells_initial = cells_in_cluster;
    max_depth = ceil(length(groups{c})/3);
    depth_initial = 1;
    patterns = recursivePatternGeneration(num_cells, active_cells_initial, potential_active_cells_initial, a_var, adjacentList, adjacentList2, Demand, epsilon_cell, h_vector_sqr, h_vector, w, gamma,...
y,sigma2FreqFrag, Ptx, max_depth, thresholds, depth_initial, d, J,eta_threshold);
    all_patterns = [all_patterns, patterns];
end
total_patterns = length(all_patterns);
disp('Recursive routing-pattern generation completed.');
% Match cell requests to satellite beam capacity for each pattern.
Rate_matrix = zeros(num_cells, total_patterns);
Rate_matrix_ideal = zeros(num_cells, total_patterns);
C = cell(total_patterns, 1);
for v = 1:total_patterns
    pattern = all_patterns{v};
    active_cells = pattern.active_cells;
    a_var_pattern = pattern.a_var;
    exchange_matrix = zeros(num_cells, num_sats);
    selected_cells_per_satellite = cell(num_sats, 1);
    for m = active_cells
        a_m = a_var_pattern{m};
        epsilon_m = epsilon_cell{m};
        temp_row = zeros(1, num_sats);
        for idx = 1:length(epsilon_m)
            sat_idx = epsilon_m(idx);
            temp_row(sat_idx) = a_m(idx);
        end
        exchange_matrix(m, :) = temp_row;
    end
    for n = 1:num_sats
        requests = exchange_matrix(:, n);
        requesting_cells = find(requests > 0);
        if ~isempty(requesting_cells)
            a_values = requests(requesting_cells);
            [~, sort_idx] = sort(a_values, 'descend');
            sorted_cells = requesting_cells(sort_idx);
            num_to_connect = min(Q{n}, length(sorted_cells));
            selected_cells = sorted_cells(1:num_to_connect);
            selected_cells_per_satellite{n} = selected_cells;
        else
            selected_cells_per_satellite{n} = [];
        end
    end
    exchange_matrix = zeros(num_cells, num_sats);
    for n = 1:num_sats
        selected_cells = selected_cells_per_satellite{n};
        if ~isempty(selected_cells)
            exchange_matrix(selected_cells, n) = 1;
        end
    end
    for m = active_cells
        epsilon_m = epsilon_cell{m};
        a_new_pattern = zeros(length(epsilon_cell{m}), 1);
        for idx = 1:length(epsilon_m)
            sat_idx = epsilon_m(idx);
            a_value = exchange_matrix(m, sat_idx);
            a_new_pattern(idx) = a_value;
        end
        a_new{m} = a_new_pattern;
    end
    for m = active_cells
        a_m = a_new{m};
        numerator(m) = dot(a_m.^2, h_vector_sqr{m, m});
        interference = 0;
        for j = active_cells
            if j ~= m
                a_j_vec = a_new{j};
                h_j_vec = h_vector_sqr{m, j};
                if isempty(a_j_vec) || isempty(h_j_vec)
                    continue;
                end
                interference = interference + dot(a_j_vec.^2, h_j_vec);
            end
        end
        denominator(m) = interference + sigma2FreqFrag(m) / Ptx;
        gamma_var(m) = numerator(m) / denominator(m);
        Rate_matrix_ideal(m, v) = w(m) * log2(1 + numerator(m) / (sigma2FreqFrag(m) / Ptx));
        Rate_matrix(m, v) = w(m) * log2(1 + gamma_var(m));
    end
    C{v} = exchange_matrix;
end
disp('Matching completed for all routing patterns.');




addpath('YOUR_GUROBI_MATLAB_PATH');
% Select pattern durations through the MILP scheduler.
gurobi_setup
N_period = 100;
model = struct();
model.modelsense = 'max';
num_tao = total_patterns;
obj = [zeros(num_tao, 1); 1];
vtype = [repmat('I', num_tao, 1); 'C'];
lb = [zeros(num_tao, 1); 0];
ub = [repmat(N_period, num_tao, 1); Inf];
A1 = [ones(1, num_tao), 0];
b1 = N_period;
sense1 = '=';
nonzero_demand_idx = find(Demand>0);
RateMatrixNZ = Rate_matrix(nonzero_demand_idx, :);
zeroRowIndicator = all(RateMatrixNZ == 0, 2);
if any(zeroRowIndicator)
    warning([num2str(sum(zeroRowIndicator)), ' positive-demand cells have zero rate and will be excluded.']);
    validIdx = ~zeroRowIndicator;
    RateMatrixNZ = RateMatrixNZ(validIdx, :);
    DemandNZ = Demand(nonzero_demand_idx(validIdx));
else
    DemandNZ = Demand(nonzero_demand_idx);
end
A2 = [RateMatrixNZ, -DemandNZ .* N_period];
b2 = zeros(size(RateMatrixNZ, 1), 1);
sense2 = repmat('>', size(RateMatrixNZ, 1), 1);
model.A = sparse([A1; A2]);
model.rhs = [b1; b2];
model.sense = [sense1; sense2];
model.obj = obj;
model.vtype = vtype;
model.lb = lb;
model.ub = ub;
params = struct();
params.OutputFlag = 0;
params.IterationLimit = 100000;
result = gurobi(model, params);
if strcmp(result.status, 'OPTIMAL') || strcmp(result.status, 'INTEGER_OPTIMAL')
    disp('MILP solved successfully.');
else
    disp(['Solver status: ', result.status]);
end
solution = result.x;
tau = solution(1:num_tao);
eta = solution(end);
disp(['Maximized performance metric eta: ', num2str(eta)]);
solution_struct = struct();
solution_struct.eta = eta;
solution_struct.tau = tau;

% Expand the selected patterns into the final slot schedule.
C_ultimate = cell(N_period, 1);
Rate_matrix_ultimate = zeros(num_cells, N_period);

current_period = 1;

for v = 1:total_patterns
    for cnt = 1:tau(v)
        if current_period > N_period
            warning('The selected pattern count exceeds N_period.');
            break;
        end
        C_ultimate{current_period} = C{v};
        
        Rate_matrix_ultimate(:, current_period) = Rate_matrix(:, v);
        
        current_period = current_period + 1;
    end
end

if current_period - 1 ~= N_period
    warning(['The number of filled slots is ', num2str(current_period-1), ', which does not match N_period = ', num2str(N_period), '.']);
end

disp('C_ultimate and Rate_matrix_ultimate generated.');

% Compute throughput, demand satisfaction, and Jain fairness.
slotLoss = 6/7;
mean_R_all = mean(Rate_matrix_ultimate, 2) * slotLoss;
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
