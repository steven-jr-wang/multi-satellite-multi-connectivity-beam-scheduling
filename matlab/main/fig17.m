%% Measure layout and scheduling runtimes for constellation 2.
clear;
clc;

demand_value = 3.0;              % Fixed demand: 3 Mbps
constellation_config_Mode = 2;   % Fixed constellation configuration: 2
cluster_config_mode = 1;
nbrOfSetups = 100;
nbrOfUserPerCell = 10;
nbrOfRB = 160;

csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);

runtime_rows = cell(nbrOfSetups * 5, 3);
row_idx = 1;

fprintf('[START] figure17 | setups=%d | demand=%.2f | constellation=%d | usersPerCell=%d | RB=%d\n', ...
    nbrOfSetups, demand_value, constellation_config_Mode, nbrOfUserPerCell, nbrOfRB);

for cntOfSetup = 1:nbrOfSetups
    fprintf('[SETUP] %d/%d | running layoutGenerate...\n', cntOfSetup, nbrOfSetups);
    tSetup = tic;
    tLayout = tic;

    [h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = ...
        layoutGenerate(demand_value, csvFile, csvFile2, csvFile3, nbrOfUserPerCell, nbrOfRB);

    layout_runtime_s = round(toc(tLayout), 3);
    runtime_rows(row_idx, :) = {cntOfSetup, 'Layout Generation', layout_runtime_s};
    row_idx = row_idx + 1;
    fprintf('[SETUP] %d/%d | layout ready | cells=%d | sats=%d | runtime=%.3fs\n', ...
        cntOfSetup, nbrOfSetups, num_cells, num_sats, layout_runtime_s);

    fprintf('  - [ALG] CPDA ...\n');
    tAlg = tic;
    [~, ~, ~, ~] = proposedMethod( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
    runtime_s = round(toc(tAlg), 3);
    runtime_rows(row_idx, :) = {cntOfSetup, 'CPDA', runtime_s};
    row_idx = row_idx + 1;
    fprintf('    [DONE] CPDA | runtime=%.3fs\n', runtime_s);

    fprintf('  - [ALG] Multi-Conn Heuristic ...\n');
    tAlg = tic;
    [~, ~, ~, ~] = Multi_Connection_Heuristic_Scheduling( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
        sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
    runtime_s = round(toc(tAlg), 3);
    runtime_rows(row_idx, :) = {cntOfSetup, 'Multi-Conn Heuristic', runtime_s};
    row_idx = row_idx + 1;
    fprintf('    [DONE] Multi-Conn Heuristic | runtime=%.3fs\n', runtime_s);

    fprintf('  - [ALG] Coordinated BH ...\n');
    tAlg = tic;
    [~, ~, ~, ~, ~] = CoordinatedBeamHopping( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, Ptx, w, Demand, wgs84, nbrOfRB);
    runtime_s = round(toc(tAlg), 3);
    runtime_rows(row_idx, :) = {cntOfSetup, 'Coordinated BH', runtime_s};
    row_idx = row_idx + 1;
    fprintf('    [DONE] Coordinated BH | runtime=%.3fs\n', runtime_s);

    fprintf('  - [ALG] Terminal-Dominated Matching ...\n');
    tAlg = tic;
    [~, ~, ~, ~] = Terminal_Dominated_Matching_Scheduling( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
        sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
    runtime_s = round(toc(tAlg), 3);
    runtime_rows(row_idx, :) = {cntOfSetup, 'Terminal-Dominated Matching', runtime_s};
    row_idx = row_idx + 1;
    fprintf('    [DONE] Terminal-Dominated Matching | runtime=%.3fs\n', runtime_s);

    fprintf('[SETUP] %d/%d finished | %.3fs\n', cntOfSetup, nbrOfSetups, toc(tSetup));
end

T_runtime = cell2table(runtime_rows, ...
    'VariableNames', {'setup_id', 'method', 'runtime_s'});
out_csv = '../../results/figure17_runtime.csv';
writetable(T_runtime, out_csv);

fprintf('[CSV] Saved: %s\n', out_csv);
disp('All simulations completed.');
