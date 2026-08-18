%% Generate the simulation data used by Figures 1-6.
clear;
clc;

demand_value = 3.0;              % Fixed demand: 3 Mbps
constellation_config_Mode = 2;   % Fixed constellation configuration: 2
nbrOfUserPerCell = 10;
nbrOfRB = 160;

cluster_config_mode = 1;
csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase1/visibility_relations_cluster%d_config1.csv', cluster_config_mode);

fprintf('[START] figure1to4 | setups=1 | demand=%.2f | usersPerCell=%d | RB=%d\n', ...
    demand_value, nbrOfUserPerCell, nbrOfRB);
fprintf('[SETUP] 1/1 | running layoutGenerate...\n');
tSetup = tic;
[h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
    num_cells, num_sats, sigma2,sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = ...
    layoutGenerate(demand_value, csvFile, csvFile2, csvFile3, nbrOfUserPerCell,nbrOfRB);
fprintf('[SETUP] 1/1 | layout ready | cells=%d | sats=%d\n', num_cells, num_sats);

fprintf('  - [ALG] CPDA ...\n');
tAlg = tic;
[C_ultimate, R_over_D_CPDA, Jain_CPDA, ~] = proposedMethod(...
    h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
    num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
fprintf('    [DONE] CPDA | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_CPDA, Jain_CPDA, toc(tAlg));

fprintf('  - [ALG] Multi-Conn Heuristic ...\n');
tAlg = tic;
[C_ultimate_compared, R_over_D_compared, Jain_compared, ~] = Multi_Connection_Heuristic_Scheduling( ...
    h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
    sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
fprintf('    [DONE] Multi-Conn Heuristic | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_compared, Jain_compared, toc(tAlg));

fprintf('  - [ALG] Coordinated BH ...\n');
tAlg = tic;
[~, C, R_over_D_CBH, Jain_CBH, ~] = CoordinatedBeamHopping(...
    h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
    num_cells, num_sats, sigma2, Ptx, w, Demand, wgs84, nbrOfRB);
fprintf('    [DONE] Coordinated BH | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_CBH, Jain_CBH, toc(tAlg));

fprintf('  - [ALG] Terminal-Dominated Matching ...\n');
tAlg = tic;
[C_ultimate_TDM, R_over_D_TDM, Jain_TDM, ~] = Terminal_Dominated_Matching_Scheduling( ...
    h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
    sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
fprintf('    [DONE] Terminal-Dominated Matching | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_TDM, Jain_TDM, toc(tAlg));
fprintf('[SETUP] 1/1 finished | %.2fs\n', toc(tSetup));




% Convert slot-level associations into connection-count matrices.
N_period=100;
connCount_C = zeros(num_cells, N_period);
connCount_C_ultimate = zeros(num_cells, N_period);
connCount_C_ultimate_TDM = zeros(num_cells, N_period);
connCount_C_ultimate_MC = zeros(num_cells, N_period);
for t = 1:N_period
    connCount_C(:, t) = sum(C{t}, 2);
    connCount_C_ultimate(:, t) = sum(C_ultimate{t}, 2);
    connCount_C_ultimate_MC(:, t) = sum(C_ultimate_compared{t}, 2);
    connCount_C_ultimate_TDM(:, t) = sum(C_ultimate_TDM{t}, 2);

end



max_conn = max([connCount_C_ultimate(:); connCount_C_ultimate_MC(:); connCount_C_ultimate_TDM(:)]);


if max_conn == 0, max_conn = 1; end

img_C_ultimate = connCount_C_ultimate / max_conn;
img_C = connCount_C / max_conn;
img_C_ultimate_MC = connCount_C_ultimate_MC / max_conn;






% Export the connection matrices and their dimensions.
writematrix(connCount_C, '../../results/figure3_connCount_C.csv');
disp('[CSV] Saved: figure3_connCount_C.csv');

T_info = table(N_period, num_cells, 'VariableNames', {'N_period', 'num_cells'});
writetable(T_info, '../../results/figure3_info.csv');
disp('[CSV] Saved: figure3_info.csv');




writematrix(connCount_C_ultimate, '../../results/figure4_connCount_C_ultimate.csv');
disp('[CSV] Saved: figure4_connCount_C_ultimate.csv');

T_info = table(N_period, num_cells, 'VariableNames', {'N_period', 'num_cells'});
writetable(T_info, '../../results/figure4_info.csv');
disp('[CSV] Saved: figure4_info.csv');


writematrix(connCount_C_ultimate_MC, '../../results/figure5_connCount_C_ultimate_MC.csv');
disp('[CSV] Saved: figure5_connCount_C_ultimate_MC.csv');


T_info5 = table(N_period, num_cells, 'VariableNames', {'N_period', 'num_cells'});
writetable(T_info5, '../../results/figure5_info.csv');
disp('[CSV] Saved: figure5_info.csv');

writematrix(connCount_C_ultimate_TDM, '../../results/figure6_connCount_C_ultimate_TDM.csv');
disp('[CSV] Saved: figure6_connCount_C_ultimate_TDM.csv');

T_info6 = table(N_period, num_cells, 'VariableNames', {'N_period', 'num_cells'});
writetable(T_info6, '../../results/figure6_info.csv');
disp('[CSV] Saved: figure6_info.csv');
disp('All simulations completed.');
