%% Generate per-cell resource data for Figures 13 and 14.
clear;
clc;

constellation_config_Mode = 2;   % Fixed: Config2
cluster_config_mode = 1;         % Heilongjiang internal area (case1)
nbrOfSetups = 100;                 % Number of Monte Carlo simulations
nbrOfUserPerCell = 10;
nbrOfRB = 160;

demand_value = 3.0;

csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);

R_over_D_CPDA_cell = cell(1, length(demand_value));
R_over_D_CBH_cell = cell(1, length(demand_value));
Jain_CPDA_cell     = cell(1, length(demand_value));
Jain_CBH_cell     = cell(1, length(demand_value));
R_over_D_FM_cell  = cell(1, length(demand_value));
Jain_FM_cell      = cell(1, length(demand_value));
R_over_D_TDM_cell = cell(1, length(demand_value));
Jain_TDM_cell     = cell(1, length(demand_value));

fprintf('[START] figure13to14 | setups=%d | demand=%.2f | usersPerCell=%d | RB=%d\n', ...
    nbrOfSetups, demand_value, nbrOfUserPerCell, nbrOfRB);

% Accumulate per-cell metrics over Monte Carlo setups.
for cntOfSetup = 1:nbrOfSetups
    fprintf('[SETUP] %d/%d | running layoutGenerate...\n', cntOfSetup, nbrOfSetups);
    tSetup = tic;

    [h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = ...
        layoutGenerate(demand_value, csvFile, csvFile2, csvFile3, nbrOfUserPerCell, nbrOfRB);

    if cntOfSetup == 1
        Demand_setup_sum = zeros(num_cells, 1);

        R_CPDA_sum        = zeros(num_cells, 1);
        R_MultiConn_sum  = zeros(num_cells, 1);
        R_CBH_sum        = zeros(num_cells, 1);
        R_TDM_sum        = zeros(num_cells, 1);

        Unmet_CPDA_sum       = zeros(num_cells, 1);
        Unmet_MultiConn_sum = zeros(num_cells, 1);
        Unmet_CBH_sum       = zeros(num_cells, 1);
        Unmet_TDM_sum       = zeros(num_cells, 1);

        Unused_CPDA_sum       = zeros(num_cells, 1);
        Unused_MultiConn_sum = zeros(num_cells, 1);
        Unused_CBH_sum       = zeros(num_cells, 1);
        Unused_TDM_sum       = zeros(num_cells, 1);
    end

    fprintf('[SETUP] %d/%d | layout ready | cells=%d | sats=%d\n', cntOfSetup, nbrOfSetups, num_cells, num_sats);

    fprintf('  - [ALG] CPDA ...\n');
    tAlg = tic;

    [C_ultimate, R_over_D_CPDA, Jain_CPDA, mean_R_CPDA] = proposedMethod( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);

    fprintf('    [DONE] CPDA | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_CPDA, Jain_CPDA, toc(tAlg));

    fprintf('  - [ALG] Multi-Conn Heuristic ...\n');
    tAlg = tic;

    [C_ultimate_compared, R_over_D_compared, Jain_compared, mean_R_compared] = Multi_Connection_Heuristic_Scheduling( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
        sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);

    fprintf('    [DONE] Multi-Conn Heuristic | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_compared, Jain_compared, toc(tAlg));

    fprintf('  - [ALG] Coordinated BH ...\n');
    tAlg = tic;

    [~, C, R_over_D_CBH, Jain_CBH, mean_R_CBH] = CoordinatedBeamHopping( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, Ptx, w, Demand, wgs84, nbrOfRB);

    fprintf('    [DONE] Coordinated BH | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_CBH, Jain_CBH, toc(tAlg));

    fprintf('  - [ALG] Terminal-Dominated Matching ...\n');
    tAlg = tic;

    [C_ultimate_tdm, R_over_D_TDM, Jain_TDM, mean_R_TDM] = Terminal_Dominated_Matching_Scheduling( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
        sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);

    fprintf('    [DONE] Terminal-Dominated Matching | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_TDM, Jain_TDM, toc(tAlg));

    Demand_setup_sum = Demand_setup_sum + Demand;

    R_CPDA_sum       = R_CPDA_sum       + mean_R_CPDA;
    R_MultiConn_sum = R_MultiConn_sum + mean_R_compared;
    R_CBH_sum       = R_CBH_sum       + mean_R_CBH;
    R_TDM_sum       = R_TDM_sum       + mean_R_TDM;

    Unmet_CPDA_sum       = Unmet_CPDA_sum       + max(Demand - mean_R_CPDA, 0);
    Unmet_MultiConn_sum = Unmet_MultiConn_sum + max(Demand - mean_R_compared, 0);
    Unmet_CBH_sum       = Unmet_CBH_sum       + max(Demand - mean_R_CBH, 0);
    Unmet_TDM_sum       = Unmet_TDM_sum       + max(Demand - mean_R_TDM, 0);

    Unused_CPDA_sum       = Unused_CPDA_sum       + max(mean_R_CPDA - Demand, 0);
    Unused_MultiConn_sum = Unused_MultiConn_sum + max(mean_R_compared - Demand, 0);
    Unused_CBH_sum       = Unused_CBH_sum       + max(mean_R_CBH - Demand, 0);
    Unused_TDM_sum       = Unused_TDM_sum       + max(mean_R_TDM - Demand, 0);
    fprintf('[SETUP] %d/%d finished | %.2fs\n', cntOfSetup, nbrOfSetups, toc(tSetup));
end

D_bar = Demand_setup_sum / nbrOfSetups;

Rbar_CPDA       = R_CPDA_sum / nbrOfSetups;
Rbar_TDM       = R_TDM_sum / nbrOfSetups;
Rbar_MultiConn = R_MultiConn_sum / nbrOfSetups;
Rbar_CBH       = R_CBH_sum / nbrOfSetups;

Unmetbar_CPDA       = Unmet_CPDA_sum / nbrOfSetups;
Unmetbar_TDM       = Unmet_TDM_sum / nbrOfSetups;
Unmetbar_MultiConn = Unmet_MultiConn_sum / nbrOfSetups;
Unmetbar_CBH       = Unmet_CBH_sum / nbrOfSetups;

Unusedbar_CPDA       = Unused_CPDA_sum / nbrOfSetups;
Unusedbar_TDM       = Unused_TDM_sum / nbrOfSetups;
Unusedbar_MultiConn = Unused_MultiConn_sum / nbrOfSetups;
Unusedbar_CBH       = Unused_CBH_sum / nbrOfSetups;

cell_id = (1:num_cells).';

T_cell_rate = table(cell_id, D_bar, Rbar_CPDA, Rbar_TDM, Rbar_MultiConn, Rbar_CBH, ...
    'VariableNames', {'cell_id','D_bar','Rbar_CPDA','Rbar_TDM','Rbar_MultiConn','Rbar_CBH'});
writetable(T_cell_rate, '../../results/figure_cell_rate.csv');

T_cell_unmet = table(cell_id, Unmetbar_CPDA, Unmetbar_TDM, Unmetbar_MultiConn, Unmetbar_CBH, ...
    'VariableNames', {'cell_id','Unmetbar_CPDA','Unmetbar_TDM','Unmetbar_MultiConn','Unmetbar_CBH'});
writetable(T_cell_unmet, '../../results/figure_cell_unmet.csv');

T_cell_unused = table(cell_id, Unusedbar_CPDA, Unusedbar_TDM, Unusedbar_MultiConn, Unusedbar_CBH, ...
    'VariableNames', {'cell_id','Unusedbar_CPDA','Unusedbar_TDM','Unusedbar_MultiConn','Unusedbar_CBH'});
writetable(T_cell_unused, '../../results/figure_cell_unused.csv');

method = {'CPDA'; 'Terminal-Dominated Matching'; 'Multi-Conn Heuristic'; 'Coordinated BH'};

total_provided_sum = [sum(R_CPDA_sum); sum(R_TDM_sum); sum(R_MultiConn_sum); sum(R_CBH_sum)];

total_unmet_sum = [sum(Unmet_CPDA_sum); sum(Unmet_TDM_sum); sum(Unmet_MultiConn_sum); sum(Unmet_CBH_sum)];

total_unused_sum = [sum(Unused_CPDA_sum); sum(Unused_TDM_sum); sum(Unused_MultiConn_sum); sum(Unused_CBH_sum)];

T_summary = table(method, total_provided_sum, total_unmet_sum, total_unused_sum, ...
    'VariableNames', {'method','total_provided_sum_Mbps','total_unmet_sum_Mbps','total_unused_sum_Mbps'});

writetable(T_summary, '../../results/figure_cell_summary.csv');

fprintf('[CSV] Saved per-cell rate, unmet, unused, and summary data.\n');
disp('All simulations completed.');
