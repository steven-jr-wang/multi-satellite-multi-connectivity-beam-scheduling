%% Generate ASE comparison data for Figure 15.
clear; clc;

cluster_config_mode = 1;
nbrOfSetups = 100;                % 
nbrOfUserPerCell = 10;
nbrOfRB = 160;                   % 30 MHz
B_MHz = 30;                      % Used directly for ASE normalization
demand_value = 3.0;              % Fixed demand: 3 Mbps

A_cell_km2 = 1770.347654491;

csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);

fprintf('[START] ASE-ThreeMethods | setups=%d | usersPerCell=%d | RB=%d | B=%.1fMHz | demand=%.1f\n', ...
    nbrOfSetups, nbrOfUserPerCell, nbrOfRB, B_MHz, demand_value);

all_ASE = {};  % cell array for table rows

for cntOfSetup = 1:nbrOfSetups
    fprintf('[SETUP] %d/%d | running layoutGenerate...\n', cntOfSetup, nbrOfSetups);
    tSetup = tic;

    [h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = ...
        layoutGenerate(demand_value, csvFile, csvFile2, csvFile3, nbrOfUserPerCell, nbrOfRB);

    fprintf('[SETUP] %d/%d | layout ready | cells=%d | sats=%d\n', cntOfSetup, nbrOfSetups, num_cells, num_sats);
    A_tot_km2 = num_cells * A_cell_km2;

    fprintf('  - [ALG] CPDA ...\n');
    tAlg = tic;
    [~, ~, ~, mean_R_CPDA] = proposedMethod( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, [], sigma2FreqFrag, Ptx, w, Demand, wgs84);
    ASE_CPDA = sum(mean_R_CPDA) / (B_MHz * A_tot_km2);
    fprintf('    [DONE] CPDA | ASE=%.6e | %.2fs\n', ASE_CPDA, toc(tAlg));

    fprintf('  - [ALG] Coordinated BH ...\n');
    tAlg = tic;
    [~, ~, ~, ~, mean_R_CBH] = CoordinatedBeamHopping( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, Ptx, [], Demand, wgs84, nbrOfRB);
    ASE_CBH = sum(mean_R_CBH) / (B_MHz * A_tot_km2);
    fprintf('    [DONE] Coordinated BH | ASE=%.6e | %.2fs\n', ASE_CBH, toc(tAlg));

    fprintf('  - [ALG] CPDA (Without Fragmentation) ...\n');
    tAlg = tic;
    [~, ~, ~, mean_R_NoFrag] = proposedMethodwoFragmentation( ...
        h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
        num_cells, num_sats, sigma2, Ptx, [], Demand, wgs84, nbrOfRB);
    ASE_NoFrag = sum(mean_R_NoFrag) / (B_MHz * A_tot_km2);
    fprintf('    [DONE] CPDA (Without Fragmentation) | ASE=%.6e | %.2fs\n', ASE_NoFrag, toc(tAlg));

    all_ASE = [all_ASE;
        {cntOfSetup, 'CPDA', ASE_CPDA};
        {cntOfSetup, 'Coordinated BH', ASE_CBH};
        {cntOfSetup, 'CPDA (Without Fragmentation)', ASE_NoFrag};
    ];

    fprintf('[SETUP] %d/%d finished | %.2fs\n', cntOfSetup, nbrOfSetups, toc(tSetup));
end

T_ASE = cell2table(all_ASE, 'VariableNames', {'scenario_id','method','ASE'});
T_ASE.B_MHz = repmat(B_MHz, height(T_ASE), 1);
T_ASE.demand_value = repmat(demand_value, height(T_ASE), 1);
T_ASE.num_cells = repmat(num_cells, height(T_ASE), 1);
T_ASE.A_cell_km2 = repmat(A_cell_km2, height(T_ASE), 1);

out_csv = '../../results/figure_ASE_three_methods.csv';
writetable(T_ASE, out_csv);

fprintf('[CSV] Saved: %s\n', out_csv);
disp('All simulations completed.');
