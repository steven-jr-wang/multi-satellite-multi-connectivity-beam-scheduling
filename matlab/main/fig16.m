%% Generate bandwidth-dependent ASE data for Figure 16.
clear; clc;

cluster_config_mode = 1;
nbrOfSetups = 100;           % 
nbrOfUserPerCell = 10;      % Fixed user count per cell

RB_list = [25, 160, 270];
B_MHz_list = [5, 30, 50];            % Used for ASE normalization
demand_list = [0.5, 3.0, 5.0];       % Demand scales with bandwidth

A_cell_km2 = 1770.347654491;

csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);

fprintf('[START] ASE-BWGap | setups=%d | usersPerCell=%d\n', nbrOfSetups, nbrOfUserPerCell);
all_rows = {}; % columns: scenario_id, B_MHz, nbrOfRB, demand_value, ASE_frag, ASE_nofrag

for k = 1:length(RB_list)
    nbrOfRB = RB_list(k);
    B_MHz = B_MHz_list(k);
    demand_value = demand_list(k);

    fprintf('\n[BW] RB=%d | B=%.1fMHz | demand=%.2f\n', nbrOfRB, B_MHz, demand_value);

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
        [~, ~, ~, mean_R_frag] = proposedMethod( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, [], sigma2FreqFrag, Ptx, w, Demand, wgs84);
        ASE_frag = sum(mean_R_frag) / (B_MHz * A_tot_km2);
        fprintf('    [DONE] CPDA | ASE=%.6e | %.2fs\n', ASE_frag, toc(tAlg));

        fprintf('  - [ALG] CPDA (Without Fragmentation) ...\n');
        tAlg = tic;
        [~, ~, ~, mean_R_nofrag] = proposedMethodwoFragmentation( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, sigma2, Ptx, [], Demand, wgs84, nbrOfRB);
        ASE_nofrag = sum(mean_R_nofrag) / (B_MHz * A_tot_km2);
        fprintf('    [DONE] CPDA (Without Fragmentation) | ASE=%.6e | %.2fs\n', ASE_nofrag, toc(tAlg));

        

        all_rows = [all_rows;
    {cntOfSetup, B_MHz, nbrOfRB, demand_value, ASE_frag, ASE_nofrag}
];


        fprintf('[SETUP] %d/%d finished | %.2fs\n', cntOfSetup, nbrOfSetups, toc(tSetup));
    end
end

T = cell2table(all_rows, 'VariableNames', ...
    {'scenario_id','B_MHz','nbrOfRB','demand_value','ASE_frag','ASE_nofrag'});

out_csv = '../../results/figure_ASE_bandwidth_gap.csv';
writetable(T, out_csv);

fprintf('[CSV] Saved: %s\n', out_csv);
disp('All simulations completed.');
