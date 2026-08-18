%% Generate constellation-comparison data for Figures 9 and 10.
clear;
clc;

demand_value = 3.0;               % demand=3Mbps
cluster_config_mode = 1;          % Heilongjiang area (case1)
nbrOfSetups = 100;                  % Monte Carlo runs

nbrOfUserPerCell = 10;
nbrOfRB = 160;

constellation_configs = [1, 2];

fprintf('[START] figure9to10 | setups=%d | demand=%.2f | usersPerCell=%d | RB=%d | constellations=%s\n', ...
    nbrOfSetups, demand_value, nbrOfUserPerCell, nbrOfRB, mat2str(constellation_configs));

csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
R_over_D_CPDA_cell = cell(1, length(constellation_configs));
R_over_D_CBH_cell = cell(1, length(constellation_configs));
R_over_D_compared_cell  = cell(1, length(constellation_configs));
R_over_D_TDM_cell = cell(1, length(constellation_configs));

Jain_CPDA_cell = cell(1, length(constellation_configs));
Jain_CBH_cell = cell(1, length(constellation_configs));
Jain_compared_cell  = cell(1, length(constellation_configs));
Jain_TDM_cell = cell(1, length(constellation_configs));


for idxConst = 1:length(constellation_configs)
    cc = constellation_configs(idxConst);  % 1 or 2
    fprintf('\n[CONSTELLATION] %d/%d | config=%d\n', idxConst, length(constellation_configs), cc);

    if cc == 1
        csvFile3 = sprintf('../../data/snapshots_Phase1/visibility_relations_cluster%d_config1.csv', cluster_config_mode);
    else
        csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);
    end

    R_over_D_CPDA_total = [];
    R_over_D_CBH_total = [];
    R_over_D_compared_total  = [];
    R_over_D_TDM_total = [];

    Jain_CPDA_total = [];
    Jain_CBH_total = [];
    Jain_compared_total  = [];
    Jain_TDM_total = [];


    for cntSetup = 1:nbrOfSetups
        fprintf('[SETUP] %d/%d | running layoutGenerate...\n', cntSetup, nbrOfSetups);
        tSetup = tic;
        [h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = layoutGenerate(demand_value, csvFile, csvFile2, csvFile3, nbrOfUserPerCell, nbrOfRB);
        fprintf('[SETUP] %d/%d | layout ready | cells=%d | sats=%d\n', cntSetup, nbrOfSetups, num_cells, num_sats);

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
        [C_ultimate_tdm, R_over_D_TDM, Jain_TDM, ~] = Terminal_Dominated_Matching_Scheduling( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
            sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
        fprintf('    [DONE] Terminal-Dominated Matching | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_TDM, Jain_TDM, toc(tAlg));


        R_over_D_CPDA_total = [R_over_D_CPDA_total; R_over_D_CPDA];
        R_over_D_CBH_total = [R_over_D_CBH_total; R_over_D_CBH];
        Jain_CPDA_total     = [Jain_CPDA_total; Jain_CPDA];
        Jain_CBH_total     = [Jain_CBH_total; Jain_CBH];

        R_over_D_compared_total  = [R_over_D_compared_total;  R_over_D_compared];
        Jain_compared_total      = [Jain_compared_total;      Jain_compared];

        R_over_D_TDM_total = [R_over_D_TDM_total; R_over_D_TDM];
        Jain_TDM_total     = [Jain_TDM_total;     Jain_TDM];

        fprintf('[SETUP] %d/%d finished | %.2fs\n', cntSetup, nbrOfSetups, toc(tSetup));
    end

    R_over_D_CPDA_cell{idxConst} = R_over_D_CPDA_total;
    R_over_D_CBH_cell{idxConst} = R_over_D_CBH_total;
    Jain_CPDA_cell{idxConst}     = Jain_CPDA_total;
    Jain_CBH_cell{idxConst}     = Jain_CBH_total;
    R_over_D_compared_cell{idxConst}  = R_over_D_compared_total;
    R_over_D_TDM_cell{idxConst} = R_over_D_TDM_total;

    Jain_compared_cell{idxConst}      = Jain_compared_total;
    Jain_TDM_cell{idxConst}     = Jain_TDM_total;

end
% Reshape setup results into long-form CSV tables.
all_RD = {};
for idxConst = 1:length(constellation_configs)
    cc = constellation_configs(idxConst);
    nCPDA = length(R_over_D_CPDA_cell{idxConst});
    tempCPDA = [repmat({cc}, nCPDA, 1), repmat({'CPDA'}, nCPDA, 1), num2cell(R_over_D_CPDA_cell{idxConst})];

    nTDM = length(R_over_D_TDM_cell{idxConst});
    tempTDM = [repmat({cc}, nTDM, 1), repmat({'Terminal-Dominated Matching'}, nTDM, 1), num2cell(R_over_D_TDM_cell{idxConst})];

    nMultiConn  = length(R_over_D_compared_cell{idxConst});
    tempMultiConn  = [repmat({cc}, nMultiConn, 1),  repmat({'Multi-Conn Heuristic'}, nMultiConn, 1), num2cell(R_over_D_compared_cell{idxConst})];

    nCBH = length(R_over_D_CBH_cell{idxConst});
    tempCBH = [repmat({cc}, nCBH, 1), repmat({'Coordinated BH'}, nCBH, 1), num2cell(R_over_D_CBH_cell{idxConst})];

    all_RD = [all_RD; tempCPDA; tempTDM; tempMultiConn; tempCBH];
end

T_RD = cell2table(all_RD, 'VariableNames', {'config','method','value'});
writetable(T_RD, '../../results/figure9.csv');
fprintf('[CSV] Saved: figure9.csv\n');
all_Jain = {};
for idxConst = 1:length(constellation_configs)
    cc = constellation_configs(idxConst);

    nCPDA = length(Jain_CPDA_cell{idxConst});
    tempCPDA = [repmat({cc}, nCPDA, 1), repmat({'CPDA'}, nCPDA, 1), num2cell(Jain_CPDA_cell{idxConst})];

    nTDM = length(Jain_TDM_cell{idxConst});
    tempTDM = [repmat({cc}, nTDM, 1), repmat({'Terminal-Dominated Matching'}, nTDM, 1), num2cell(Jain_TDM_cell{idxConst})];

    nMultiConn  = length(Jain_compared_cell{idxConst});
    tempMultiConn  = [repmat({cc}, nMultiConn, 1), repmat({'Multi-Conn Heuristic'}, nMultiConn, 1), num2cell(Jain_compared_cell{idxConst})];

    nCBH = length(Jain_CBH_cell{idxConst});
    tempCBH = [repmat({cc}, nCBH, 1), repmat({'Coordinated BH'}, nCBH, 1), num2cell(Jain_CBH_cell{idxConst})];

    all_Jain = [all_Jain; tempCPDA; tempTDM; tempMultiConn; tempCBH];
end

T_Jain = cell2table(all_Jain, 'VariableNames', {'config','method','value'});
writetable(T_Jain, '../../results/figure10.csv');
fprintf('[CSV] Saved: figure10.csv\n');
disp('All simulations completed.');

