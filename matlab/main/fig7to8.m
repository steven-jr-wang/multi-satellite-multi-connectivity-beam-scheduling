%% Generate demand-sweep data for Figures 7 and 8.
clear;
clc;

constellation_config_Mode = 2;   % Fixed: Config2
cluster_config_mode = 1;           % Heilongjiang internal area (case1)
nbrOfSetups = 100;                 % Number of Monte Carlo simulations
nbrOfUserPerCell = 10;
nbrOfRB = 160;
demand_values = [1.5,3.0,6.0];
colors = lines(length(demand_values));
csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);

R_over_D_CPDA_cell = cell(1, length(demand_values));
R_over_D_CBH_cell = cell(1, length(demand_values));
Jain_CPDA_cell     = cell(1, length(demand_values));
Jain_CBH_cell     = cell(1, length(demand_values));
R_over_D_FM_cell = cell(1, length(demand_values));
Jain_FM_cell     = cell(1, length(demand_values));
R_over_D_TDM_cell = cell(1, length(demand_values));
Jain_TDM_cell     = cell(1, length(demand_values));
fprintf('[START] figure7to8 | setups=%d | usersPerCell=%d | RB=%d | demands=%s\n', ...
    nbrOfSetups, nbrOfUserPerCell, nbrOfRB, mat2str(demand_values));

for idx = 1:length(demand_values)
    current_demand = demand_values(idx);
    fprintf('\n[DEMAND] idx=%d/%d | demand=%.2f Mbps\n', idx, length(demand_values), current_demand);

    R_over_D_CPDA_total = [];
    R_over_D_CBH_total = [];
    Jain_CPDA_total = [];
    Jain_CBH_total = [];

    R_over_D_FM_total = [];
    Jain_FM_total     = [];
    R_over_D_TDM_total = [];
    Jain_TDM_total     = [];



    
    for cntOfSetup = 1:nbrOfSetups
        fprintf('[SETUP] %d/%d | running layoutGenerate...\n', cntOfSetup, nbrOfSetups);
tSetup = tic;

        [h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, sigma2,sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = ...
            layoutGenerate(current_demand, csvFile, csvFile2, csvFile3, nbrOfUserPerCell,nbrOfRB);

        fprintf('[SETUP] %d/%d | layout ready | cells=%d | sats=%d\n', cntOfSetup, nbrOfSetups, num_cells, num_sats);
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

        R_over_D_FM_total = [R_over_D_FM_total; R_over_D_compared];
        Jain_FM_total     = [Jain_FM_total;     Jain_compared];
        R_over_D_TDM_total = [R_over_D_TDM_total; R_over_D_TDM];
        Jain_TDM_total     = [Jain_TDM_total;     Jain_TDM];

        fprintf('[SETUP] %d/%d finished | %.2fs\n', cntOfSetup, nbrOfSetups, toc(tSetup));


   


    end

    R_over_D_CPDA_cell{idx} = R_over_D_CPDA_total;
    R_over_D_CBH_cell{idx} = R_over_D_CBH_total;
    Jain_CPDA_cell{idx}     = Jain_CPDA_total;
    Jain_CBH_cell{idx}     = Jain_CBH_total;

    R_over_D_FM_cell{idx} = R_over_D_FM_total;
    Jain_FM_cell{idx}     = Jain_FM_total;
    R_over_D_TDM_cell{idx} = R_over_D_TDM_total;
    Jain_TDM_cell{idx}     = Jain_TDM_total;


end

all_RD = {};
% Reshape setup results into long-form CSV tables.
for idx = 1:length(demand_values)
    current_demand = demand_values(idx);
    nCPDA = length(R_over_D_CPDA_cell{idx});
    tempCPDA = [num2cell(repmat(current_demand, nCPDA, 1)), repmat({'CPDA'}, nCPDA, 1), num2cell(R_over_D_CPDA_cell{idx})];
    nCBH = length(R_over_D_CBH_cell{idx});
    tempCBH = [num2cell(repmat(current_demand, nCBH, 1)), repmat({'Coordinated BH'}, nCBH, 1), num2cell(R_over_D_CBH_cell{idx})];


    nFM = length(R_over_D_FM_cell{idx});
    tempFM = [num2cell(repmat(current_demand, nFM, 1)), repmat({'Multi-Conn Heuristic'}, nFM, 1), num2cell(R_over_D_FM_cell{idx})];


    nTDM = length(R_over_D_TDM_cell{idx});
    tempTDM = [num2cell(repmat(current_demand, nTDM, 1)), repmat({'Terminal-Dominated Matching'}, nTDM, 1), ...
           num2cell(R_over_D_TDM_cell{idx})];

    all_RD = [all_RD; tempCPDA; tempTDM; tempFM; tempCBH];







end
T_RD = cell2table(all_RD, 'VariableNames', {'demand','method','value'});
writetable(T_RD, '../../results/figure7.csv');
fprintf('[CSV] Saved: figure7.csv\n');

all_Jain = {};
for idx = 1:length(demand_values)
    current_demand = demand_values(idx);
    nCPDA = length(Jain_CPDA_cell{idx});
    tempCPDA = [num2cell(repmat(current_demand, nCPDA, 1)), repmat({'CPDA'}, nCPDA, 1), num2cell(Jain_CPDA_cell{idx})];
    nCBH = length(Jain_CBH_cell{idx});
    tempCBH = [num2cell(repmat(current_demand, nCBH, 1)), repmat({'Coordinated BH'}, nCBH, 1), num2cell(Jain_CBH_cell{idx})];

    nFM = length(Jain_FM_cell{idx});
    tempFM = [num2cell(repmat(current_demand, nFM, 1)), repmat({'Multi-Conn Heuristic'}, nFM, 1), num2cell(Jain_FM_cell{idx})];
    nTDM = length(Jain_TDM_cell{idx});
    tempTDM = [num2cell(repmat(current_demand, nTDM, 1)), repmat({'Terminal-Dominated Matching'}, nTDM, 1), ...
           num2cell(Jain_TDM_cell{idx})];

    all_Jain = [all_Jain; tempCPDA; tempTDM; tempFM; tempCBH];




end
T_Jain = cell2table(all_Jain, 'VariableNames', {'demand','method','value'});
writetable(T_Jain, '../../results/figure8.csv');
fprintf('[CSV] Saved: figure8.csv\n');
disp('All simulations completed.');
