%% Generate traffic-distribution data for Figures 11 and 12.
clear;
clc;

demand_value = 3.0;                 % Fixed demand = 3 Mbps
cluster_config_mode = 1;            % Heilongjiang internal area (case1)
nbrOfSetups = 100;                    % Number of Monte Carlo simulations
nbrOfUserPerCell = 10;
nbrOfRB = 160;

csvFile  = '../../data/cell_clusters/output_cell_cluster1.csv';
csvFile2 = '../../data/relative_density.csv';
csvFile3 = sprintf('../../data/snapshots_Phase2/visibility_relations_cluster%d_config2.csv', cluster_config_mode);

scenarioNames = {'Uniform', 'Moderate', 'Extreme'};
numScenarios = length(scenarioNames);

R_over_D_CPDA_cell = cell(1, numScenarios);
R_over_D_CBH_cell = cell(1, numScenarios);
R_over_D_compared_cell = cell(1, numScenarios);
R_over_D_TDM_cell = cell(1, numScenarios);

Jain_CPDA_cell = cell(1, numScenarios);
Jain_CBH_cell = cell(1, numScenarios);
Jain_compared_cell = cell(1, numScenarios);
Jain_TDM_cell = cell(1, numScenarios);


relDensityParams.useGenerated = true;

R_over_D_CPDA_total = [];
R_over_D_CBH_total = [];
R_over_D_compared_total = [];
R_over_D_TDM_total = [];

Jain_CPDA_total = [];
Jain_CBH_total = [];
Jain_compared_total = [];
Jain_TDM_total = [];

fprintf('[START] figure11to12 | demand=%.2f Mbps | scenarios=%d | setups=%d | usersPerCell=%d | RB=%d\n', ...
    demand_value, numScenarios, nbrOfSetups, nbrOfUserPerCell, nbrOfRB);
tAll = tic;
for sc = 1:numScenarios
    fprintf('\n[SCENARIO] %d/%d | %s\n', sc, numScenarios, scenarioNames{sc});
    tScenario = tic;
    switch scenarioNames{sc}
        case 'Uniform'
            relDensityParams.distType = 0;  % Uniform distribution
        case 'Moderate'
            relDensityParams.distType = 2;  % Lognormal distribution
            relDensityParams.mu = 0;
            relDensityParams.sigma = 0.8;   % Moderate nonuniformity
        case 'Extreme'
            relDensityParams.distType = 2;  % Lognormal distribution
            relDensityParams.mu = 0;
            relDensityParams.sigma = 2.0;   % Extreme nonuniformity
    end

    fprintf('[SCENARIO] relDensityParams: useGenerated=%d | distType=%d\n', relDensityParams.useGenerated, relDensityParams.distType);
    for cntOfSetup = 1:nbrOfSetups
        fprintf('[SETUP] %d/%d | running layoutGenerate...\n', cntOfSetup, nbrOfSetups);
        tSetup = tic;

        [h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, sigma2,sigma2FreqFrag, Ptx, w, Demand, wgs84, satellite_positions, norm_var_users] = ...
            layoutGenerate(demand_value, csvFile, csvFile2, csvFile3, nbrOfUserPerCell,nbrOfRB, relDensityParams);
        fprintf('[SETUP] %d/%d | layout ready | cells=%d | sats=%d\n', cntOfSetup, nbrOfSetups, num_cells, num_sats);

        fprintf('  - [ALG] CPDA ...\n');
        tAlg = tic;
        [~, R_over_D_CPDA, Jain_CPDA, ~] = proposedMethod( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
        fprintf('    [DONE] CPDA | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_CPDA, Jain_CPDA, toc(tAlg));

        fprintf('  - [ALG] Coordinated BH ...\n');
        tAlg = tic;
        [~, ~, R_over_D_CBH, Jain_CBH, ~] = CoordinatedBeamHopping( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, ...
            num_cells, num_sats, sigma2, Ptx, w, Demand, wgs84, nbrOfRB);
        fprintf('    [DONE] Coordinated BH | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_CBH, Jain_CBH, toc(tAlg));

        fprintf('  - [ALG] Multi-Conn Heuristic ...\n');
        tAlg = tic;
        [~, R_over_D_compared, Jain_compared, ~] = Multi_Connection_Heuristic_Scheduling( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
            sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
        fprintf('    [DONE] Multi-Conn Heuristic | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_compared, Jain_compared, toc(tAlg));

        fprintf('  - [ALG] Terminal-Dominated Matching ...\n');
        tAlg = tic;
        [~, R_over_D_TDM, Jain_TDM, ~] = Terminal_Dominated_Matching_Scheduling( ...
            h_vector, h_vector_sqr, epsilon_cell, data, center_data, num_cells, num_sats, ...
            sigma2, sigma2FreqFrag, Ptx, w, Demand, wgs84);
        fprintf('    [DONE] Terminal-Dominated Matching | R/D=%.4f | Jain=%.4f | %.2fs\n', R_over_D_TDM, Jain_TDM, toc(tAlg));


        R_over_D_CPDA_total = [R_over_D_CPDA_total; R_over_D_CPDA];
        R_over_D_CBH_total = [R_over_D_CBH_total; R_over_D_CBH];
        Jain_CPDA_total     = [Jain_CPDA_total; Jain_CPDA];
        Jain_CBH_total     = [Jain_CBH_total; Jain_CBH];
        R_over_D_compared_total = [R_over_D_compared_total; R_over_D_compared];
        Jain_compared_total     = [Jain_compared_total;     Jain_compared];

        R_over_D_TDM_total = [R_over_D_TDM_total; R_over_D_TDM];
        Jain_TDM_total     = [Jain_TDM_total;     Jain_TDM];
        fprintf('[SETUP] %d/%d finished | %.2fs\n', cntOfSetup, nbrOfSetups, toc(tSetup));

    end
    R_over_D_CPDA_cell{sc} = R_over_D_CPDA_total;
    R_over_D_CBH_cell{sc} = R_over_D_CBH_total;
    Jain_CPDA_cell{sc}     = Jain_CPDA_total;
    Jain_CBH_cell{sc}     = Jain_CBH_total;
    R_over_D_compared_cell{sc} = R_over_D_compared_total;
    R_over_D_TDM_cell{sc}      = R_over_D_TDM_total;

    Jain_compared_cell{sc}     = Jain_compared_total;
    Jain_TDM_cell{sc}          = Jain_TDM_total;

end


all_RD = {};
% Reshape setup results into long-form CSV tables.
for sc = 1:numScenarios
    scenarioName = scenarioNames{sc};

    nCPDA = length(R_over_D_CPDA_cell{sc});
    tempCPDA = [repmat({scenarioName}, nCPDA, 1), repmat({'CPDA'}, nCPDA, 1), num2cell(R_over_D_CPDA_cell{sc})];

    nTDM = length(R_over_D_TDM_cell{sc});
    tempTDM = [repmat({scenarioName}, nTDM, 1), repmat({'Terminal-Dominated Matching'}, nTDM, 1), num2cell(R_over_D_TDM_cell{sc})];

    nMC = length(R_over_D_compared_cell{sc});
    tempMC = [repmat({scenarioName}, nMC, 1), repmat({'Multi-Conn Heuristic'}, nMC, 1), num2cell(R_over_D_compared_cell{sc})];

    nCBH = length(R_over_D_CBH_cell{sc});
    tempCBH = [repmat({scenarioName}, nCBH, 1), repmat({'Coordinated BH'}, nCBH, 1), num2cell(R_over_D_CBH_cell{sc})];

    all_RD = [all_RD; tempCPDA; tempTDM; tempMC; tempCBH];
end
T_RD = cell2table(all_RD, 'VariableNames', {'scenario','method','value'});
writetable(T_RD, '../../results/figure11.csv');
fprintf('[CSV] Saved: figure11.csv\n');
all_Jain = {};
for sc = 1:numScenarios
    scenarioName = scenarioNames{sc};

    nCPDA = length(Jain_CPDA_cell{sc});
    tempCPDA = [repmat({scenarioName}, nCPDA, 1), repmat({'CPDA'}, nCPDA, 1), num2cell(Jain_CPDA_cell{sc})];

    nTDM = length(Jain_TDM_cell{sc});
    tempTDM = [repmat({scenarioName}, nTDM, 1), repmat({'Terminal-Dominated Matching'}, nTDM, 1), num2cell(Jain_TDM_cell{sc})];

    nMC = length(Jain_compared_cell{sc});
    tempMC = [repmat({scenarioName}, nMC, 1), repmat({'Multi-Conn Heuristic'}, nMC, 1), num2cell(Jain_compared_cell{sc})];

    nCBH = length(Jain_CBH_cell{sc});
    tempCBH = [repmat({scenarioName}, nCBH, 1), repmat({'Coordinated BH'}, nCBH, 1), num2cell(Jain_CBH_cell{sc})];

    all_Jain = [all_Jain; tempCPDA; tempTDM; tempMC; tempCBH];
end
T_Jain = cell2table(all_Jain, 'VariableNames', {'scenario','method','value'});
writetable(T_Jain, '../../results/figure12.csv');
fprintf('[CSV] Saved: figure12.csv\n');
disp('All simulations completed.');
