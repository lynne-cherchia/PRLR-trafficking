% Lynne Cherchia
% 26 August 2025
% Run downstream kinetic rates tuning (k6, k10, k21)

type = 'weighted';

% Model prediction timepoints %
predTime = [0:60:24*3600];

% Set length of parameter variations
N = 8;

% Define parameter guesses for each rate %
k6Baseline = 0.4;
k6Min = k6Baseline/10;
k6Max = k6Baseline*10;
k6tuningRange = logspace(log10(k6Min),log10(k6Max),N);
k6tuningRange = [k6tuningRange(1:4), k6Baseline, k6tuningRange(5:end)];
% Save k6tuningRange to .mat
save k6tuningRange.mat k6tuningRange

% k6tuningRange_manual = [k6Min k6Baseline/5 k6Baseline/3 k6Baseline/1.5...
%     k6Baseline k6Baseline*1.5 k6Baseline*3 k6Baseline*5 k6Max];

k10Baseline = 0.003;
k10Min = k10Baseline/10;
k10Max = k10Baseline*10;
k10tuningRange = logspace(log10(k10Min),log10(k10Max),N);
k10tuningRange = [k10tuningRange(1:4), k10Baseline, k10tuningRange(5:end)];
save k10tuningRange.mat k10tuningRange

k21Baseline = 0.0163;
k21Min = k21Baseline/10;
k21Max = k21Baseline*10;
k21tuningRange = logspace(log10(k21Min),log10(k21Max),N);
k21tuningRange = [k21tuningRange(1:4), k21Baseline, k21tuningRange(5:end)];
save k21tuningRange.mat k21tuningRange

% Define parameters %
endosomal_params = importdata('traffickingParams_initValues/endosomal_param_sets.mat');

baselineParams = zeros(1,66);
baselineParams(1,1:65) = endosomal_params(1,1:65);

% Import and define initial species values, & 5 RJ initial ratios %
PRLR_initvals = importdata('traffickingParams_initValues/endosomal_initvals.mat');
allMembrane_PRLR_initvals = PRLR_initvals(1,:);
membrane_PRLR_75_25_initvals = PRLR_initvals(2,:);
membrane_PRLR_50_50_initvals = PRLR_initvals(3,:);
membrane_PRLR_25_75_initvals = PRLR_initvals(4,:);
allEndosomal_PRLR_initvals = PRLR_initvals(5,:);

% Define initial conditions %
initialConditions = zeros(5,67);
initialConditions(1,:) = allEndosomal_PRLR_initvals;
initialConditions(2,:) = membrane_PRLR_25_75_initvals;
initialConditions(3,:) = membrane_PRLR_50_50_initvals;
initialConditions(4,:) = membrane_PRLR_75_25_initvals;
initialConditions(5,:) = allMembrane_PRLR_initvals;

numInitialConditions = 5;

% Define storage for outputs %
k6_timeCourses = cell(numInitialConditions,length(k6tuningRange));
k10_timeCourses = cell(numInitialConditions,length(k10tuningRange));
k21_timeCourses = cell(numInitialConditions,length(k21tuningRange));

%numSims = numInitialConditions * length(k6_timeCourses);

% Set up calculation indices from the 67 species %
% k6 %
k6_totCytSTATA = cell(numInitialConditions,length(k6_timeCourses));
k6_totCytSTATB = cell(numInitialConditions,length(k6_timeCourses));
k6_totNuclSTATA = cell(numInitialConditions,length(k6_timeCourses));
k6_totNuclSTATB = cell(numInitialConditions,length(k6_timeCourses));
k6_nuclCytRatioA = cell(numInitialConditions,length(k6_timeCourses));
k6_nuclCytRatioB = cell(numInitialConditions,length(k6_timeCourses));
k6_surface_pJAK = cell(numInitialConditions,length(k6_timeCourses));
k6_internal_pJAK = cell(numInitialConditions,length(k6_timeCourses));
k6_total_pJAK = cell(numInitialConditions,length(k6_timeCourses));
k6_BCL = cell(numInitialConditions,length(k6_timeCourses));

k6_totCytSTATA_norm = cell(numInitialConditions,length(k6_timeCourses));
k6_totCytSTATB_norm = cell(numInitialConditions,length(k6_timeCourses));
k6_totNuclSTATA_norm = cell(numInitialConditions,length(k6_timeCourses));
k6_totNuclSTATB_norm = cell(numInitialConditions,length(k6_timeCourses));
k6_surface_pJAK_norm = cell(numInitialConditions,length(k6_timeCourses));
k6_internal_pJAK_norm = cell(numInitialConditions,length(k6_timeCourses));
k6_total_pJAK_norm = cell(numInitialConditions,length(k6_timeCourses));

% k10 %
k10_totCytSTATA = cell(numInitialConditions,length(k10_timeCourses));
k10_totCytSTATB = cell(numInitialConditions,length(k10_timeCourses));
k10_totNuclSTATA = cell(numInitialConditions,length(k10_timeCourses));
k10_totNuclSTATB = cell(numInitialConditions,length(k10_timeCourses));
k10_nuclCytRatioA = cell(numInitialConditions,length(k10_timeCourses));
k10_nuclCytRatioB = cell(numInitialConditions,length(k10_timeCourses));
k10_surface_pJAK = cell(numInitialConditions,length(k10_timeCourses));
k10_internal_pJAK = cell(numInitialConditions,length(k10_timeCourses));
k10_total_pJAK = cell(numInitialConditions,length(k10_timeCourses));
k10_BCL = cell(numInitialConditions,length(k10_timeCourses));

k10_totCytSTATA_norm = cell(numInitialConditions,length(k10_timeCourses));
k10_totCytSTATB_norm = cell(numInitialConditions,length(k10_timeCourses));
k10_totNuclSTATA_norm = cell(numInitialConditions,length(k10_timeCourses));
k10_totNuclSTATB_norm = cell(numInitialConditions,length(k10_timeCourses));
k10_surface_pJAK_norm = cell(numInitialConditions,length(k10_timeCourses));
k10_internal_pJAK_norm = cell(numInitialConditions,length(k10_timeCourses));
k10_total_pJAK_norm = cell(numInitialConditions,length(k10_timeCourses));

% k21 %
k21_totCytSTATA = cell(numInitialConditions,length(k21_timeCourses));
k21_totCytSTATB = cell(numInitialConditions,length(k21_timeCourses));
k21_totNuclSTATA = cell(numInitialConditions,length(k21_timeCourses));
k21_totNuclSTATB = cell(numInitialConditions,length(k21_timeCourses));
k21_nuclCytRatioA = cell(numInitialConditions,length(k21_timeCourses));
k21_nuclCytRatioB = cell(numInitialConditions,length(k21_timeCourses));
k21_surface_pJAK = cell(numInitialConditions,length(k21_timeCourses));
k21_internal_pJAK = cell(numInitialConditions,length(k21_timeCourses));
k21_total_pJAK = cell(numInitialConditions,length(k21_timeCourses));
k21_BCL = cell(numInitialConditions,length(k21_timeCourses));

k21_totCytSTATA_norm = cell(numInitialConditions,length(k21_timeCourses));
k21_totCytSTATB_norm = cell(numInitialConditions,length(k21_timeCourses));
k21_totNuclSTATA_norm = cell(numInitialConditions,length(k21_timeCourses));
k21_totNuclSTATB_norm = cell(numInitialConditions,length(k21_timeCourses));
k21_surface_pJAK_norm = cell(numInitialConditions,length(k21_timeCourses));
k21_internal_pJAK_norm = cell(numInitialConditions,length(k21_timeCourses));
k21_total_pJAK_norm = cell(numInitialConditions,length(k21_timeCourses));

% Set options %
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'NonNegative',[1:size(initialConditions,2)]);

%% Simulate range of k6 values %%

for i = 1:numInitialConditions
    for j = 1:length(k6tuningRange)

        params = baselineParams;
        params(66) = k6tuningRange(j);
        initvalues = initialConditions(i,:);

        [~, predConc_k6] = ode15s(@core_tune_k6_endosomalCompartment,predTime,initvalues,options,params);

        % Store this struct at the correct index of the cell
        k6_timeCourses{i,j} = predConc_k6;

        data.PRL = predConc_k6(:,1);
        data.RJ = predConc_k6(:,2);
        data.iRJ = predConc_k6(:,3);
        data.S5Ac = predConc_k6(:,4);
        data.S5Bc = predConc_k6(:,5);
        data.SHP2 = predConc_k6(:,6);
        data.PPX = predConc_k6(:,7);
        data.PPN = predConc_k6(:,8);
        data.PRLRJ = predConc_k6(:,9);
        data.iPRLRJ = predConc_k6(:,10);
        data.PRLRJ2 = predConc_k6(:,11);
        data.iPRLRJ2 = predConc_k6(:,12);
        data.PRLRJ2a = predConc_k6(:,13);
        data.iPRLRJ2a = predConc_k6(:,14);
        data.PRLRJ2aS5Ac = predConc_k6(:,15);
        data.iPRLRJ2aS5Ac = predConc_k6(:,16);
        data.PRLRJ2aS5Bc = predConc_k6(:,17);
        data.iPRLRJ2aS5Bc = predConc_k6(:,18);
        data.pS5Ac = predConc_k6(:,19);
        data.pS5Bc = predConc_k6(:,20);
        data.pS5AcpS5Ac = predConc_k6(:,21);
        data.pS5AcpS5Bc = predConc_k6(:,22);
        data.pS5BcpS5Bc = predConc_k6(:,23);
        data.PRLRJ2aSHP2 = predConc_k6(:,24);
        data.iPRLRJ2aSHP2 = predConc_k6(:,25);
        data.PPXpS5Ac = predConc_k6(:,26);
        data.PPXpS5Bc = predConc_k6(:,27);
        data.PPXpS5AcpS5Ac = predConc_k6(:,28);
        data.PPXpS5BcpS5Bc = predConc_k6(:,29);
        data.PPXpS5AcpS5Bc = predConc_k6(:,30);
        data.pS5AcS5Ac = predConc_k6(:,31);
        data.pS5BcS5Bc = predConc_k6(:,32);
        data.S5AcpS5Bc = predConc_k6(:,33);
        data.pS5AcS5Bc = predConc_k6(:,34);
        data.pS5AnpS5An = predConc_k6(:,35);
        data.pS5AnpS5Bn = predConc_k6(:,36);
        data.pS5BnpS5Bn = predConc_k6(:,37);
        data.pS5An = predConc_k6(:,38);
        data.pS5Bn = predConc_k6(:,39);
        data.PPNpS5An = predConc_k6(:,40);
        data.PPNpS5Bn = predConc_k6(:,41);
        data.S5An = predConc_k6(:,42);
        data.S5Bn = predConc_k6(:,43);
        data.PPNpS5AnpS5An = predConc_k6(:,44);
        data.PPNpS5AnpS5Bn = predConc_k6(:,45);
        data.PPNpS5BnpS5Bn = predConc_k6(:,46);
        data.pS5AnS5An = predConc_k6(:,47);
        data.pS5BnS5Bn = predConc_k6(:,48);
        data.S5AnpS5Bn = predConc_k6(:,49);
        data.pS5AnS5Bn = predConc_k6(:,50);
        data.mRNAn = predConc_k6(:,51);
        data.mRNAc = predConc_k6(:,52);
        data.SOCS1 = predConc_k6(:,53);
        data.SOCS1PRLRJ2a = predConc_k6(:,54);
        data.iSOCS1PRLRJ2a = predConc_k6(:,55);
        data.PRLRJ2aS5AcSHP2 = predConc_k6(:,56);
        data.iPRLRJ2aS5AcSHP2 = predConc_k6(:,57);
        data.PRLRJ2aS5BcSHP2 = predConc_k6(:,58);
        data.iPRLRJ2aS5BcSHP2 = predConc_k6(:,59);
        data.SOCS1PRLRJ2aSHP2 = predConc_k6(:,60);
        data.iSOCS1PRLRJ2aSHP2 = predConc_k6(:,61);
        data.mRn = predConc_k6(:,62);
        data.mRc = predConc_k6(:,63);
        data.Rc = predConc_k6(:,64);
        data.mBCLn = predConc_k6(:,65);
        data.mBCLc = predConc_k6(:,66);
        data.BCL = predConc_k6(:,67);

        % Calculate outputs & store %
        k6_totCytSTATA{i,j} = data.S5Ac + data.PRLRJ2aS5Ac + data.iPRLRJ2aS5Ac + data.pS5Ac + 2.*data.pS5AcpS5Ac + data.pS5AcpS5Bc + data.PPXpS5Ac + 2.*data.PPXpS5AcpS5Ac + data.PPXpS5AcpS5Bc + 2.*data.pS5AcS5Ac + data.S5AcpS5Bc + data.pS5AcS5Bc + data.PRLRJ2aS5AcSHP2 + data.iPRLRJ2aS5AcSHP2;
        k6_totCytSTATB{i,j} = data.S5Bc + data.PRLRJ2aS5Bc + data.iPRLRJ2aS5Bc + data.pS5Bc + 2.*data.pS5BcpS5Bc + data.pS5AcpS5Bc + data.PPXpS5Bc + 2.*data.PPXpS5BcpS5Bc + data.PPXpS5AcpS5Bc + 2.*data.pS5BcS5Bc + data.S5AcpS5Bc + data.pS5AcS5Bc + data.PRLRJ2aS5BcSHP2 + data.iPRLRJ2aS5BcSHP2;

        k6_totNuclSTATA{i,j} = 2.*data.pS5AnpS5An + data.pS5AnpS5Bn + data.pS5An + data.PPNpS5An + data.S5An + 2.*data.PPNpS5AnpS5An + data.PPNpS5AnpS5Bn + 2.*data.pS5AnS5An + data.S5AnpS5Bn + data.pS5AnS5Bn;
        k6_totNuclSTATB{i,j} = 2.*data.pS5BnpS5Bn + data.pS5AnpS5Bn + data.pS5Bn + data.PPNpS5Bn + data.S5Bn + 2.*data.PPNpS5BnpS5Bn + data.PPNpS5AnpS5Bn + 2.*data.pS5BnS5Bn + data.S5AnpS5Bn + data.pS5AnS5Bn;

        k6_nuclCytRatioA{i,j} = k6_totNuclSTATA{i,j}./k6_totCytSTATA{i,j};
        k6_nuclCytRatioB{i,j} = k6_totNuclSTATB{i,j}./k6_totCytSTATB{i,j};

        k6_surface_pJAK{i,j} = 2.*data.PRLRJ2a + 2.*data.PRLRJ2aS5Ac + 2.*data.PRLRJ2aS5Bc + 2.*data.PRLRJ2aSHP2 + 2.*data.PRLRJ2aS5AcSHP2 + 2.*data.PRLRJ2aS5BcSHP2 + 2.*data.SOCS1PRLRJ2a + 2.*data.SOCS1PRLRJ2aSHP2;
        k6_surface_pJAK_norm{i,j} = k6_surface_pJAK{i,j}(:,1)./k6_surface_pJAK{i,j}(31,1);
        k6_internal_pJAK{i,j} = 2.*data.iPRLRJ2a + 2.*data.iPRLRJ2aS5Ac + 2.*data.iPRLRJ2aS5Bc + 2.*data.iPRLRJ2aSHP2 + 2.*data.iPRLRJ2aS5AcSHP2 + 2.*data.iPRLRJ2aS5BcSHP2 + 2.*data.iSOCS1PRLRJ2a + 2.*data.iSOCS1PRLRJ2aSHP2;
        k6_internal_pJAK_norm{i,j} = k6_internal_pJAK{i,j}(:,1)./k6_internal_pJAK{i,j}(31,1);

        k6_total_pJAK{i,j} = k6_surface_pJAK{i,j}(:,1) + k6_internal_pJAK{i,j}(:,1);
        k6_total_pJAK_norm{i,j} = k6_total_pJAK{i,j}(:,1)./k6_total_pJAK{i,j}(31,1);

        % Bcl fold change - Bcl(:,ii) = predConc(:,56)./predConc(1,56); -relative to starting amount
        k6_BCL{i,j} = data.BCL./data.BCL(67); 
    end
end

%% Simulate range of k10 values %%

for i = 1:numInitialConditions
    for j = 1:length(k10tuningRange)

        params = baselineParams;
        params(66) = k10tuningRange(j);
        initvalues = initialConditions(i,:);

        [~, predConc_k10] = ode15s(@core_tune_k10_endosomalCompartment,predTime,initvalues,options,params);

        % Store this struct at the correct index of the cell
        k10_timeCourses{i,j} = predConc_k10;

        data.PRL = predConc_k10(:,1);
        data.RJ = predConc_k10(:,2);
        data.iRJ = predConc_k10(:,3);
        data.S5Ac = predConc_k10(:,4);
        data.S5Bc = predConc_k10(:,5);
        data.SHP2 = predConc_k10(:,6);
        data.PPX = predConc_k10(:,7);
        data.PPN = predConc_k10(:,8);
        data.PRLRJ = predConc_k10(:,9);
        data.iPRLRJ = predConc_k10(:,10);
        data.PRLRJ2 = predConc_k10(:,11);
        data.iPRLRJ2 = predConc_k10(:,12);
        data.PRLRJ2a = predConc_k10(:,13);
        data.iPRLRJ2a = predConc_k10(:,14);
        data.PRLRJ2aS5Ac = predConc_k10(:,15);
        data.iPRLRJ2aS5Ac = predConc_k10(:,16);
        data.PRLRJ2aS5Bc = predConc_k10(:,17);
        data.iPRLRJ2aS5Bc = predConc_k10(:,18);
        data.pS5Ac = predConc_k10(:,19);
        data.pS5Bc = predConc_k10(:,20);
        data.pS5AcpS5Ac = predConc_k10(:,21);
        data.pS5AcpS5Bc = predConc_k10(:,22);
        data.pS5BcpS5Bc = predConc_k10(:,23);
        data.PRLRJ2aSHP2 = predConc_k10(:,24);
        data.iPRLRJ2aSHP2 = predConc_k10(:,25);
        data.PPXpS5Ac = predConc_k10(:,26);
        data.PPXpS5Bc = predConc_k10(:,27);
        data.PPXpS5AcpS5Ac = predConc_k10(:,28);
        data.PPXpS5BcpS5Bc = predConc_k10(:,29);
        data.PPXpS5AcpS5Bc = predConc_k10(:,30);
        data.pS5AcS5Ac = predConc_k10(:,31);
        data.pS5BcS5Bc = predConc_k10(:,32);
        data.S5AcpS5Bc = predConc_k10(:,33);
        data.pS5AcS5Bc = predConc_k10(:,34);
        data.pS5AnpS5An = predConc_k10(:,35);
        data.pS5AnpS5Bn = predConc_k10(:,36);
        data.pS5BnpS5Bn = predConc_k10(:,37);
        data.pS5An = predConc_k10(:,38);
        data.pS5Bn = predConc_k10(:,39);
        data.PPNpS5An = predConc_k10(:,40);
        data.PPNpS5Bn = predConc_k10(:,41);
        data.S5An = predConc_k10(:,42);
        data.S5Bn = predConc_k10(:,43);
        data.PPNpS5AnpS5An = predConc_k10(:,44);
        data.PPNpS5AnpS5Bn = predConc_k10(:,45);
        data.PPNpS5BnpS5Bn = predConc_k10(:,46);
        data.pS5AnS5An = predConc_k10(:,47);
        data.pS5BnS5Bn = predConc_k10(:,48);
        data.S5AnpS5Bn = predConc_k10(:,49);
        data.pS5AnS5Bn = predConc_k10(:,50);
        data.mRNAn = predConc_k10(:,51);
        data.mRNAc = predConc_k10(:,52);
        data.SOCS1 = predConc_k10(:,53);
        data.SOCS1PRLRJ2a = predConc_k10(:,54);
        data.iSOCS1PRLRJ2a = predConc_k10(:,55);
        data.PRLRJ2aS5AcSHP2 = predConc_k10(:,56);
        data.iPRLRJ2aS5AcSHP2 = predConc_k10(:,57);
        data.PRLRJ2aS5BcSHP2 = predConc_k10(:,58);
        data.iPRLRJ2aS5BcSHP2 = predConc_k10(:,59);
        data.SOCS1PRLRJ2aSHP2 = predConc_k10(:,60);
        data.iSOCS1PRLRJ2aSHP2 = predConc_k10(:,61);
        data.mRn = predConc_k10(:,62);
        data.mRc = predConc_k10(:,63);
        data.Rc = predConc_k10(:,64);
        data.mBCLn = predConc_k10(:,65);
        data.mBCLc = predConc_k10(:,66);
        data.BCL = predConc_k10(:,67);

        % Calculate outputs & store %
        k10_totCytSTATA{i,j} = data.S5Ac + data.PRLRJ2aS5Ac + data.iPRLRJ2aS5Ac + data.pS5Ac + 2.*data.pS5AcpS5Ac + data.pS5AcpS5Bc + data.PPXpS5Ac + 2.*data.PPXpS5AcpS5Ac + data.PPXpS5AcpS5Bc + 2.*data.pS5AcS5Ac + data.S5AcpS5Bc + data.pS5AcS5Bc + data.PRLRJ2aS5AcSHP2 + data.iPRLRJ2aS5AcSHP2;
        k10_totCytSTATB{i,j} = data.S5Bc + data.PRLRJ2aS5Bc + data.iPRLRJ2aS5Bc + data.pS5Bc + 2.*data.pS5BcpS5Bc + data.pS5AcpS5Bc + data.PPXpS5Bc + 2.*data.PPXpS5BcpS5Bc + data.PPXpS5AcpS5Bc + 2.*data.pS5BcS5Bc + data.S5AcpS5Bc + data.pS5AcS5Bc + data.PRLRJ2aS5BcSHP2 + data.iPRLRJ2aS5BcSHP2;

        k10_totNuclSTATA{i,j} = 2.*data.pS5AnpS5An + data.pS5AnpS5Bn + data.pS5An + data.PPNpS5An + data.S5An + 2.*data.PPNpS5AnpS5An + data.PPNpS5AnpS5Bn + 2.*data.pS5AnS5An + data.S5AnpS5Bn + data.pS5AnS5Bn;
        k10_totNuclSTATB{i,j} = 2.*data.pS5BnpS5Bn + data.pS5AnpS5Bn + data.pS5Bn + data.PPNpS5Bn + data.S5Bn + 2.*data.PPNpS5BnpS5Bn + data.PPNpS5AnpS5Bn + 2.*data.pS5BnS5Bn + data.S5AnpS5Bn + data.pS5AnS5Bn;

        k10_nuclCytRatioA{i,j} = k10_totNuclSTATA{i,j}./k10_totCytSTATA{i,j};
        k10_nuclCytRatioB{i,j} = k10_totNuclSTATB{i,j}./k10_totCytSTATB{i,j};

        k10_surface_pJAK{i,j} = 2.*data.PRLRJ2a + 2.*data.PRLRJ2aS5Ac + 2.*data.PRLRJ2aS5Bc + 2.*data.PRLRJ2aSHP2 + 2.*data.PRLRJ2aS5AcSHP2 + 2.*data.PRLRJ2aS5BcSHP2 + 2.*data.SOCS1PRLRJ2a + 2.*data.SOCS1PRLRJ2aSHP2;
        k10_surface_pJAK_norm{i,j} = k10_surface_pJAK{i,j}(:,1)./k10_surface_pJAK{i,j}(31,1);
        k10_internal_pJAK{i,j} = 2.*data.iPRLRJ2a + 2.*data.iPRLRJ2aS5Ac + 2.*data.iPRLRJ2aS5Bc + 2.*data.iPRLRJ2aSHP2 + 2.*data.iPRLRJ2aS5AcSHP2 + 2.*data.iPRLRJ2aS5BcSHP2 + 2.*data.iSOCS1PRLRJ2a + 2.*data.iSOCS1PRLRJ2aSHP2;
        k10_internal_pJAK_norm{i,j} = k10_internal_pJAK{i,j}(:,1)./k10_internal_pJAK{i,j}(31,1);

        k10_total_pJAK{i,j} = k10_surface_pJAK{i,j}(:,1) + k10_internal_pJAK{i,j}(:,1);
        k10_total_pJAK_norm{i,j} = k10_total_pJAK{i,j}(:,1)./k10_total_pJAK{i,j}(31,1);

        % Bcl fold change - Bcl(:,ii) = predConc(:,56)./predConc(1,56); -relative to starting amount
        k10_BCL{i,j} = data.BCL./data.BCL(67);
    end
end

%% Simulate range of k21 values %%

for i = 1:numInitialConditions
    for j = 1:length(k21tuningRange)

        params = baselineParams;
        params(66) = k21tuningRange(j);
        initvalues = initialConditions(i,:);

        [~, predConc_k21] = ode15s(@core_tune_k21_endosomalCompartment,predTime,initvalues,options,params);

        % Store this struct at the correct index of the cell
        k21_timeCourses{i,j} = predConc_k21;

        data.PRL = predConc_k21(:,1);
        data.RJ = predConc_k21(:,2);
        data.iRJ = predConc_k21(:,3);
        data.S5Ac = predConc_k21(:,4);
        data.S5Bc = predConc_k21(:,5);
        data.SHP2 = predConc_k21(:,6);
        data.PPX = predConc_k21(:,7);
        data.PPN = predConc_k21(:,8);
        data.PRLRJ = predConc_k21(:,9);
        data.iPRLRJ = predConc_k21(:,10);
        data.PRLRJ2 = predConc_k21(:,11);
        data.iPRLRJ2 = predConc_k21(:,12);
        data.PRLRJ2a = predConc_k21(:,13);
        data.iPRLRJ2a = predConc_k21(:,14);
        data.PRLRJ2aS5Ac = predConc_k21(:,15);
        data.iPRLRJ2aS5Ac = predConc_k21(:,16);
        data.PRLRJ2aS5Bc = predConc_k21(:,17);
        data.iPRLRJ2aS5Bc = predConc_k21(:,18);
        data.pS5Ac = predConc_k21(:,19);
        data.pS5Bc = predConc_k21(:,20);
        data.pS5AcpS5Ac = predConc_k21(:,21);
        data.pS5AcpS5Bc = predConc_k21(:,22);
        data.pS5BcpS5Bc = predConc_k21(:,23);
        data.PRLRJ2aSHP2 = predConc_k21(:,24);
        data.iPRLRJ2aSHP2 = predConc_k21(:,25);
        data.PPXpS5Ac = predConc_k21(:,26);
        data.PPXpS5Bc = predConc_k21(:,27);
        data.PPXpS5AcpS5Ac = predConc_k21(:,28);
        data.PPXpS5BcpS5Bc = predConc_k21(:,29);
        data.PPXpS5AcpS5Bc = predConc_k21(:,30);
        data.pS5AcS5Ac = predConc_k21(:,31);
        data.pS5BcS5Bc = predConc_k21(:,32);
        data.S5AcpS5Bc = predConc_k21(:,33);
        data.pS5AcS5Bc = predConc_k21(:,34);
        data.pS5AnpS5An = predConc_k21(:,35);
        data.pS5AnpS5Bn = predConc_k21(:,36);
        data.pS5BnpS5Bn = predConc_k21(:,37);
        data.pS5An = predConc_k21(:,38);
        data.pS5Bn = predConc_k21(:,39);
        data.PPNpS5An = predConc_k21(:,40);
        data.PPNpS5Bn = predConc_k21(:,41);
        data.S5An = predConc_k21(:,42);
        data.S5Bn = predConc_k21(:,43);
        data.PPNpS5AnpS5An = predConc_k21(:,44);
        data.PPNpS5AnpS5Bn = predConc_k21(:,45);
        data.PPNpS5BnpS5Bn = predConc_k21(:,46);
        data.pS5AnS5An = predConc_k21(:,47);
        data.pS5BnS5Bn = predConc_k21(:,48);
        data.S5AnpS5Bn = predConc_k21(:,49);
        data.pS5AnS5Bn = predConc_k21(:,50);
        data.mRNAn = predConc_k21(:,51);
        data.mRNAc = predConc_k21(:,52);
        data.SOCS1 = predConc_k21(:,53);
        data.SOCS1PRLRJ2a = predConc_k21(:,54);
        data.iSOCS1PRLRJ2a = predConc_k21(:,55);
        data.PRLRJ2aS5AcSHP2 = predConc_k21(:,56);
        data.iPRLRJ2aS5AcSHP2 = predConc_k21(:,57);
        data.PRLRJ2aS5BcSHP2 = predConc_k21(:,58);
        data.iPRLRJ2aS5BcSHP2 = predConc_k21(:,59);
        data.SOCS1PRLRJ2aSHP2 = predConc_k21(:,60);
        data.iSOCS1PRLRJ2aSHP2 = predConc_k21(:,61);
        data.mRn = predConc_k21(:,62);
        data.mRc = predConc_k21(:,63);
        data.Rc = predConc_k21(:,64);
        data.mBCLn = predConc_k21(:,65);
        data.mBCLc = predConc_k21(:,66);
        data.BCL = predConc_k21(:,67);

        % Calculate outputs & store %
        k21_totCytSTATA{i,j} = data.S5Ac + data.PRLRJ2aS5Ac + data.iPRLRJ2aS5Ac + data.pS5Ac + 2.*data.pS5AcpS5Ac + data.pS5AcpS5Bc + data.PPXpS5Ac + 2.*data.PPXpS5AcpS5Ac + data.PPXpS5AcpS5Bc + 2.*data.pS5AcS5Ac + data.S5AcpS5Bc + data.pS5AcS5Bc + data.PRLRJ2aS5AcSHP2 + data.iPRLRJ2aS5AcSHP2;
        k21_totCytSTATB{i,j} = data.S5Bc + data.PRLRJ2aS5Bc + data.iPRLRJ2aS5Bc + data.pS5Bc + 2.*data.pS5BcpS5Bc + data.pS5AcpS5Bc + data.PPXpS5Bc + 2.*data.PPXpS5BcpS5Bc + data.PPXpS5AcpS5Bc + 2.*data.pS5BcS5Bc + data.S5AcpS5Bc + data.pS5AcS5Bc + data.PRLRJ2aS5BcSHP2 + data.iPRLRJ2aS5BcSHP2;

        k21_totNuclSTATA{i,j} = 2.*data.pS5AnpS5An + data.pS5AnpS5Bn + data.pS5An + data.PPNpS5An + data.S5An + 2.*data.PPNpS5AnpS5An + data.PPNpS5AnpS5Bn + 2.*data.pS5AnS5An + data.S5AnpS5Bn + data.pS5AnS5Bn;
        k21_totNuclSTATB{i,j} = 2.*data.pS5BnpS5Bn + data.pS5AnpS5Bn + data.pS5Bn + data.PPNpS5Bn + data.S5Bn + 2.*data.PPNpS5BnpS5Bn + data.PPNpS5AnpS5Bn + 2.*data.pS5BnS5Bn + data.S5AnpS5Bn + data.pS5AnS5Bn;

        k21_nuclCytRatioA{i,j} = k21_totNuclSTATA{i,j}./k21_totCytSTATA{i,j};
        k21_nuclCytRatioB{i,j} = k21_totNuclSTATB{i,j}./k21_totCytSTATB{i,j};

        k21_surface_pJAK{i,j} = 2.*data.PRLRJ2a + 2.*data.PRLRJ2aS5Ac + 2.*data.PRLRJ2aS5Bc + 2.*data.PRLRJ2aSHP2 + 2.*data.PRLRJ2aS5AcSHP2 + 2.*data.PRLRJ2aS5BcSHP2 + 2.*data.SOCS1PRLRJ2a + 2.*data.SOCS1PRLRJ2aSHP2;
        k21_surface_pJAK_norm{i,j} = k21_surface_pJAK{i,j}(:,1)./k21_surface_pJAK{i,j}(31,1);
        k21_internal_pJAK{i,j} = 2.*data.iPRLRJ2a + 2.*data.iPRLRJ2aS5Ac + 2.*data.iPRLRJ2aS5Bc + 2.*data.iPRLRJ2aSHP2 + 2.*data.iPRLRJ2aS5AcSHP2 + 2.*data.iPRLRJ2aS5BcSHP2 + 2.*data.iSOCS1PRLRJ2a + 2.*data.iSOCS1PRLRJ2aSHP2;
        k21_internal_pJAK_norm{i,j} = k21_internal_pJAK{i,j}(:,1)./k21_internal_pJAK{i,j}(31,1);

        k21_total_pJAK{i,j} = k21_surface_pJAK{i,j}(:,1) + k21_internal_pJAK{i,j}(:,1);
        k21_total_pJAK_norm{i,j} = k21_total_pJAK{i,j}(:,1)./k21_total_pJAK{i,j}(31,1);

        % Bcl fold change - Bcl(:,ii) = predConc(:,56)./predConc(1,56); -relative to starting amount
        k21_BCL{i,j} = data.BCL./data.BCL(67);
    end
end

%% Import simulation outputs %%

k6_nuclCytRatioA = importdata('downstreamRateTuning/k6_nuclCytRatioA.mat');
k6_nuclCytRatioB = importdata('downstreamRateTuning/k6_nuclCytRatioB.mat');
k6_total_pJAK_norm = importdata('downstreamRateTuning/k6_total_pJAK_norm.mat');
k6_BCL = importdata('downstreamRateTuning/k6_BCL.mat');

k10_nuclCytRatioA = importdata('downstreamRateTuning/k10_nuclCytRatioA.mat');
k10_nuclCytRatioB = importdata('downstreamRateTuning/k10_nuclCytRatioB.mat');
k10_total_pJAK_norm = importdata('downstreamRateTuning/k10_total_pJAK_norm.mat');
k10_BCL = importdata('downstreamRateTuning/k10_BCL.mat');

k21_nuclCytRatioA = importdata('downstreamRateTuning/k21_nuclCytRatioA.mat');
k21_nuclCytRatioB = importdata('downstreamRateTuning/k21_nuclCytRatioB.mat');
k21_total_pJAK_norm = importdata('downstreamRateTuning/k21_total_pJAK_norm.mat');
k21_BCL = importdata('downstreamRateTuning/k21_BCL.mat');

%%

% Import base model outputs & experimental data %
baseModel_outputs = importdata('baseModelOutputs/predConc_baseModel.mat');
baseModel_time = importdata('baseModelOutputs/predTime_baseModel.mat');
baseModel_StatA_nucleus = importdata('baseModelOutputs/Stat_nucleusA_baseModel.mat');
baseModel_StatA_cyto = importdata('baseModelOutputs/Stat_cytoA_baseModel.mat');
baseModel_pStatA_norm = importdata('baseModelOutputs/pStatA_norm_baseModel.mat');
baseModel_StatB_nucleus = importdata('baseModelOutputs/Stat_nucleusB_baseModel.mat');
baseModel_StatB_cyto = importdata('baseModelOutputs/Stat_cytoB_baseModel.mat');
baseModel_pStatB_norm = importdata('baseModelOutputs/pStatB_norm_baseModel.mat');
baseModel_BCL = importdata('baseModelOutputs_t24h/Bcl_24h_baseModel.mat');
baseModel_rec_total_allMembrane = importdata('baseModelOutputs/rec_total_baseModel.mat');
baseModel_pJAK_norm = importdata('baseModelOutputs/pJAK_norm_baseModel.mat');
baseModel_nucleus_cyto_ratioA = baseModel_StatA_nucleus./baseModel_StatA_cyto;
baseModel_nucleus_cyto_ratioB = baseModel_StatB_nucleus./baseModel_StatB_cyto;

% Extended time course experimental data were shared privately %

% Plot titles for saving files % 
plotTitles = cell(1,numInitialConditions);
plotTitles{1,1} = 'Initial condition 1: RJ_0 all internal';
plotTitles{1,2} = 'Initial condition 2: RJ_0 internal 75:25 surface';
plotTitles{1,3} = 'Initial condition 3: RJ_0 internal 50:50 surface';
plotTitles{1,4} = 'Initial condition 4: RJ_0 internal 25:75 surface';
plotTitles{1,5} = 'Initial condition 5: RJ_0 all surface';

fileExtension = '.pdf';

%% Nucleus:Cytosol STAT5A with k6 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k6_nuclCytRatioA{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]); % darkest shade of red
    plot(predTime/3600, k6_nuclCytRatioA{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k6_nuclCytRatioA{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k6_nuclCytRatioA{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k6_nuclCytRatioA{k,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k6_nuclCytRatioA{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k6_nuclCytRatioA{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k6_nuclCytRatioA{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k6_nuclCytRatioA{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]); % darkest shade of blue
    
    plot(baseModel_time/3600, baseModel_nucleus_cyto_ratioA, '--', 'LineWidth', 2, 'color','k');
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'fontname','arial',...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Ratio of nuclear/cytosolic STAT5A');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% Nucleus:Cytosol STAT5B with k6 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k6_nuclCytRatioB{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k6_nuclCytRatioB{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k6_nuclCytRatioB{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k6_nuclCytRatioB{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k6_nuclCytRatioB{1,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k6_nuclCytRatioB{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k6_nuclCytRatioB{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k6_nuclCytRatioB{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k6_nuclCytRatioB{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_nucleus_cyto_ratioB, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_translocb(1,:), exp_mu_translocb(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 7])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Ratio of nuclear/cytosolic STAT5B');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% Relative pJAK with k6 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k6_total_pJAK_norm{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k6_total_pJAK_norm{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k6_total_pJAK_norm{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k6_total_pJAK_norm{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k6_total_pJAK_norm{k,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k6_total_pJAK_norm{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k6_total_pJAK_norm{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k6_total_pJAK_norm{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k6_total_pJAK_norm{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_pJAK_norm, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_pjak(1,:), exp_mu_pjak(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Relative pJAK');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% BCL fold change from t0 with k6 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k6_BCL{1,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k6_BCL{1,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k6_BCL{1,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k6_BCL{1,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k6_BCL{1,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k6_BCL{1,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k6_BCL{1,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k6_BCL{1,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k6_BCL{1,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(predTime/3600, baseModel_BCL, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_bcl(1,:), exp_mu_bcl(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 24 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('BCL-xL fold change');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end

%% Nucleus:Cytosol STAT5A with k10 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k10_nuclCytRatioA{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k10_nuclCytRatioA{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k10_nuclCytRatioA{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k10_nuclCytRatioA{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k10_nuclCytRatioA{k,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k10_nuclCytRatioA{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k10_nuclCytRatioA{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k10_nuclCytRatioA{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k10_nuclCytRatioA{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_nucleus_cyto_ratioA, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_transloca(1,:), exp_mu_transloca(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Ratio of nuclear/cytosolic STAT5A');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% Nucleus:Cytosol STAT5B with k10 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k10_nuclCytRatioB{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k10_nuclCytRatioB{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k10_nuclCytRatioB{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k10_nuclCytRatioB{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k10_nuclCytRatioB{1,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k10_nuclCytRatioB{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k10_nuclCytRatioB{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k10_nuclCytRatioB{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k10_nuclCytRatioB{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_nucleus_cyto_ratioB, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_translocb(1,:), exp_mu_translocb(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 7])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Ratio of nuclear/cytosolic STAT5B');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% Relative pJAK with k10 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k10_total_pJAK_norm{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k10_total_pJAK_norm{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k10_total_pJAK_norm{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k10_total_pJAK_norm{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k10_total_pJAK_norm{k,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k10_total_pJAK_norm{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k10_total_pJAK_norm{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k10_total_pJAK_norm{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k10_total_pJAK_norm{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_pJAK_norm, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_pjak(1,:), exp_mu_pjak(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Relative pJAK');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% BCL fold change from t0 with k10 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k10_BCL{1,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k10_BCL{1,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k10_BCL{1,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k10_BCL{1,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k10_BCL{1,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k10_BCL{1,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k10_BCL{1,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k10_BCL{1,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k10_BCL{1,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(predTime/3600, baseModel_BCL, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_bcl(1,:), exp_mu_bcl(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 24 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('BCL-xL fold change');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end

%% Nucleus:Cytosol STAT5A with k21 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k21_nuclCytRatioA{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k21_nuclCytRatioA{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k21_nuclCytRatioA{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k21_nuclCytRatioA{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k21_nuclCytRatioA{k,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k21_nuclCytRatioA{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k21_nuclCytRatioA{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k21_nuclCytRatioA{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k21_nuclCytRatioA{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_nucleus_cyto_ratioA, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_transloca(1,:), exp_mu_transloca(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 3])
    % xticks([0.5 1.0 1.5 3.0 6.0])
    % xtickangle(45)
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Ratio of nuclear/cytosolic STAT5A');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% Nucleus:Cytosol STAT5B with k21 tuning %% 

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k21_nuclCytRatioB{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k21_nuclCytRatioB{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k21_nuclCytRatioB{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k21_nuclCytRatioB{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k21_nuclCytRatioB{1,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k21_nuclCytRatioB{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k21_nuclCytRatioB{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k21_nuclCytRatioB{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k21_nuclCytRatioB{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_nucleus_cyto_ratioB, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_translocb(1,:), exp_mu_translocb(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 7])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Ratio of nuclear/cytosolic STAT5B');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% Relative pJAK with k21 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k21_total_pJAK_norm{k,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k21_total_pJAK_norm{k,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k21_total_pJAK_norm{k,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k21_total_pJAK_norm{k,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k21_total_pJAK_norm{k,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k21_total_pJAK_norm{k,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k21_total_pJAK_norm{k,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k21_total_pJAK_norm{k,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k21_total_pJAK_norm{k,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(baseModel_time/3600, baseModel_pJAK_norm, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_pjak(1,:), exp_mu_pjak(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 6 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('Relative pJAK');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end
%% BCL fold change from t0 with k21 tuning %%

filePath = '/...';

for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]);
    hold on

    % below baseline %
    plot(predTime/3600, k21_BCL{1,1}, 'LineWidth', 2, 'color', [0.6471 0 0.1490]);
    plot(predTime/3600, k21_BCL{1,2}, 'LineWidth', 2, 'color', [0.8431 0.1882 0.1529]);
    plot(predTime/3600, k21_BCL{1,3}, 'LineWidth', 2, 'color', [0.9569 0.4275 0.2627]);
    plot(predTime/3600, k21_BCL{1,4}, 'LineWidth', 2, 'color', [0.9922 0.6824 0.3804]);
    
    plot(predTime/3600, k21_BCL{1,5}, 'LineWidth', 2, 'color', 'k');
    % above baseline %
    plot(predTime/3600, k21_BCL{1,6}, 'LineWidth', 2, 'color', [0.6706 0.8510 0.9137]);
    plot(predTime/3600, k21_BCL{1,7}, 'LineWidth', 2, 'color', [0.4549 0.6784 0.8196]);
    plot(predTime/3600, k21_BCL{1,8}, 'LineWidth', 2, 'color', [0.2706 0.4588 0.7059]);
    plot(predTime/3600, k21_BCL{1,9}, 'LineWidth', 2, 'color', [0.1922 0.2118 0.5843]);
    
    plot(predTime/3600, baseModel_BCL, '--', 'LineWidth', 2, 'color','k');
    scatter(exp_mu_bcl(1,:), exp_mu_bcl(2,:),100,'MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
    
    fontSize = 18;
    set(gca,'YDir','normal',...
        'fontsize',fontSize,...
        'linewidth',1,...
        'layer','top',...
        'clipping','on',...
        'box','on',...
        'xscale','linear');
    set(gca,'TickLength',[0.025, 0.01])
    axis([0 24 0 3])
    title([plotTitles(k)])
    xlabel('Time (h)');
    ylabel('BCL-xL fold change');
    fileName = char(strcat(plotTitles(k),fileExtension));
    exportgraphics(gcf,fullfile(filePath, fileName),'ContentType','vector');
    hold off
end