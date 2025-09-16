% Lynne Cherchia
% 26 August 2025
% Run meshgridSims & save heatmap plots

% FOR USE ON HPC (high-performance computing). 
% Run internalization model using meshgrid rate values & save output plots only.

type = 'weighted';

% Model prediction timepoints %
predTime = [0:60:24*3600];

% Set length of parameter variations
N = 250;

% Define meshgrid for each rate %
kintfreeBaseline = 0.000774;
kintfreeMin = kintfreeBaseline/10;
kintfreeMax = kintfreeBaseline*10;
kintfreeMesh = linspace(kintfreeMin,kintfreeMax,N);

krecfreeBaseline = 0.029;
krecfreeMin = krecfreeBaseline/10;
krecfreeMax = krecfreeBaseline*10;
krecfreeMesh = linspace(krecfreeMin,krecfreeMax,N);

kintboundBaseline = 0.0205;
kintboundMin = kintboundBaseline/10;
kintboundMax = kintboundBaseline*10;
kintboundMesh = linspace(kintboundMin,kintboundMax,N);

krecboundBaseline = 0.0142;
krecboundMin = krecboundBaseline/10;
krecboundMax = krecboundBaseline;
krecboundMesh = linspace(krecboundMin,krecboundMax,N);

% Define parameters %
endosomal_params = [...
0.00623687525965315	5.60000000000000e-05	0.00560000000000000	0.443315639726432	0.200000000000000	0.0286921194929444	0.00254591929512555	0.800000000000000	0.400000000000000	0.000256721177985165	6.51757781357932	0.00284100355864639	0.100000000000000	0.165244239510054	0.100000000000000	0.00143646308044359	0.100000000000000	0.000754106930384025	0.200000000000000	0.00300000000000000	0.00395806205704936	0.200000000000000	0.00300000000000000	2.00000000000000e-07	0.687901957753156	0.121168221998327	0.313027298092336	0.0489952673652474	0.00100001869260733	0.200000000000000	0.0152239690045907	0.0241408974486121	0.0402348290810201	0.00629008702628322	0.00925012797982826	0.0100000000000000	400	0.000486863909742820	0.0100000000000000	0.0162962614462024	0.100000000000000	0.000231015457716723	0.000132157361500284	0.0338205094268733	0.0128240121564962	400	0.00100000000000000	9.62202117475019e-06	0.0112685379318887	0.0100000000000000	0.500000000000000	1.20000000000000	1.36000000000000	96.3215089063843	0.00177973790137018	400	0.00100000000000000	0.000500000000000000	0.0100000000000000	9.58200113926253e-05	0.00133001957198602	0.000774000000000000	0.0290000000000000	0.0205000000000000	0.0142000000000000	0;...
0.00623687500000000	5.60000000000000e-05	0.00560000000000000	0.443315640000000	0.200000000000000	0.0286921190000000	0.00254591900000000	0.800000000000000	0.400000000000000	0.000256721000000000	6.51757781400000	0.00284100400000000	0.100000000000000	0.165244240000000	0.100000000000000	0.00143646300000000	0.100000000000000	0.000754107000000000	0.200000000000000	0.00300000000000000	0.00395806200000000	0.200000000000000	0.00300000000000000	2.00000000000000e-07	0.687901958000000	0.121168222000000	0.313027298000000	0.0489952670000000	0.00100001900000000	0.200000000000000	0.0152239690000000	0.0241408970000000	0.0402348290000000	0.00629008700000000	0.00925012800000000	0.0100000000000000	400	0.000486864000000000	0.0100000000000000	0.0162962610000000	0.100000000000000	0.000231015000000000	0.000132157000000000	0.0338205090000000	0.0128240120000000	400	0.00100000000000000	9.62000000000000e-06	0.0112685380000000	0.0100000000000000	0.500000000000000	1.20000000000000	1.36000000000000	96.3215089100000	0.00177973800000000	400	0.00100000000000000	0.000500000000000000	0.0100000000000000	9.58000000000000e-05	0.00133002000000000	0.00155000000000000	0.0781000000000000	0.0830000000000000	0.000777000000000000	0;...
0.00623687500000000	5.60000000000000e-05	0.00560000000000000	0.443315640000000	0.200000000000000	0.0286921190000000	0.00254591900000000	0.800000000000000	0.400000000000000	0.000256721000000000	6.51757781400000	0.00284100400000000	0.100000000000000	0.165244240000000	0.100000000000000	0.00143646300000000	0.100000000000000	0.000754107000000000	0.200000000000000	0.00300000000000000	0.00395806200000000	0.200000000000000	0.00300000000000000	2.00000000000000e-07	0.687901958000000	0.121168222000000	0.313027298000000	0.0489952670000000	0.00100001900000000	0.200000000000000	0.0152239690000000	0.0241408970000000	0.0402348290000000	0.00629008700000000	0.00925012800000000	0.0100000000000000	400	0.000486864000000000	0.0100000000000000	0.0162962610000000	0.100000000000000	0.000231015000000000	0.000132157000000000	0.0338205090000000	0.0128240120000000	400	0.00100000000000000	9.62000000000000e-06	0.0112685380000000	0.0100000000000000	0.500000000000000	1.20000000000000	1.36000000000000	96.3215089100000	0.00177973800000000	400	0.00100000000000000	0.000500000000000000	0.0100000000000000	9.58000000000000e-05	0.00133002000000000	0.00155000000000000	0.0781000000000000	0.0172000000000000	0.0942000000000000	0;...
0.00623687525965315	5.60000000000000e-05	0.00560000000000000	0.443315639726432	0.200000000000000	0.0286921194929444	0.00254591929512555	0.800000000000000	0.400000000000000	0.000256721177985165	6.51757781357932	0.00284100355864639	0.100000000000000	0.165244239510054	0.100000000000000	0.00143646308044359	0.100000000000000	0.000754106930384025	0.200000000000000	0.00300000000000000	0.00395806205704936	0.200000000000000	0.00300000000000000	2.00000000000000e-07	0.687901957753156	0.121168221998327	0.313027298092336	0.0489952673652474	0.00100001869260733	0.200000000000000	0.0152239690045907	0.0241408974486121	0.0402348290810201	0.00629008702628322	0.00925012797982826	0.0100000000000000	400	0.000486863909742820	0.0100000000000000	0.0162962614462024	0.100000000000000	0.000231015457716723	0.000132157361500284	0.0338205094268733	0.0128240121564962	400	0.00100000000000000	9.62202117475019e-06	0.0112685379318887	0.0100000000000000	0.500000000000000	1.20000000000000	1.36000000000000	96.3215089063843	0.00177973790137018	400	0.00100000000000000	0.000500000000000000	0.0100000000000000	9.58200113926253e-05	0.00133001957198602	0.000774000000000000	0.0290000000000000	0.0205000000000000	0.0142000000000000	0.400000000000000;...
0.00623687525965315	5.60000000000000e-05	0.00560000000000000	0.443315639726432	0.200000000000000	0.0286921194929444	0.00254591929512555	0.800000000000000	0.400000000000000	0.000256721177985165	6.51757781357932	0.00284100355864639	0.100000000000000	0.165244239510054	0.100000000000000	0.00143646308044359	0.100000000000000	0.000754106930384025	0.200000000000000	0.00300000000000000	0.00395806205704936	0.200000000000000	0.00300000000000000	2.00000000000000e-07	0.687901957753156	0.121168221998327	0.313027298092336	0.0489952673652474	0.00100001869260733	0.200000000000000	0.0152239690045907	0.0241408974486121	0.0402348290810201	0.00629008702628322	0.00925012797982826	0.0100000000000000	400	0.000486863909742820	0.0100000000000000	0.0162962614462024	0.100000000000000	0.000231015457716723	0.000132157361500284	0.0338205094268733	0.0128240121564962	400	0.00100000000000000	9.62202117475019e-06	0.0112685379318887	0.0100000000000000	0.500000000000000	1.20000000000000	1.36000000000000	96.3215089063843	0.00177973790137018	400	0.00100000000000000	0.000500000000000000	0.0100000000000000	9.58200113926253e-05	0.00133001957198602	0.000774000000000000	0.0290000000000000	0.0205000000000000	0.0142000000000000	0.00300000000000000;...
0.00623687525965315	5.60000000000000e-05	0.00560000000000000	0.443315639726432	0.200000000000000	0.0286921194929444	0.00254591929512555	0.800000000000000	0.400000000000000	0.000256721177985165	6.51757781357932	0.00284100355864639	0.100000000000000	0.165244239510054	0.100000000000000	0.00143646308044359	0.100000000000000	0.000754106930384025	0.200000000000000	0.00300000000000000	0.00395806205704936	0.200000000000000	0.00300000000000000	2.00000000000000e-07	0.687901957753156	0.121168221998327	0.313027298092336	0.0489952673652474	0.00100001869260733	0.200000000000000	0.0152239690045907	0.0241408974486121	0.0402348290810201	0.00629008702628322	0.00925012797982826	0.0100000000000000	400	0.000486863909742820	0.0100000000000000	0.0162962614462024	0.100000000000000	0.000231015457716723	0.000132157361500284	0.0338205094268733	0.0128240121564962	400	0.00100000000000000	9.62202117475019e-06	0.0112685379318887	0.0100000000000000	0.500000000000000	1.20000000000000	1.36000000000000	96.3215089063843	0.00177973790137018	400	0.00100000000000000	0.000500000000000000	0.0100000000000000	9.58200113926253e-05	0.00133001957198602	0.000774000000000000	0.0290000000000000	0.0205000000000000	0.0142000000000000	0.0162962614462024];

baselineParams = zeros(1,65);
baselineParams(1,1:61) = endosomal_params(1,1:61);

% Define initial values %
PRLR_initvals = [...
9.09000000000000	24.2943543200000	0	60.2009430700000	57.3342314900000	85.8571154300000	55.0106982900000	102.239163600000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	72.2411316800000	77.9745548300000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	13.8803946300000;...
9.09000000000000	18.2207657400000	6.07358858000000	60.2009430700000	57.3342314900000	85.8571154300000	55.0106982900000	102.239163600000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	72.2411316800000	77.9745548300000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	13.8803946300000;...
9.09000000000000	12.1471771600000	12.1471771600000	60.2009430700000	57.3342314900000	85.8571154300000	55.0106982900000	102.239163600000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	72.2411316800000	77.9745548300000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	13.8803946300000;...
9.09000000000000	6.07358858000000	18.2207657400000	60.2009431000000	57.3342315000000	85.8571154000000	55.0106983000000	102.239164000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	72.2411317000000	77.9745548000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	13.8803946000000;...
9.09000000000000	0	24.2943543200000	60.2009430700000	57.3342314900000	85.8571154300000	55.0106982900000	102.239163600000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	72.2411316800000	77.9745548300000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	13.8803946300000];

%PRLR_initvals = importdata('traffickingParams_initValues/endosomal_initvals.mat');
allMembrane_PRLR_initvals = PRLR_initvals(1,:);
membrane_PRLR_75_25_initvals = PRLR_initvals(2,:);
membrane_PRLR_50_50_initvals = PRLR_initvals(3,:);
membrane_PRLR_25_75_initvals = PRLR_initvals(4,:);
allEndosomal_PRLR_initvals = PRLR_initvals(5,:);

numInitialConditions = 5;
initialConditions = zeros(5,67);
initialConditions(1,:) = allEndosomal_PRLR_initvals;
initialConditions(2,:) = membrane_PRLR_25_75_initvals;
initialConditions(3,:) = membrane_PRLR_50_50_initvals;
initialConditions(4,:) = membrane_PRLR_75_25_initvals;
initialConditions(5,:) = allMembrane_PRLR_initvals;

% selected_initvals = allMembrane_PRLR_initvals;
% selected_initvals = membrane_PRLR_75_25_initvals;
% selected_initvals = membrane_PRLR_50_50_initvals;
% selected_initvals = membrane_PRLR_25_75_initvals;
selected_initvals = allEndosomal_PRLR_initvals;

numCalculations = 13;

% Set options %
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'NonNegative',[1:size(selected_initvals,2)]);


% Define storage for outputs %
% Z holds the desired feature to plot as the third dimension – end (6 h) value of PRLRJ2a/iPRLRJ2a
Z_kfree_PRLRJ2aRatio_endValues = zeros(length(kintfreeMesh),length(krecfreeMesh));
Z_kbound_PRLRJ2aRatio_endValues = zeros(length(kintboundMesh),length(krecboundMesh));
Z_kint_PRLRJ2aRatio_endValues = zeros(length(kintfreeMesh),length(kintboundMesh));
Z_krec_PRLRJ2aRatio_endValues = zeros(length(krecfreeMesh),length(krecboundMesh));
Z_kintfree_krecbound_PRLRJ2aRatio_endValues = zeros(length(kintfreeMesh),length(krecboundMesh));
Z_krecfree_kintbound_PRLRJ2aRatio_endValues = zeros(length(krecfreeMesh),length(kintboundMesh));

calculations_kfree = zeros(length(predTime),length(numCalculations));
calculations_kbound = zeros(length(predTime),length(numCalculations));
calculations_kint = zeros(length(predTime),length(numCalculations));
calculations_krec = zeros(length(predTime),length(numCalculations));
calculations_kintfree_krecbound = zeros(length(predTime),length(numCalculations));
calculations_krecfree_kintbound = zeros(length(predTime),length(numCalculations));

% Define cells to store all output structs %
outputCell_kfree = cell(length(kintfreeMesh),length(krecfreeMesh),length(numInitialConditions));    % store all prediction outputs (1441x67 for each param value)
outputCell_kbound = cell(length(kintboundMesh),length(krecboundMesh),length(numInitialConditions)); 
outputCell_kint = cell(length(kintfreeMesh),length(kintboundMesh),length(numInitialConditions)); 
outputCell_krec = cell(length(krecfreeMesh),length(krecboundMesh),length(numInitialConditions));
outputCell_kintfree_krecbound = cell(length(kintfreeMesh),length(krecboundMesh),length(numInitialConditions)); 
outputCell_krecfree_kintbound = cell(length(krecfreeMesh),length(kintboundMesh),length(numInitialConditions)); 

%% Define @core_endosomalCompartment function for running the solver %%
function dydt = core_endosomalCompartment(t,y,params)

% Lynne Cherchia
% 26 August 2025
% Adds endosomal PRLR signaling compartment to RM's base model. Addresses
% ALL terms affected by iPRLRJ species. 
% ** Uses RM'S degradation terms **

% Specify variable names
PRL = y(1);
RJ = y(2);
iRJ = y(3);             % internal (endosomal) RJ
S5Ac = y(4);
S5Bc = y(5);
SHP2 = y(6);
PPX = y(7);
PPN = y(8);
PRLRJ = y(9);
iPRLRJ = y(10);
PRLRJ2 = y(11);
iPRLRJ2 = y(12);
PRLRJ2a = y(13);
iPRLRJ2a = y(14);
PRLRJ2aS5Ac = y(15);
iPRLRJ2aS5Ac = y(16);
PRLRJ2aS5Bc = y(17);
iPRLRJ2aS5Bc = y(18);
pS5Ac = y(19);
pS5Bc = y(20);
pS5AcpS5Ac = y(21);
pS5AcpS5Bc = y(22);
pS5BcpS5Bc = y(23);
PRLRJ2aSHP2 = y(24);
iPRLRJ2aSHP2 = y(25);
PPXpS5Ac = y(26);
PPXpS5Bc = y(27);
PPXpS5AcpS5Ac = y(28);
PPXpS5BcpS5Bc = y(29);
PPXpS5AcpS5Bc = y(30);
pS5AcS5Ac = y(31);
pS5BcS5Bc = y(32);
S5AcpS5Bc = y(33);
pS5AcS5Bc = y(34);
pS5AnpS5An = y(35);
pS5AnpS5Bn = y(36);
pS5BnpS5Bn = y(37);
pS5An = y(38);
pS5Bn = y(39);
PPNpS5An = y(40);
PPNpS5Bn = y(41);
S5An = y(42);
S5Bn = y(43);
PPNpS5AnpS5An = y(44);
PPNpS5AnpS5Bn = y(45);
PPNpS5BnpS5Bn = y(46);
pS5AnS5An = y(47);
pS5BnS5Bn = y(48);
S5AnpS5Bn = y(49);
pS5AnS5Bn = y(50);
mRNAn = y(51);
mRNAc = y(52);
SOCS1 = y(53);
SOCS1PRLRJ2a = y(54);
iSOCS1PRLRJ2a = y(55);
PRLRJ2aS5AcSHP2 = y(56);
iPRLRJ2aS5AcSHP2 = y(57);
PRLRJ2aS5BcSHP2 = y(58);
iPRLRJ2aS5BcSHP2 = y(59);
SOCS1PRLRJ2aSHP2 = y(60);
iSOCS1PRLRJ2aSHP2 = y(61);
mRn = y(62);
mRc = y(63);
Rc = y(64);
mBCLn = y(65);
mBCLc = y(66);
BCL = y(67);


% Specify Parameter names
k1 = params(1);
k2 = params(2);
k_2 = params(3);
k3 = params(4);
k_3 = params(5);
k4 = params(6);
k5 = params(7);
k_5 = params(8);
k6 = params(9);
kdeg = params(10);
deg_ratio = params(11);
k8A = params(12);
k_8A = params(13);
k8B = params(14);
k_8B = params(15);
k8AB = params(16);
k_8AB = params(17);
k9 = params(18);
k_9 = params(19);
k10 = params(20);
k11 = params(21);
k_11 = params(22);
k12 = params(23);
k13 = params(24);
k_13 = params(25);
k14A = params(26);
k14B = params(27);
k14AB = params(28);
k15 = params(29);
k_15 = params(30);
k16 = params(31);
k17inA = params(32);
k17outA = params(33);
k17inB = params(34);
k17outB = params(35);
k18a = params(36);
k18b = params(37);
k19 = params(38);
k20 = params(39);
k21 = params(40);
k_21 = params(41);
k22 = params(42);
k23 = params(43);
k24 = params(44);
k25a = params(45);
k25b = params(46);
k26 = params(47);
k27 = params(48);
k28 = params(49);
k29 = params(50);
Vratio = params(51);
ncratioA = params(52);
ncratioB = params(53);
totalSTAT = params(54);
k30a = params(55);
k30b = params(56);
k31 = params(57);
k32 = params(58);
k33 = params(59);
k34 = params(60);
k35 = params(61);
kint_free = params(62);            % internalization of unbound receptor complexes to endosomes
krec_free = params(63);            % recycling of unbound receptor complexes to membrane
kint_bound = params(64);           % internalization of PRL-bound receptor complexes to endosomes
krec_bound = params(65);           % recycling of PRL-bound receptor complexes to membrane



% Specify ODEs
dydt(1,1) = (k_2*PRLRJ - k2*PRL*RJ)*1.39E-4 - (1 * 0.693 / (6 * 3600) * y(1)); % Adjusted for volume of cell to volume of cytoplasm ratio, 1st order deg with 6 hour half life

dydt(2,1) = k1 + k_2*PRLRJ - k2*PRL*RJ - k3*PRLRJ*RJ + k_3*PRLRJ2 - kdeg*RJ + k29*Rc - kint_free*RJ + krec_free*iRJ;
dydt(3,1) = k_3*iPRLRJ2 - k3*iPRLRJ*iRJ - kdeg*iRJ + kint_free*RJ - krec_free*iRJ;

dydt(4,1) = k_5*PRLRJ2aS5Ac + k_5*iPRLRJ2aS5Ac - k5*S5Ac*PRLRJ2a - k5*iPRLRJ2a*S5Ac + k12*PPXpS5Ac - k13*S5Ac*pS5Ac + k_13*pS5AcS5Ac - k13*S5Ac*pS5Bc + k_13*S5AcpS5Bc + k10*PRLRJ2aS5AcSHP2 + k10*iPRLRJ2aS5AcSHP2 + kdeg*deg_ratio*(PRLRJ2aS5Ac + PRLRJ2aS5AcSHP2) + kdeg*deg_ratio*(iPRLRJ2aS5Ac + iPRLRJ2aS5AcSHP2) - k17inA*S5Ac + k17outA*S5An.*Vratio;
dydt(5,1) = k_5*PRLRJ2aS5Bc + k_5*iPRLRJ2aS5Bc - k5*S5Bc*PRLRJ2a - k5*iPRLRJ2a*S5Bc + k12*PPXpS5Bc - k13*S5Bc*pS5Bc + k_13*pS5BcS5Bc - k13*S5Bc*pS5Ac + k_13*pS5AcS5Bc + k10*PRLRJ2aS5BcSHP2 + k10*iPRLRJ2aS5BcSHP2 + kdeg*deg_ratio*(PRLRJ2aS5Bc + PRLRJ2aS5BcSHP2) + kdeg*deg_ratio*(iPRLRJ2aS5Bc + iPRLRJ2aS5BcSHP2) - k17inB*S5Bc + k17outB*S5Bn.*Vratio;
dydt(6,1) = k_9*PRLRJ2aSHP2 + k_9*iPRLRJ2aSHP2 - k9*PRLRJ2a*SHP2 - k9*iPRLRJ2a*SHP2 + k10*PRLRJ2aSHP2 + k10*iPRLRJ2aSHP2 - k9*SHP2*SOCS1PRLRJ2a - k9*SHP2*iSOCS1PRLRJ2a + k_9*SOCS1PRLRJ2aSHP2 + k_9*iSOCS1PRLRJ2aSHP2 + k10*SOCS1PRLRJ2aSHP2 + k10*iSOCS1PRLRJ2aSHP2 - k9*SHP2*PRLRJ2aS5Ac - k9*SHP2*iPRLRJ2aS5Ac + k_9*PRLRJ2aS5AcSHP2 + k_9*iPRLRJ2aS5AcSHP2 - k9*SHP2*PRLRJ2aS5Bc - k9*SHP2*iPRLRJ2aS5Bc + k_9*PRLRJ2aS5BcSHP2 + k_9*iPRLRJ2aS5BcSHP2 + k10*PRLRJ2aS5AcSHP2 + k10*iPRLRJ2aS5AcSHP2 + k10*PRLRJ2aS5BcSHP2 + k10*iPRLRJ2aS5BcSHP2 + kdeg*deg_ratio*(PRLRJ2aSHP2 + PRLRJ2aS5AcSHP2 + PRLRJ2aS5BcSHP2 + SOCS1PRLRJ2aSHP2) + kdeg*deg_ratio*(iPRLRJ2aSHP2 + iPRLRJ2aS5AcSHP2 + iPRLRJ2aS5BcSHP2 + iSOCS1PRLRJ2aSHP2) + k24*SOCS1PRLRJ2aSHP2 + k24*iSOCS1PRLRJ2aSHP2;

dydt(7,1) = k_11*PPXpS5Ac - k11*PPX*pS5Ac + k12*PPXpS5Ac + k_11*PPXpS5Bc - k11*PPX*pS5Bc + k12*PPXpS5Bc - k11*PPX*pS5AcpS5Ac + k_11*PPXpS5AcpS5Ac + k12*PPXpS5AcpS5Ac - k11*PPX*pS5BcpS5Bc + k_11*PPXpS5BcpS5Bc + k12*PPXpS5BcpS5Bc - k11*PPX*pS5AcpS5Bc + k_11*PPXpS5AcpS5Bc + k12*PPXpS5AcpS5Bc;
dydt(8,1) = k_15*PPNpS5An - k15*PPN*pS5An + k16*PPNpS5An + k_15*PPNpS5Bn - k15*PPN*pS5Bn + k16*PPNpS5Bn - k15*PPN*pS5AnpS5An + k_15*PPNpS5AnpS5An + k16*PPNpS5AnpS5An - k15*PPN*pS5BnpS5Bn + k_15*PPNpS5BnpS5Bn + k16*PPNpS5BnpS5Bn - k15*PPN*pS5AnpS5Bn + k_15*PPNpS5AnpS5Bn + k16*PPNpS5AnpS5Bn;

dydt(9,1) =  k2*PRL*RJ - k_2*PRLRJ - k3*PRLRJ*RJ + k_3*PRLRJ2 - kdeg*deg_ratio*PRLRJ - kint_bound*PRLRJ + krec_bound*iPRLRJ;
dydt(10,1) = k_3*iPRLRJ2 - k3*iPRLRJ*iRJ - kdeg*deg_ratio*iPRLRJ + kint_bound*PRLRJ - krec_bound*iPRLRJ;

dydt(11,1) = k3*PRLRJ*RJ - k_3*PRLRJ2 - k4*PRLRJ2 + k10*PRLRJ2aSHP2 + k10*SOCS1PRLRJ2aSHP2 + k10*PRLRJ2aS5AcSHP2 + k10*PRLRJ2aS5BcSHP2 - kdeg*deg_ratio*PRLRJ2 - kint_bound*PRLRJ2 + krec_bound*iPRLRJ2;
dydt(12,1) = k3*iPRLRJ*iRJ - k_3*iPRLRJ2 - k4*iPRLRJ2 + k10*iPRLRJ2aSHP2 + k10*iSOCS1PRLRJ2aSHP2 + k10*iPRLRJ2aS5AcSHP2 + k10*iPRLRJ2aS5BcSHP2 - kdeg*deg_ratio*iPRLRJ2 + kint_bound*PRLRJ2 - krec_bound*iPRLRJ2;

dydt(13,1) = k4*PRLRJ2 - k5*PRLRJ2a*S5Ac + k_5*PRLRJ2aS5Ac + k6*PRLRJ2aS5Ac - k5*PRLRJ2a*S5Bc + k_5*PRLRJ2aS5Bc + k6*PRLRJ2aS5Bc - k9*PRLRJ2a*SHP2 + k_9*PRLRJ2aSHP2 - k21*PRLRJ2a*SOCS1 + k_21*SOCS1PRLRJ2a + k23*SOCS1PRLRJ2a - kdeg*deg_ratio*PRLRJ2a - kint_bound*PRLRJ2a + krec_bound*iPRLRJ2a;
dydt(14,1) = k4*iPRLRJ2 - k5*iPRLRJ2a*S5Ac + k_5*iPRLRJ2aS5Ac + k6*iPRLRJ2aS5Ac - k5*iPRLRJ2a*S5Bc + k_5*iPRLRJ2aS5Bc + k6*iPRLRJ2aS5Bc - k9*iPRLRJ2a*SHP2 + k_9*iPRLRJ2aSHP2 - k21*iPRLRJ2a*SOCS1 + k_21*iSOCS1PRLRJ2a + k23*iSOCS1PRLRJ2a - kdeg*deg_ratio*iPRLRJ2a + kint_bound*PRLRJ2a - krec_bound*iPRLRJ2a;

dydt(15,1) = k5*PRLRJ2a*S5Ac - k_5*PRLRJ2aS5Ac - k6*PRLRJ2aS5Ac - k9*SHP2*PRLRJ2aS5Ac + k_9*PRLRJ2aS5AcSHP2 - kdeg*deg_ratio*PRLRJ2aS5Ac - kint_bound*PRLRJ2aS5Ac + krec_bound*iPRLRJ2aS5Ac;
dydt(16,1) = k5*iPRLRJ2a*S5Ac - k_5*iPRLRJ2aS5Ac - k6*iPRLRJ2aS5Ac - k9*SHP2*iPRLRJ2aS5Ac + k_9*iPRLRJ2aS5AcSHP2 - kdeg*deg_ratio*iPRLRJ2aS5Ac + kint_bound*PRLRJ2aS5Ac - krec_bound*iPRLRJ2aS5Ac;

dydt(17,1) = k5*PRLRJ2a*S5Bc - k_5*PRLRJ2aS5Bc - k6*PRLRJ2aS5Bc - k9*SHP2*PRLRJ2aS5Bc + k_9*PRLRJ2aS5BcSHP2 - kdeg*deg_ratio*PRLRJ2aS5Bc - kint_bound*PRLRJ2aS5Bc + krec_bound*iPRLRJ2aS5Bc;
dydt(18,1) = k5*iPRLRJ2a*S5Bc - k_5*iPRLRJ2aS5Bc - k6*iPRLRJ2aS5Bc - k9*SHP2*iPRLRJ2aS5Bc + k_9*iPRLRJ2aS5BcSHP2 - kdeg*deg_ratio*iPRLRJ2aS5Bc + kint_bound*PRLRJ2aS5Bc - krec_bound*iPRLRJ2aS5Bc;

dydt(19,1) = k6*PRLRJ2aS5Ac + k6*iPRLRJ2aS5Ac - 2*k8A*pS5Ac*pS5Ac + 2*k_8A*pS5AcpS5Ac - k8AB*pS5Ac*pS5Bc + k_8AB*pS5AcpS5Bc - k11*PPX*pS5Ac + k_11*PPXpS5Ac - k13*pS5Ac*S5Ac + k_13*pS5AcS5Ac - k13*pS5Ac*S5Bc + k_13*pS5AcS5Bc;
dydt(20,1) = k6*PRLRJ2aS5Bc + k6*iPRLRJ2aS5Bc - 2*k8B*pS5Bc*pS5Bc + 2*k_8B*pS5BcpS5Bc - k8AB*pS5Ac*pS5Bc + k_8AB*pS5AcpS5Bc - k11*PPX*pS5Bc + k_11*PPXpS5Bc - k13*pS5Bc*S5Bc + k_13*pS5BcS5Bc - k13*pS5Bc*S5Ac + k_13*S5AcpS5Bc;

dydt(21,1) = k8A*pS5Ac*pS5Ac - k_8A*pS5AcpS5Ac - k11*PPX*pS5AcpS5Ac + k_11*PPXpS5AcpS5Ac - k14A*pS5AcpS5Ac;
dydt(22,1) = k8AB*pS5Ac*pS5Bc - k_8AB*pS5AcpS5Bc - k11*PPX*pS5AcpS5Bc + k_11*PPXpS5AcpS5Bc - k14AB*pS5AcpS5Bc;
dydt(23,1) = k8B*pS5Bc*pS5Bc - k_8B*pS5BcpS5Bc - k11*PPX*pS5BcpS5Bc + k_11*PPXpS5BcpS5Bc - k14B*pS5BcpS5Bc;

dydt(24,1) = k9*PRLRJ2a*SHP2 - k_9*PRLRJ2aSHP2 - k10*PRLRJ2aSHP2 + k23*SOCS1PRLRJ2aSHP2 - kdeg*deg_ratio*PRLRJ2aSHP2 - kint_bound*PRLRJ2aSHP2 + krec_bound*iPRLRJ2aSHP2;
dydt(25,1) = k9*iPRLRJ2a*SHP2 - k_9*iPRLRJ2aSHP2 - k10*iPRLRJ2aSHP2 + k23*iSOCS1PRLRJ2aSHP2 - kdeg*deg_ratio*iPRLRJ2aSHP2 + kint_bound*PRLRJ2aSHP2 - krec_bound*iPRLRJ2aSHP2;

dydt(26,1) = k11*PPX*pS5Ac - k_11*PPXpS5Ac - k12*PPXpS5Ac;
dydt(27,1) = k11*PPX*pS5Bc - k_11*PPXpS5Bc - k12*PPXpS5Bc;
dydt(28,1) = k11*PPX*pS5AcpS5Ac - k_11*PPXpS5AcpS5Ac - k12*PPXpS5AcpS5Ac;
dydt(29,1) = k11*PPX*pS5BcpS5Bc - k_11*PPXpS5BcpS5Bc - k12*PPXpS5BcpS5Bc;
dydt(30,1) = k11*PPX*pS5AcpS5Bc - k_11*PPXpS5AcpS5Bc - k12*PPXpS5AcpS5Bc;
dydt(31,1) = k13*pS5Ac*S5Ac - k_13*pS5AcS5Ac + k12*PPXpS5AcpS5Ac;
dydt(32,1) = k13*pS5Bc*S5Bc - k_13*pS5BcS5Bc + k12*PPXpS5BcpS5Bc;
dydt(33,1) = k13*S5Ac*pS5Bc - k_13*S5AcpS5Bc + 0.5*k12*PPXpS5AcpS5Bc;
dydt(34,1) = k13*pS5Ac*S5Bc - k_13*pS5AcS5Bc + 0.5*k12*PPXpS5AcpS5Bc;
dydt(35,1) = k14A*pS5AcpS5Ac./Vratio + k8A*pS5An*pS5An - k_8A*pS5AnpS5An - k15*PPN*pS5AnpS5An + k_15*PPNpS5AnpS5An;
dydt(36,1) = k14AB*pS5AcpS5Bc./Vratio + k8AB*pS5An*pS5Bn - k_8AB*pS5AnpS5Bn - k15*PPN*pS5AnpS5Bn + k_15*PPNpS5AnpS5Bn;
dydt(37,1) = k14B*pS5BcpS5Bc./Vratio + k8B*pS5Bn*pS5Bn - k_8B*pS5BnpS5Bn - k15*PPN*pS5BnpS5Bn + k_15*PPNpS5BnpS5Bn;
dydt(38,1) = 2*k_8A*pS5AnpS5An - 2*k8A*pS5An*pS5An + k_8AB*pS5AnpS5Bn - k8AB*pS5An*pS5Bn - k15*PPN*pS5An + k_15*PPNpS5An - k13*pS5An*S5An + k_13*pS5AnS5An - k13*pS5An*S5Bn + k_13*pS5AnS5Bn;
dydt(39,1) = 2*k_8B*pS5BnpS5Bn - 2*k8B*pS5Bn*pS5Bn + k_8AB*pS5AnpS5Bn - k8AB*pS5An*pS5Bn - k15*PPN*pS5Bn + k_15*PPNpS5Bn - k13*pS5Bn*S5Bn + k_13*pS5BnS5Bn - k13*pS5Bn*S5An + k_13*S5AnpS5Bn;
dydt(40,1) = k15*PPN*pS5An - k_15*PPNpS5An - k16*PPNpS5An;
dydt(41,1) = k15*PPN*pS5Bn - k_15*PPNpS5Bn - k16*PPNpS5Bn;
dydt(42,1) = k16*PPNpS5An - k13*pS5An*S5An + k_13*pS5AnS5An - k13*pS5Bn*S5An + k_13*S5AnpS5Bn + k17inA*S5Ac./Vratio - k17outA*S5An;
dydt(43,1) = k16*PPNpS5Bn - k13*pS5Bn*S5Bn + k_13*pS5BnS5Bn - k13*pS5An*S5Bn + k_13*pS5AnS5Bn + k17inB*S5Bc./Vratio - k17outB*S5Bn;
dydt(44,1) = k15*PPN*pS5AnpS5An - k_15*PPNpS5AnpS5An - k16*PPNpS5AnpS5An;
dydt(45,1) = k15*PPN*pS5AnpS5Bn - k_15*PPNpS5AnpS5Bn - k16*PPNpS5AnpS5Bn;
dydt(46,1) = k15*PPN*pS5BnpS5Bn - k_15*PPNpS5BnpS5Bn - k16*PPNpS5BnpS5Bn;
dydt(47,1) = k16*PPNpS5AnpS5An + k13*pS5An*S5An - k_13*pS5AnS5An;
dydt(48,1) = k16*PPNpS5BnpS5Bn + k13*pS5Bn*S5Bn - k_13*pS5BnS5Bn;
dydt(49,1) = 0.5*k16*PPNpS5AnpS5Bn + k13*S5An*pS5Bn - k_13*S5AnpS5Bn;
dydt(50,1) = 0.5*k16*PPNpS5AnpS5Bn + k13*pS5An*S5Bn - k_13*pS5AnS5Bn;
dydt(51,1) = k18a*(pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn)/(k18b+(pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn)) - k19*mRNAn;
dydt(52,1) = k19*mRNAn.*Vratio - k22*mRNAc;

dydt(53,1) = k20*mRNAc- k21*SOCS1*PRLRJ2a - k21*SOCS1*iPRLRJ2a + k_21*SOCS1PRLRJ2a + k_21*iSOCS1PRLRJ2a - k23*SOCS1 + k10*SOCS1PRLRJ2aSHP2 + k10*iSOCS1PRLRJ2aSHP2 + kdeg*deg_ratio*(SOCS1PRLRJ2a + SOCS1PRLRJ2aSHP2) + kdeg*deg_ratio*(iSOCS1PRLRJ2a + iSOCS1PRLRJ2aSHP2);

dydt(54,1) = k21*SOCS1*PRLRJ2a - k_21*SOCS1PRLRJ2a - k9*SHP2*SOCS1PRLRJ2a + k_9*SOCS1PRLRJ2aSHP2 - k23*SOCS1PRLRJ2a - kdeg*deg_ratio*SOCS1PRLRJ2a - k24*SOCS1PRLRJ2a - kint_bound*SOCS1PRLRJ2a + krec_bound*iSOCS1PRLRJ2a;
dydt(55,1) = k21*SOCS1*iPRLRJ2a - k_21*iSOCS1PRLRJ2a - k9*SHP2*iSOCS1PRLRJ2a + k_9*iSOCS1PRLRJ2aSHP2 - k23*iSOCS1PRLRJ2a - kdeg*deg_ratio*iSOCS1PRLRJ2a - k24*iSOCS1PRLRJ2a + kint_bound*SOCS1PRLRJ2a - krec_bound*iSOCS1PRLRJ2a;

dydt(56,1) = k9*SHP2*PRLRJ2aS5Ac - k_9*PRLRJ2aS5AcSHP2 - k10*PRLRJ2aS5AcSHP2 - kdeg*deg_ratio*PRLRJ2aS5AcSHP2 - kint_bound*PRLRJ2aS5AcSHP2 + krec_bound*PRLRJ2aS5AcSHP2;
dydt(57,1) = k9*SHP2*iPRLRJ2aS5Ac - k_9*iPRLRJ2aS5AcSHP2 - k10*iPRLRJ2aS5AcSHP2 - kdeg*deg_ratio*iPRLRJ2aS5AcSHP2 + kint_bound*PRLRJ2aS5AcSHP2 - krec_bound*PRLRJ2aS5AcSHP2;

dydt(58,1) = k9*SHP2*PRLRJ2aS5Bc - k_9*PRLRJ2aS5BcSHP2 - k10*PRLRJ2aS5BcSHP2 - kdeg*deg_ratio*PRLRJ2aS5BcSHP2 - kint_bound*PRLRJ2aS5BcSHP2 + krec_bound*PRLRJ2aS5BcSHP2;
dydt(59,1) = k9*SHP2*iPRLRJ2aS5Bc - k_9*iPRLRJ2aS5BcSHP2 - k10*iPRLRJ2aS5BcSHP2 - kdeg*deg_ratio*iPRLRJ2aS5BcSHP2 + kint_bound*PRLRJ2aS5BcSHP2 - krec_bound*PRLRJ2aS5BcSHP2;

dydt(60,1) = k9*SHP2*SOCS1PRLRJ2a - k_9*SOCS1PRLRJ2aSHP2 - k10*SOCS1PRLRJ2aSHP2 - k23*SOCS1PRLRJ2aSHP2 - kdeg*deg_ratio*SOCS1PRLRJ2aSHP2 - k24*SOCS1PRLRJ2aSHP2 - kint_bound*SOCS1PRLRJ2aSHP2 + krec_bound*iSOCS1PRLRJ2aSHP2;
dydt(61,1) = k9*SHP2*iSOCS1PRLRJ2a - k_9*iSOCS1PRLRJ2aSHP2 - k10*iSOCS1PRLRJ2aSHP2 - k23*iSOCS1PRLRJ2aSHP2 - kdeg*deg_ratio*iSOCS1PRLRJ2aSHP2 - k24*iSOCS1PRLRJ2aSHP2 + kint_bound*SOCS1PRLRJ2aSHP2 - krec_bound*iSOCS1PRLRJ2aSHP2;

dydt(62,1) = k25a*(pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn)/(k25b+(pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn)) - k26*mRn;
dydt(63,1) = k26*mRn.*Vratio - k27*mRc;
dydt(64,1) = k28*mRc - k29*Rc;
dydt(65,1) = k30a*(pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn)/(k30b+(pS5AnpS5An + pS5AnpS5Bn + pS5BnpS5Bn)) - k31*mBCLn;
dydt(66,1) = k31*mBCLn.*Vratio - k32*mBCLc;
dydt(67,1) = k33*mBCLc - k34*BCL + k35;

end


%% Loop initial conditions %%

% Simulate varying kfree cycling rate values %
tic
for k = 1:numInitialConditions
    for i = 1:length(kintfreeMesh)
        for j = 1:length(krecfreeMesh)
            params = baselineParams;
            params(62) = kintfreeMesh(i);
            params(63) = krecfreeMesh(j);
            params(64) = kintboundBaseline;
            params(65) = krecboundBaseline;
            initvalues = initialConditions(k,:);

            [~, predConc_kfree] = ode15s(@core_endosomalCompartment,predTime,initvalues,options,params);

            PRL = predConc_kfree(:,1);
            RJ = predConc_kfree(:,2);
            iRJ = predConc_kfree(:,3);
            S5Ac = predConc_kfree(:,4);
            S5Bc = predConc_kfree(:,5);
            SHP2 = predConc_kfree(:,6);
            PPX = predConc_kfree(:,7);
            PPN = predConc_kfree(:,8);
            PRLRJ = predConc_kfree(:,9);
            iPRLRJ = predConc_kfree(:,10);
            PRLRJ2 = predConc_kfree(:,11);
            iPRLRJ2 = predConc_kfree(:,12);
            PRLRJ2a = predConc_kfree(:,13);
            iPRLRJ2a = predConc_kfree(:,14);
            PRLRJ2aS5Ac = predConc_kfree(:,15);
            iPRLRJ2aS5Ac = predConc_kfree(:,16);
            PRLRJ2aS5Bc = predConc_kfree(:,17);
            iPRLRJ2aS5Bc = predConc_kfree(:,18);
            pS5Ac = predConc_kfree(:,19);
            pS5Bc = predConc_kfree(:,20);
            pS5AcpS5Ac = predConc_kfree(:,21);
            pS5AcpS5Bc = predConc_kfree(:,22);
            pS5BcpS5Bc = predConc_kfree(:,23);
            PRLRJ2aSHP2 = predConc_kfree(:,24);
            iPRLRJ2aSHP2 = predConc_kfree(:,25);
            PPXpS5Ac = predConc_kfree(:,26);
            PPXpS5Bc = predConc_kfree(:,27);
            PPXpS5AcpS5Ac = predConc_kfree(:,28);
            PPXpS5BcpS5Bc = predConc_kfree(:,29);
            PPXpS5AcpS5Bc = predConc_kfree(:,30);
            pS5AcS5Ac = predConc_kfree(:,31);
            pS5BcS5Bc = predConc_kfree(:,32);
            S5AcpS5Bc = predConc_kfree(:,33);
            pS5AcS5Bc = predConc_kfree(:,34);
            pS5AnpS5An = predConc_kfree(:,35);
            pS5AnpS5Bn = predConc_kfree(:,36);
            pS5BnpS5Bn = predConc_kfree(:,37);
            pS5An = predConc_kfree(:,38);
            pS5Bn = predConc_kfree(:,39);
            PPNpS5An = predConc_kfree(:,40);
            PPNpS5Bn = predConc_kfree(:,41);
            S5An = predConc_kfree(:,42);
            S5Bn = predConc_kfree(:,43);
            PPNpS5AnpS5An = predConc_kfree(:,44);
            PPNpS5AnpS5Bn = predConc_kfree(:,45);
            PPNpS5BnpS5Bn = predConc_kfree(:,46);
            pS5AnS5An = predConc_kfree(:,47);
            pS5BnS5Bn = predConc_kfree(:,48);
            S5AnpS5Bn = predConc_kfree(:,49);
            pS5AnS5Bn = predConc_kfree(:,50);
            mRNAn = predConc_kfree(:,51);
            mRNAc = predConc_kfree(:,52);
            SOCS1 = predConc_kfree(:,53);
            SOCS1PRLRJ2a = predConc_kfree(:,54);
            iSOCS1PRLRJ2a = predConc_kfree(:,55);
            PRLRJ2aS5AcSHP2 = predConc_kfree(:,56);
            iPRLRJ2aS5AcSHP2 = predConc_kfree(:,57);
            PRLRJ2aS5BcSHP2 = predConc_kfree(:,58);
            iPRLRJ2aS5BcSHP2 = predConc_kfree(:,59);
            SOCS1PRLRJ2aSHP2 = predConc_kfree(:,60);
            iSOCS1PRLRJ2aSHP2 = predConc_kfree(:,61);
            mRn = predConc_kfree(:,62);
            mRc = predConc_kfree(:,63);
            Rc = predConc_kfree(:,64);
            mBCLn = predConc_kfree(:,65);
            mBCLc = predConc_kfree(:,66);
            BCL = predConc_kfree(:,67);

            % Store each predConc entry at the correct index of the cell
            outputCell_kfree{i,j,k} = predConc_kfree(:,:);

            % Calculate and store end value ratio of PRLRJ2a/iPRLRJ2a
            ratio_kfree = PRLRJ2a./iPRLRJ2a;
            ratio_endValue_kfree = ratio_kfree(end);
            Z_kfree_PRLRJ2aRatio_endValues(i,j,k) = ratio_endValue_kfree;
        end
    end
end
toc

Z_kfree_1 = Z_kfree_PRLRJ2aRatio_endValues(:,:,1);
Z_kfree_2 = Z_kfree_PRLRJ2aRatio_endValues(:,:,2);
Z_kfree_3 = Z_kfree_PRLRJ2aRatio_endValues(:,:,3);
Z_kfree_4 = Z_kfree_PRLRJ2aRatio_endValues(:,:,4);
Z_kfree_5 = Z_kfree_PRLRJ2aRatio_endValues(:,:,5);

% Simulate varying kbound cycling rate values %
tic
for k = 1:numInitialConditions
    for i = 1:length(kintboundMesh)
        for j = 1:length(krecboundMesh)
            params = baselineParams;
            params(62) = kintfreeBaseline;
            params(63) = krecfreeBaseline;
            params(64) = kintboundMesh(i);
            params(65) = krecboundMesh(j);
            initvalues = initialConditions(k,:);

            [~, predConc_kbound] = ode15s(@core_endosomalCompartment,predTime,initvalues,options,params);

            PRL = predConc_kbound(:,1);
            RJ = predConc_kbound(:,2);
            iRJ = predConc_kbound(:,3);
            S5Ac = predConc_kbound(:,4);
            S5Bc = predConc_kbound(:,5);
            SHP2 = predConc_kbound(:,6);
            PPX = predConc_kbound(:,7);
            PPN = predConc_kbound(:,8);
            PRLRJ = predConc_kbound(:,9);
            iPRLRJ = predConc_kbound(:,10);
            PRLRJ2 = predConc_kbound(:,11);
            iPRLRJ2 = predConc_kbound(:,12);
            PRLRJ2a = predConc_kbound(:,13);
            iPRLRJ2a = predConc_kbound(:,14);
            PRLRJ2aS5Ac = predConc_kbound(:,15);
            iPRLRJ2aS5Ac = predConc_kbound(:,16);
            PRLRJ2aS5Bc = predConc_kbound(:,17);
            iPRLRJ2aS5Bc = predConc_kbound(:,18);
            pS5Ac = predConc_kbound(:,19);
            pS5Bc = predConc_kbound(:,20);
            pS5AcpS5Ac = predConc_kbound(:,21);
            pS5AcpS5Bc = predConc_kbound(:,22);
            pS5BcpS5Bc = predConc_kbound(:,23);
            PRLRJ2aSHP2 = predConc_kbound(:,24);
            iPRLRJ2aSHP2 = predConc_kbound(:,25);
            PPXpS5Ac = predConc_kbound(:,26);
            PPXpS5Bc = predConc_kbound(:,27);
            PPXpS5AcpS5Ac = predConc_kbound(:,28);
            PPXpS5BcpS5Bc = predConc_kbound(:,29);
            PPXpS5AcpS5Bc = predConc_kbound(:,30);
            pS5AcS5Ac = predConc_kbound(:,31);
            pS5BcS5Bc = predConc_kbound(:,32);
            S5AcpS5Bc = predConc_kbound(:,33);
            pS5AcS5Bc = predConc_kbound(:,34);
            pS5AnpS5An = predConc_kbound(:,35);
            pS5AnpS5Bn = predConc_kbound(:,36);
            pS5BnpS5Bn = predConc_kbound(:,37);
            pS5An = predConc_kbound(:,38);
            pS5Bn = predConc_kbound(:,39);
            PPNpS5An = predConc_kbound(:,40);
            PPNpS5Bn = predConc_kbound(:,41);
            S5An = predConc_kbound(:,42);
            S5Bn = predConc_kbound(:,43);
            PPNpS5AnpS5An = predConc_kbound(:,44);
            PPNpS5AnpS5Bn = predConc_kbound(:,45);
            PPNpS5BnpS5Bn = predConc_kbound(:,46);
            pS5AnS5An = predConc_kbound(:,47);
            pS5BnS5Bn = predConc_kbound(:,48);
            S5AnpS5Bn = predConc_kbound(:,49);
            pS5AnS5Bn = predConc_kbound(:,50);
            mRNAn = predConc_kbound(:,51);
            mRNAc = predConc_kbound(:,52);
            SOCS1 = predConc_kbound(:,53);
            SOCS1PRLRJ2a = predConc_kbound(:,54);
            iSOCS1PRLRJ2a = predConc_kbound(:,55);
            PRLRJ2aS5AcSHP2 = predConc_kbound(:,56);
            iPRLRJ2aS5AcSHP2 = predConc_kbound(:,57);
            PRLRJ2aS5BcSHP2 = predConc_kbound(:,58);
            iPRLRJ2aS5BcSHP2 = predConc_kbound(:,59);
            SOCS1PRLRJ2aSHP2 = predConc_kbound(:,60);
            iSOCS1PRLRJ2aSHP2 = predConc_kbound(:,61);
            mRn = predConc_kbound(:,62);
            mRc = predConc_kbound(:,63);
            Rc = predConc_kbound(:,64);
            mBCLn = predConc_kbound(:,65);
            mBCLc = predConc_kbound(:,66);
            BCL = predConc_kbound(:,67);

            % Store each predConc entry at the correct index of the cell
            outputCell_kbound{i,j,k} = predConc_kbound(:,:);
            
            % Calculate and store end value ratio of PRLRJ2a/iPRLRJ2a
            ratio_kbound = PRLRJ2a./iPRLRJ2a;
            ratio_endValue_kbound = ratio_kbound(end);
            Z_kbound_PRLRJ2aRatio_endValues(i,j,k) = ratio_endValue_kbound;
        end
    end
end
toc

Z_kbound_1 = Z_kbound_PRLRJ2aRatio_endValues(:,:,1);
Z_kbound_2 = Z_kbound_PRLRJ2aRatio_endValues(:,:,2);
Z_kbound_3 = Z_kbound_PRLRJ2aRatio_endValues(:,:,3);
Z_kbound_4 = Z_kbound_PRLRJ2aRatio_endValues(:,:,4);
Z_kbound_5 = Z_kbound_PRLRJ2aRatio_endValues(:,:,5);

% Simulate varying kint cycling rate values %
tic
for k = 1:numInitialConditions
    for i = 1:length(kintfreeMesh)
        for j = 1:length(kintboundMesh)
            params = baselineParams;
            params(62) = kintfreeMesh(i);
            params(63) = krecfreeBaseline;
            params(64) = kintboundMesh(j);
            params(65) = krecboundBaseline;
            initvalues = initialConditions(k,:);

            [~, predConc_kint] = ode15s(@core_endosomalCompartment,predTime,initvalues,options,params);

            PRL = predConc_kint(:,1);
            RJ = predConc_kint(:,2);
            iRJ = predConc_kint(:,3);
            S5Ac = predConc_kint(:,4);
            S5Bc = predConc_kint(:,5);
            SHP2 = predConc_kint(:,6);
            PPX = predConc_kint(:,7);
            PPN = predConc_kint(:,8);
            PRLRJ = predConc_kint(:,9);
            iPRLRJ = predConc_kint(:,10);
            PRLRJ2 = predConc_kint(:,11);
            iPRLRJ2 = predConc_kint(:,12);
            PRLRJ2a = predConc_kint(:,13);
            iPRLRJ2a = predConc_kint(:,14);
            PRLRJ2aS5Ac = predConc_kint(:,15);
            iPRLRJ2aS5Ac = predConc_kint(:,16);
            PRLRJ2aS5Bc = predConc_kint(:,17);
            iPRLRJ2aS5Bc = predConc_kint(:,18);
            pS5Ac = predConc_kint(:,19);
            pS5Bc = predConc_kint(:,20);
            pS5AcpS5Ac = predConc_kint(:,21);
            pS5AcpS5Bc = predConc_kint(:,22);
            pS5BcpS5Bc = predConc_kint(:,23);
            PRLRJ2aSHP2 = predConc_kint(:,24);
            iPRLRJ2aSHP2 = predConc_kint(:,25);
            PPXpS5Ac = predConc_kint(:,26);
            PPXpS5Bc = predConc_kint(:,27);
            PPXpS5AcpS5Ac = predConc_kint(:,28);
            PPXpS5BcpS5Bc = predConc_kint(:,29);
            PPXpS5AcpS5Bc = predConc_kint(:,30);
            pS5AcS5Ac = predConc_kint(:,31);
            pS5BcS5Bc = predConc_kint(:,32);
            S5AcpS5Bc = predConc_kint(:,33);
            pS5AcS5Bc = predConc_kint(:,34);
            pS5AnpS5An = predConc_kint(:,35);
            pS5AnpS5Bn = predConc_kint(:,36);
            pS5BnpS5Bn = predConc_kint(:,37);
            pS5An = predConc_kint(:,38);
            pS5Bn = predConc_kint(:,39);
            PPNpS5An = predConc_kint(:,40);
            PPNpS5Bn = predConc_kint(:,41);
            S5An = predConc_kint(:,42);
            S5Bn = predConc_kint(:,43);
            PPNpS5AnpS5An = predConc_kint(:,44);
            PPNpS5AnpS5Bn = predConc_kint(:,45);
            PPNpS5BnpS5Bn = predConc_kint(:,46);
            pS5AnS5An = predConc_kint(:,47);
            pS5BnS5Bn = predConc_kint(:,48);
            S5AnpS5Bn = predConc_kint(:,49);
            pS5AnS5Bn = predConc_kint(:,50);
            mRNAn = predConc_kint(:,51);
            mRNAc = predConc_kint(:,52);
            SOCS1 = predConc_kint(:,53);
            SOCS1PRLRJ2a = predConc_kint(:,54);
            iSOCS1PRLRJ2a = predConc_kint(:,55);
            PRLRJ2aS5AcSHP2 = predConc_kint(:,56);
            iPRLRJ2aS5AcSHP2 = predConc_kint(:,57);
            PRLRJ2aS5BcSHP2 = predConc_kint(:,58);
            iPRLRJ2aS5BcSHP2 = predConc_kint(:,59);
            SOCS1PRLRJ2aSHP2 = predConc_kint(:,60);
            iSOCS1PRLRJ2aSHP2 = predConc_kint(:,61);
            mRn = predConc_kint(:,62);
            mRc = predConc_kint(:,63);
            Rc = predConc_kint(:,64);
            mBCLn = predConc_kint(:,65);
            mBCLc = predConc_kint(:,66);
            BCL = predConc_kint(:,67);

            % Store each predConc entry at the correct index of the cell
            outputCell_kint{i,j,k} = predConc_kint(:,:);
            
            % Calculate and store end value ratio of PRLRJ2a/iPRLRJ2a
            ratio_kint = PRLRJ2a./iPRLRJ2a;
            ratio_endValue_kint = ratio_kint(end);
            Z_kint_PRLRJ2aRatio_endValues(i,j,k) = ratio_endValue_kint;
        end
    end
end
toc

Z_kint_1 = Z_kint_PRLRJ2aRatio_endValues(:,:,1);
Z_kint_2 = Z_kint_PRLRJ2aRatio_endValues(:,:,2);
Z_kint_3 = Z_kint_PRLRJ2aRatio_endValues(:,:,3);
Z_kint_4 = Z_kint_PRLRJ2aRatio_endValues(:,:,4);
Z_kint_5 = Z_kint_PRLRJ2aRatio_endValues(:,:,5);

% Simulate varying krec cycling rate values %
tic
for k = 1:numInitialConditions
    for i = 1:length(krecfreeMesh)
        for j = 1:length(krecboundMesh)
            params = baselineParams;
            params(62) = kintfreeBaseline;
            params(63) = krecfreeMesh(i);
            params(64) = kintboundBaseline;
            params(65) = krecboundMesh(j);
            initvalues = initialConditions(k,:);

            [~, predConc_krec] = ode15s(@core_endosomalCompartment,predTime,initvalues,options,params);
            
            PRL = predConc_krec(:,1);
            RJ = predConc_krec(:,2);
            iRJ = predConc_krec(:,3);
            S5Ac = predConc_krec(:,4);
            S5Bc = predConc_krec(:,5);
            SHP2 = predConc_krec(:,6);
            PPX = predConc_krec(:,7);
            PPN = predConc_krec(:,8);
            PRLRJ = predConc_krec(:,9);
            iPRLRJ = predConc_krec(:,10);
            PRLRJ2 = predConc_krec(:,11);
            iPRLRJ2 = predConc_krec(:,12);
            PRLRJ2a = predConc_krec(:,13);
            iPRLRJ2a = predConc_krec(:,14);
            PRLRJ2aS5Ac = predConc_krec(:,15);
            iPRLRJ2aS5Ac = predConc_krec(:,16);
            PRLRJ2aS5Bc = predConc_krec(:,17);
            iPRLRJ2aS5Bc = predConc_krec(:,18);
            pS5Ac = predConc_krec(:,19);
            pS5Bc = predConc_krec(:,20);
            pS5AcpS5Ac = predConc_krec(:,21);
            pS5AcpS5Bc = predConc_krec(:,22);
            pS5BcpS5Bc = predConc_krec(:,23);
            PRLRJ2aSHP2 = predConc_krec(:,24);
            iPRLRJ2aSHP2 = predConc_krec(:,25);
            PPXpS5Ac = predConc_krec(:,26);
            PPXpS5Bc = predConc_krec(:,27);
            PPXpS5AcpS5Ac = predConc_krec(:,28);
            PPXpS5BcpS5Bc = predConc_krec(:,29);
            PPXpS5AcpS5Bc = predConc_krec(:,30);
            pS5AcS5Ac = predConc_krec(:,31);
            pS5BcS5Bc = predConc_krec(:,32);
            S5AcpS5Bc = predConc_krec(:,33);
            pS5AcS5Bc = predConc_krec(:,34);
            pS5AnpS5An = predConc_krec(:,35);
            pS5AnpS5Bn = predConc_krec(:,36);
            pS5BnpS5Bn = predConc_krec(:,37);
            pS5An = predConc_krec(:,38);
            pS5Bn = predConc_krec(:,39);
            PPNpS5An = predConc_krec(:,40);
            PPNpS5Bn = predConc_krec(:,41);
            S5An = predConc_krec(:,42);
            S5Bn = predConc_krec(:,43);
            PPNpS5AnpS5An = predConc_krec(:,44);
            PPNpS5AnpS5Bn = predConc_krec(:,45);
            PPNpS5BnpS5Bn = predConc_krec(:,46);
            pS5AnS5An = predConc_krec(:,47);
            pS5BnS5Bn = predConc_krec(:,48);
            S5AnpS5Bn = predConc_krec(:,49);
            pS5AnS5Bn = predConc_krec(:,50);
            mRNAn = predConc_krec(:,51);
            mRNAc = predConc_krec(:,52);
            SOCS1 = predConc_krec(:,53);
            SOCS1PRLRJ2a = predConc_krec(:,54);
            iSOCS1PRLRJ2a = predConc_krec(:,55);
            PRLRJ2aS5AcSHP2 = predConc_krec(:,56);
            iPRLRJ2aS5AcSHP2 = predConc_krec(:,57);
            PRLRJ2aS5BcSHP2 = predConc_krec(:,58);
            iPRLRJ2aS5BcSHP2 = predConc_krec(:,59);
            SOCS1PRLRJ2aSHP2 = predConc_krec(:,60);
            iSOCS1PRLRJ2aSHP2 = predConc_krec(:,61);
            mRn = predConc_krec(:,62);
            mRc = predConc_krec(:,63);
            Rc = predConc_krec(:,64);
            mBCLn = predConc_krec(:,65);
            mBCLc = predConc_krec(:,66);
            BCL = predConc_krec(:,67);

            % Store each predConc entry at the correct index of the cell
            outputCell_krec{i,j,k} = predConc_krec(:,:);
            
            % Calculate and store end value ratio of PRLRJ2a/iPRLRJ2a
            ratio_krec = PRLRJ2a./iPRLRJ2a;
            ratio_endValue_krec = ratio_krec(end);
            Z_krec_PRLRJ2aRatio_endValues(i,j,k) = ratio_endValue_krec;
        end
    end
end
toc

Z_krec_1 = Z_krec_PRLRJ2aRatio_endValues(:,:,1);
Z_krec_2 = Z_krec_PRLRJ2aRatio_endValues(:,:,2);
Z_krec_3 = Z_krec_PRLRJ2aRatio_endValues(:,:,3);
Z_krec_4 = Z_krec_PRLRJ2aRatio_endValues(:,:,4);
Z_krec_5 = Z_krec_PRLRJ2aRatio_endValues(:,:,5);

% Simulate varying kintfree,krecbound cycling rate values %
tic
for k = 1:numInitialConditions
    for i = 1:length(kintfreeMesh)
        for j = 1:length(krecboundMesh)
            params = baselineParams;
            params(62) = kintfreeMesh(i);
            params(63) = krecfreeBaseline;
            params(64) = kintboundBaseline;
            params(65) = krecboundMesh(j);
            initvalues = initialConditions(k,:);

            [~, predConc_kintfree_krecbound] = ode15s(@core_endosomalCompartment,predTime,initvalues,options,params);

            PRL = predConc_kintfree_krecbound(:,1);
            RJ = predConc_kintfree_krecbound(:,2);
            iRJ = predConc_kintfree_krecbound(:,3);
            S5Ac = predConc_kintfree_krecbound(:,4);
            S5Bc = predConc_kintfree_krecbound(:,5);
            SHP2 = predConc_kintfree_krecbound(:,6);
            PPX = predConc_kintfree_krecbound(:,7);
            PPN = predConc_kintfree_krecbound(:,8);
            PRLRJ = predConc_kintfree_krecbound(:,9);
            iPRLRJ = predConc_kintfree_krecbound(:,10);
            PRLRJ2 = predConc_kintfree_krecbound(:,11);
            iPRLRJ2 = predConc_kintfree_krecbound(:,12);
            PRLRJ2a = predConc_kintfree_krecbound(:,13);
            iPRLRJ2a = predConc_kintfree_krecbound(:,14);
            PRLRJ2aS5Ac = predConc_kintfree_krecbound(:,15);
            iPRLRJ2aS5Ac = predConc_kintfree_krecbound(:,16);
            PRLRJ2aS5Bc = predConc_kintfree_krecbound(:,17);
            iPRLRJ2aS5Bc = predConc_kintfree_krecbound(:,18);
            pS5Ac = predConc_kintfree_krecbound(:,19);
            pS5Bc = predConc_kintfree_krecbound(:,20);
            pS5AcpS5Ac = predConc_kintfree_krecbound(:,21);
            pS5AcpS5Bc = predConc_kintfree_krecbound(:,22);
            pS5BcpS5Bc = predConc_kintfree_krecbound(:,23);
            PRLRJ2aSHP2 = predConc_kintfree_krecbound(:,24);
            iPRLRJ2aSHP2 = predConc_kintfree_krecbound(:,25);
            PPXpS5Ac = predConc_kintfree_krecbound(:,26);
            PPXpS5Bc = predConc_kintfree_krecbound(:,27);
            PPXpS5AcpS5Ac = predConc_kintfree_krecbound(:,28);
            PPXpS5BcpS5Bc = predConc_kintfree_krecbound(:,29);
            PPXpS5AcpS5Bc = predConc_kintfree_krecbound(:,30);
            pS5AcS5Ac = predConc_kintfree_krecbound(:,31);
            pS5BcS5Bc = predConc_kintfree_krecbound(:,32);
            S5AcpS5Bc = predConc_kintfree_krecbound(:,33);
            pS5AcS5Bc = predConc_kintfree_krecbound(:,34);
            pS5AnpS5An = predConc_kintfree_krecbound(:,35);
            pS5AnpS5Bn = predConc_kintfree_krecbound(:,36);
            pS5BnpS5Bn = predConc_kintfree_krecbound(:,37);
            pS5An = predConc_kintfree_krecbound(:,38);
            pS5Bn = predConc_kintfree_krecbound(:,39);
            PPNpS5An = predConc_kintfree_krecbound(:,40);
            PPNpS5Bn = predConc_kintfree_krecbound(:,41);
            S5An = predConc_kintfree_krecbound(:,42);
            S5Bn = predConc_kintfree_krecbound(:,43);
            PPNpS5AnpS5An = predConc_kintfree_krecbound(:,44);
            PPNpS5AnpS5Bn = predConc_kintfree_krecbound(:,45);
            PPNpS5BnpS5Bn = predConc_kintfree_krecbound(:,46);
            pS5AnS5An = predConc_kintfree_krecbound(:,47);
            pS5BnS5Bn = predConc_kintfree_krecbound(:,48);
            S5AnpS5Bn = predConc_kintfree_krecbound(:,49);
            pS5AnS5Bn = predConc_kintfree_krecbound(:,50);
            mRNAn = predConc_kintfree_krecbound(:,51);
            mRNAc = predConc_kintfree_krecbound(:,52);
            SOCS1 = predConc_kintfree_krecbound(:,53);
            SOCS1PRLRJ2a = predConc_kintfree_krecbound(:,54);
            iSOCS1PRLRJ2a = predConc_kintfree_krecbound(:,55);
            PRLRJ2aS5AcSHP2 = predConc_kintfree_krecbound(:,56);
            iPRLRJ2aS5AcSHP2 = predConc_kintfree_krecbound(:,57);
            PRLRJ2aS5BcSHP2 = predConc_kintfree_krecbound(:,58);
            iPRLRJ2aS5BcSHP2 = predConc_kintfree_krecbound(:,59);
            SOCS1PRLRJ2aSHP2 = predConc_kintfree_krecbound(:,60);
            iSOCS1PRLRJ2aSHP2 = predConc_kintfree_krecbound(:,61);
            mRn = predConc_kintfree_krecbound(:,62);
            mRc = predConc_kintfree_krecbound(:,63);
            Rc = predConc_kintfree_krecbound(:,64);
            mBCLn = predConc_kintfree_krecbound(:,65);
            mBCLc = predConc_kintfree_krecbound(:,66);
            BCL = predConc_kintfree_krecbound(:,67);

            % Store each predConc entry at the correct index of the cell
            outputCell_kintfree_krecbound{i,j,k} = predConc_kintfree_krecbound(:,:);
            
            % Calculate and store end value ratio of PRLRJ2a/iPRLRJ2a
            ratio_kintfree_krecbound = PRLRJ2a./iPRLRJ2a;
            ratio_endValue_kintfree_krecbound = ratio_kintfree_krecbound(end);
            Z_kintfree_krecbound_PRLRJ2aRatio_endValues(i,j,k) = ratio_endValue_kintfree_krecbound;
        end
    end
end
toc

Z_kintfree_krecbound_1 = Z_kintfree_krecbound_PRLRJ2aRatio_endValues(:,:,1);
Z_kintfree_krecbound_2 = Z_kintfree_krecbound_PRLRJ2aRatio_endValues(:,:,2);
Z_kintfree_krecbound_3 = Z_kintfree_krecbound_PRLRJ2aRatio_endValues(:,:,3);
Z_kintfree_krecbound_4 = Z_kintfree_krecbound_PRLRJ2aRatio_endValues(:,:,4);
Z_kintfree_krecbound_5 = Z_kintfree_krecbound_PRLRJ2aRatio_endValues(:,:,5);

% Simulate varying krecfree,kintbound cycling rate values %
tic
for k = 1:numInitialConditions
    for i = 1:length(krecfreeMesh)
        for j = 1:length(kintboundMesh)
            params = baselineParams;
            params(62) = kintfreeBaseline;
            params(63) = krecfreeMesh(i);
            params(64) = kintboundMesh(j);
            params(65) = krecboundBaseline;
            initvalues = initialConditions(k,:);

            [~, predConc_krecfree_kintbound] = ode15s(@core_endosomalCompartment,predTime,initvalues,options,params);

            PRL = predConc_krecfree_kintbound(:,1);
            RJ = predConc_krecfree_kintbound(:,2);
            iRJ = predConc_krecfree_kintbound(:,3);
            S5Ac = predConc_krecfree_kintbound(:,4);
            S5Bc = predConc_krecfree_kintbound(:,5);
            SHP2 = predConc_krecfree_kintbound(:,6);
            PPX = predConc_krecfree_kintbound(:,7);
            PPN = predConc_krecfree_kintbound(:,8);
            PRLRJ = predConc_krecfree_kintbound(:,9);
            iPRLRJ = predConc_krecfree_kintbound(:,10);
            PRLRJ2 = predConc_krecfree_kintbound(:,11);
            iPRLRJ2 = predConc_krecfree_kintbound(:,12);
            PRLRJ2a = predConc_krecfree_kintbound(:,13);
            iPRLRJ2a = predConc_krecfree_kintbound(:,14);
            PRLRJ2aS5Ac = predConc_krecfree_kintbound(:,15);
            iPRLRJ2aS5Ac = predConc_krecfree_kintbound(:,16);
            PRLRJ2aS5Bc = predConc_krecfree_kintbound(:,17);
            iPRLRJ2aS5Bc = predConc_krecfree_kintbound(:,18);
            pS5Ac = predConc_krecfree_kintbound(:,19);
            pS5Bc = predConc_krecfree_kintbound(:,20);
            pS5AcpS5Ac = predConc_krecfree_kintbound(:,21);
            pS5AcpS5Bc = predConc_krecfree_kintbound(:,22);
            pS5BcpS5Bc = predConc_krecfree_kintbound(:,23);
            PRLRJ2aSHP2 = predConc_krecfree_kintbound(:,24);
            iPRLRJ2aSHP2 = predConc_krecfree_kintbound(:,25);
            PPXpS5Ac = predConc_krecfree_kintbound(:,26);
            PPXpS5Bc = predConc_krecfree_kintbound(:,27);
            PPXpS5AcpS5Ac = predConc_krecfree_kintbound(:,28);
            PPXpS5BcpS5Bc = predConc_krecfree_kintbound(:,29);
            PPXpS5AcpS5Bc = predConc_krecfree_kintbound(:,30);
            pS5AcS5Ac = predConc_krecfree_kintbound(:,31);
            pS5BcS5Bc = predConc_krecfree_kintbound(:,32);
            S5AcpS5Bc = predConc_krecfree_kintbound(:,33);
            pS5AcS5Bc = predConc_krecfree_kintbound(:,34);
            pS5AnpS5An = predConc_krecfree_kintbound(:,35);
            pS5AnpS5Bn = predConc_krecfree_kintbound(:,36);
            pS5BnpS5Bn = predConc_krecfree_kintbound(:,37);
            pS5An = predConc_krecfree_kintbound(:,38);
            pS5Bn = predConc_krecfree_kintbound(:,39);
            PPNpS5An = predConc_krecfree_kintbound(:,40);
            PPNpS5Bn = predConc_krecfree_kintbound(:,41);
            S5An = predConc_krecfree_kintbound(:,42);
            S5Bn = predConc_krecfree_kintbound(:,43);
            PPNpS5AnpS5An = predConc_krecfree_kintbound(:,44);
            PPNpS5AnpS5Bn = predConc_krecfree_kintbound(:,45);
            PPNpS5BnpS5Bn = predConc_krecfree_kintbound(:,46);
            pS5AnS5An = predConc_krecfree_kintbound(:,47);
            pS5BnS5Bn = predConc_krecfree_kintbound(:,48);
            S5AnpS5Bn = predConc_krecfree_kintbound(:,49);
            pS5AnS5Bn = predConc_krecfree_kintbound(:,50);
            mRNAn = predConc_krecfree_kintbound(:,51);
            mRNAc = predConc_krecfree_kintbound(:,52);
            SOCS1 = predConc_krecfree_kintbound(:,53);
            SOCS1PRLRJ2a = predConc_krecfree_kintbound(:,54);
            iSOCS1PRLRJ2a = predConc_krecfree_kintbound(:,55);
            PRLRJ2aS5AcSHP2 = predConc_krecfree_kintbound(:,56);
            iPRLRJ2aS5AcSHP2 = predConc_krecfree_kintbound(:,57);
            PRLRJ2aS5BcSHP2 = predConc_krecfree_kintbound(:,58);
            iPRLRJ2aS5BcSHP2 = predConc_krecfree_kintbound(:,59);
            SOCS1PRLRJ2aSHP2 = predConc_krecfree_kintbound(:,60);
            iSOCS1PRLRJ2aSHP2 = predConc_krecfree_kintbound(:,61);
            mRn = predConc_krecfree_kintbound(:,62);
            mRc = predConc_krecfree_kintbound(:,63);
            Rc = predConc_krecfree_kintbound(:,64);
            mBCLn = predConc_krecfree_kintbound(:,65);
            mBCLc = predConc_krecfree_kintbound(:,66);
            BCL = predConc_krecfree_kintbound(:,67);

            % Store each predConc entry at the correct index of the cell
            outputCell_krecfree_kintbound{i,j,k} = predConc_krecfree_kintbound(:,:);

            % Calculate and store end value ratio of PRLRJ2a/iPRLRJ2a
            ratio_krecfree_kintbound = PRLRJ2a./iPRLRJ2a;
            ratio_endValue_krecfree_kintbound = ratio_krecfree_kintbound(end);
            Z_krecfree_kintbound_PRLRJ2aRatio_endValues(i,j,k) = ratio_endValue_krecfree_kintbound;
        end
    end
end
toc

Z_krecfree_kintbound_1 = Z_krecfree_kintbound_PRLRJ2aRatio_endValues(:,:,1);
Z_krecfree_kintbound_2 = Z_krecfree_kintbound_PRLRJ2aRatio_endValues(:,:,2);
Z_krecfree_kintbound_3 = Z_krecfree_kintbound_PRLRJ2aRatio_endValues(:,:,3);
Z_krecfree_kintbound_4 = Z_krecfree_kintbound_PRLRJ2aRatio_endValues(:,:,4);
Z_krecfree_kintbound_5 = Z_krecfree_kintbound_PRLRJ2aRatio_endValues(:,:,5);

%% Heat Map Plots %%

fontSize = 18;
N = 250;
numInitialConditions = 5;

% Find the minimum and maximum of each dataset
bottom = zeros(1,numInitialConditions);
top = zeros(1,numInitialConditions);
Z_kfree = zeros(N,N,numInitialConditions);
Z_kbound = zeros(N,N,numInitialConditions);
Z_kint = zeros(N,N,numInitialConditions);
Z_krec = zeros(N,N,numInitialConditions);
Z_kintfree_krecbound = zeros(N,N,numInitialConditions);
Z_krecfree_kintbound = zeros(N,N,numInitialConditions);

for i = 1:numInitialConditions
    Z_kfree(:,:,i) = eval(strcat('Z_kfree_', int2str(i)));
    Z_kbound(:,:,i) = eval(strcat('Z_kbound_', int2str(i)));
    Z_kint(:,:,i) = eval(strcat('Z_kint_', int2str(i)));
    Z_krec(:,:,i) = eval(strcat('Z_krec_', int2str(i)));
    Z_kintfree_krecbound(:,:,i) = eval(strcat('Z_kintfree_krecbound_', int2str(i)));
    Z_krecfree_kintbound(:,:,i) = eval(strcat('Z_krecfree_kintbound_', int2str(i)));
    % disp(strcat('Z_kfree_', int2str(i)))
end

% min can only take the minimum value from two arrays, min(A,B)

for i = 1:numInitialConditions
    bottom(i) = min([min(Z_kfree(:,:,i),[],"all"),...
        min(Z_kbound(:,:,i),[],"all"),...
        min(Z_kint(:,:,i),[],"all"),...
        min(Z_krec(:,:,i),[],"all"),...
        min(Z_kintfree_krecbound(:,:,i),[],"all"),...
        min(Z_krecfree_kintbound(:,:,i),[],"all")]); % will be the minimum value of the colorbar
    top(i) = max([max(Z_kfree(:,:,i),[],"all"),...
        max(Z_kbound(:,:,i),[],"all"),...
        max(Z_kint(:,:,i),[],"all"),...
        max(Z_krec(:,:,i),[],"all"),...
        max(Z_kintfree_krecbound(:,:,i),[],"all"),...
        max(Z_krecfree_kintbound(:,:,i),[],"all")]); % will be the maximum value of the colorbar
end

% surface plots with contours for PRLRJ2a/iPRLRJ2a at 6 h %
clim_min = min(bottom);
clim_max = max(top);

%% 1. Plot Z_kfree ICs 1-5 %%

% IC1 %
figure
surf(krecfreeMesh,kintfreeMesh,Z_kfree_1(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]);  
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2 0.25])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecfreeMesh,kintfreeMesh,Z_kfree_1(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecfreeBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecfreeMin krecfreeMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecfree")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kfree: t_{0} = internal 100:0 surface');
fileName = 'kfree_allEndo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC2 %
figure
surf(krecfreeMesh,kintfreeMesh,Z_kfree_2(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]);  
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2 0.25])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecfreeMesh,kintfreeMesh,Z_kfree_2(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecfreeBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecfreeMin krecfreeMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecfree")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kfree: t_{0} = internal 75:25 surface');
fileName = 'kfree_75Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC3 %
figure
surf(krecfreeMesh,kintfreeMesh,Z_kfree_3(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2 0.25])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecfreeMesh,kintfreeMesh,Z_kfree_3(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecfreeBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecfreeMin krecfreeMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecfree")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kfree: t_{0} = internal 50:50 surface');
fileName = 'kfree_50Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC4 %
figure
surf(krecfreeMesh,kintfreeMesh,Z_kfree_4(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2 0.25])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecfreeMesh,kintfreeMesh,Z_kfree_4(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecfreeBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecfreeMin krecfreeMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecfree")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kfree: t_{0} = internal 25:75 surface');
fileName = 'kfree_25Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC5 %
figure
surf(krecfreeMesh,kintfreeMesh,Z_kfree_5(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2 0.25])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecfreeMesh,kintfreeMesh,Z_kfree_5(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecfreeBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecfreeMin krecfreeMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecfree")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kfree: t_{0} = internal 0:100 surface');
fileName = 'kfree_0Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

%% 2. Plot Z_kbound ICs 1-5 %%

% IC1 %
figure
surf(krecboundMesh,kintboundMesh,Z_kbound_1(:,:,1)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2])
contour(krecboundMesh,kintboundMesh,Z_kbound_1(:,:,1)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintboundBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintboundMin kintboundMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintbound")
fontsize(14,'points')
axis square
title('vary kbound: t_{0} = internal 100:0 surface');
fileName = 'kbound_allEndo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');


% IC2 %
figure
surf(krecboundMesh,kintboundMesh,Z_kbound_2(:,:,1)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2])
contour(krecboundMesh,kintboundMesh,Z_kbound_2(:,:,1)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintboundBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintboundMin kintboundMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintbound")
fontsize(14,'points')
axis square
title('vary kbound: t_{0} = internal 75:25 surface');
fileName = 'kbound_75Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');


% IC3 %
figure
surf(krecboundMesh,kintboundMesh,Z_kbound_3(:,:,1)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2])
contour(krecboundMesh,kintboundMesh,Z_kbound_3(:,:,1)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintboundBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintboundMin kintboundMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintbound")
fontsize(14,'points')
axis square
title('vary kbound: t_{0} = internal 50:50 surface');
fileName = 'kbound_50Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');


% IC4 %
figure
surf(krecboundMesh,kintboundMesh,Z_kbound_4(:,:,1)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2])
contour(krecboundMesh,kintboundMesh,Z_kbound_4(:,:,1)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintboundBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintboundMin kintboundMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintbound")
fontsize(14,'points')
axis square
title('vary kbound: t_{0} = internal 25:75 surface');
fileName = 'kbound_25Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC5 %
figure
surf(krecboundMesh,kintboundMesh,Z_kbound_5(:,:,1)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2])
contour(krecboundMesh,kintboundMesh,Z_kbound_5(:,:,1)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintboundBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintboundMin kintboundMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintbound")
fontsize(14,'points')
axis square
title('vary kbound: t_{0} = internal 0:100 surface');
fileName = 'kbound_0Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

%% 3. Plot Z_kint ICs 1-5 %%

% IC1 %
figure
surf(kintboundMesh,kintfreeMesh,Z_kint_1(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(kintboundMesh,kintfreeMesh,Z_kint_1(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kint: t_{0} = internal 100:0 surface');
fileName = 'kint_allEndo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC2 %
figure
surf(kintboundMesh,kintfreeMesh,Z_kint_2(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(kintboundMesh,kintfreeMesh,Z_kint_2(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kint: t_{0} = internal 75:25 surface');
fileName = 'kint_75Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC3 %
figure
surf(kintboundMesh,kintfreeMesh,Z_kint_3(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(kintboundMesh,kintfreeMesh,Z_kint_3(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kint: t_{0} = internal 50:50 surface');
fileName = 'kint_50Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC4 %
figure
surf(kintboundMesh,kintfreeMesh,Z_kint_4(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(kintboundMesh,kintfreeMesh,Z_kint_4(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kint: t_{0} = internal 25:75 surface');
fileName = 'kint_25Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC5 %
figure
surf(kintboundMesh,kintfreeMesh,Z_kint_5(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(kintboundMesh,kintfreeMesh,Z_kint_5(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kint: t_{0} = internal 0:100 surface');
fileName = 'kint_0Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

%% 4. Plot Z_krec ICs 1-5 %%

% IC1 %
figure
surf(krecboundMesh,krecfreeMesh,Z_krec_1(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(krecboundMesh,krecfreeMesh,Z_krec_1(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krec: t_{0} = internal 100:0 surface');
fileName = 'krec_allEndo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC2 %
figure
surf(krecboundMesh,krecfreeMesh,Z_krec_2(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(krecboundMesh,krecfreeMesh,Z_krec_2(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krec: t_{0} = internal 75:25 surface');
fileName = 'krec_75Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC3 %
figure
surf(krecboundMesh,krecfreeMesh,Z_krec_3(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(krecboundMesh,krecfreeMesh,Z_krec_3(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krec: t_{0} = internal 50:50 surface');
fileName = 'krec_50Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC4 %
figure
surf(krecboundMesh,krecfreeMesh,Z_krec_4(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(krecboundMesh,krecfreeMesh,Z_krec_4(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krec: t_{0} = internal 25:75 surface');
fileName = 'krec_25Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC5 %
figure
surf(krecboundMesh,krecfreeMesh,Z_krec_5(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(krecboundMesh,krecfreeMesh,Z_krec_5(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krec: t_{0} = internal 0:100 surface');
fileName = 'krec_0Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

%% 5. Plot Z_kintfree_krecbound ICs 1-5 %%

% IC1 %
figure
surf(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_1(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_1(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kintf,krecb: t_{0} = internal 100:0 surface');
fileName = 'kintf_krecb_allEndo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC2 %
figure
surf(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_2(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_2(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kintf,krecb: t_{0} = internal 75:25 surface');
fileName = 'kintf_krecb_75Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC3 %
figure
surf(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_3(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_3(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kintf,krecb: t_{0} = internal 50:50 surface');
fileName = 'kintf_krecb_50Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC4 %
figure
surf(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_4(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_4(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kintf,krecb: t_{0} = internal 25:75 surface');
fileName = 'kintf_krecb_25Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC5 %
figure
surf(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_5(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.005 0.01 0.015])
% yticks([0 1e-3 2e-3 3e-3 4e-3 5e-3 6e-3 7e-3])
contour(krecboundMesh,kintfreeMesh,Z_kintfree_krecbound_5(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(krecboundBaseline,kintfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([krecboundMin krecboundMax]);
ylim([kintfreeMin kintfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("krecbound")
ylabel("kintfree")
fontsize(14,'points')
axis square
title('vary kintf,krecb: t_{0} = internal 0:100 surface');
fileName = 'kintf_krecb_0Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

%% 6. Plot Z_krecfree_kintbound ICs 1-5 %%

% IC1 %
figure
surf(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_1(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_1(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krecf,kintb: t_{0} = internal 100:0 surface');
fileName = 'krecf_kintb_allEndo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC2 %
figure
surf(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_2(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_2(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krecf,kintb: t_{0} = internal 75:25 surface');
fileName = 'krecf_kintb_75Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC3 %
figure
surf(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_3(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_3(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krecf,kintb: t_{0} = internal 50:50 surface');
fileName = 'krecf_kintb_50Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC4 %
figure
surf(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_4(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_4(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krecf,kintb: t_{0} = internal 25:75 surface');
fileName = 'krecf_kintb_25Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');

% IC5 %
figure
surf(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_5(:,:)','LineWidth',1,'EdgeColor','none')
view(2)
hold on
colormap spring
clim("manual") 
clim([clim_min clim_max]); 
colorbar
set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'ColorScale','log','TickLength',[0.025, 0.01]);
% xticks([0 0.05 0.1 0.15 0.2])
% yticks([0 0.05 0.1 0.15 0.2 0.25])
contour(kintboundMesh,krecfreeMesh,Z_krecfree_kintbound_5(:,:)','LineWidth',1,'color','k','ZLocation','zmax')

% plot baseline param values %
scatter3(kintboundBaseline,krecfreeBaseline,10,40,'pentagram','MarkerEdgeColor',[0.3255 0.4078 0.4706],...
                  'MarkerFaceColor',[0.3255 0.4078 0.4706],...
                  'LineWidth',1.5);
xlim([kintboundMin kintboundMax]);
ylim([krecfreeMin krecfreeMax]);
xscale('log');
yscale('log');
grid off
hold off
xlabel("kintbound")
ylabel("krecfree")
fontsize(14,'points')
axis square
title('vary krecf,kintb: t_{0} = internal 0:100 surface');
fileName = 'krecf_kintb_0Endo_RJ_2D_contours.pdf';
exportgraphics(gcf,fullfile(fileName),'ContentType','vector');