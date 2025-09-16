% Lynne Cherchia
% 26 August 2025
% Run meshgridSims & plot 6h time course for PRLRJ2a/iPRLRJ2a

% FOR USE ON HPC (high-performance computing). 
% Run internalization model using meshgrid rate values.

% kfree %

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

% Set options %
options = odeset('RelTol',1e-9,'AbsTol',1e-12,'NonNegative',[1:size(selected_initvals,2)]);

% Define storage for outputs %
timeCourse_kfree_PRLRJ2aRatio = cell(length(kintfreeMesh),length(krecfreeMesh),length(numInitialConditions));    % store all prediction outputs (1441x67 for each param value)

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

            % Calculate and store each timepoint for PRLRJ2a/iPRLRJ2a
            ratio_kfree = PRLRJ2a./iPRLRJ2a;
            timeCourse_kfree_PRLRJ2aRatio{i,j,k} = ratio_kfree;
            % % Save ratio values to .mat
            % save timeCourse_kfree_PRLRJ2aRatio.mat timeCourse_kfree_PRLRJ2aRatio
        end
    end
end

%% Plot time course for PRLRJ2a/iPRLRJ2a 0-6 h %%

predTime = [0:60:24*3600];
cm = colormap(bone(250));
fontSize = 18;
fileExtension = '.pdf';

fileNameVector = cell(1,numInitialConditions);
fileNameVector{1,1} = 'kfree_allEndo_RJ_receptorRatioTimeCourse';
fileNameVector{1,2} = 'kfree_75Endo_RJ_receptorRatioTimeCourse';
fileNameVector{1,3} = 'kfree_50Endo_RJ_receptorRatioTimeCourse';
fileNameVector{1,4} = 'kfree_25Endo_RJ_receptorRatioTimeCourse';
fileNameVector{1,5} = 'kfree_allSurface_RJ_receptorRatioTimeCourse';


for k = 1:numInitialConditions
    figure('Position', [10 10 900 450]); 
    set(gca,'YDir','normal',...
    'fontsize',fontSize,...
    'fontname','arial',...
    'linewidth',1,...
    'layer','top',...
    'clipping','on',...
    'box','on',...
    'xscale','linear',...
    'TickLength',[0.025, 0.01]);
    xlim([0 6]);
    ylim([0 2]);
    xlabel('Time (h)');
    ylabel('^{PRLRJ2a}/_{iPRLRJ2a} complex concentration');
    hold on
    for i = 1:10:250
        for j = 1:10:250
            plot(predTime(1:361)/3600,timeCourse_kfree_PRLRJ2aRatio{i,j,k}(1:361),'Color',cm(j,:),'LineWidth',0.3);
        end
    end
    fileName = char(strcat(fileNameVector(k), fileExtension));
    exportgraphics(gcf,fullfile(fileName),'ContentType','vector');
end