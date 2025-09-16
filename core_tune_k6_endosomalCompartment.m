function dydt = core_tune_k6_endosomalCompartment(t,y,params )

% Lynne Cherchia
% 26 August 2025
% Adds endosomal PRLR signaling compartment to RM's base model. Addresses
% ALL terms affected by iPRLRJ species. 
% ** Uses RM'S degradation terms **
 
% ***Tunes parameter k6, to explore the effects of internalized PRLR on
% signaling strength***

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
k6_internal = params(66);


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
dydt(14,1) = k4*iPRLRJ2 - k5*iPRLRJ2a*S5Ac + k_5*iPRLRJ2aS5Ac + k6_internal*iPRLRJ2aS5Ac - k5*iPRLRJ2a*S5Bc + k_5*iPRLRJ2aS5Bc + k6_internal*iPRLRJ2aS5Bc - k9*iPRLRJ2a*SHP2 + k_9*iPRLRJ2aSHP2 - k21*iPRLRJ2a*SOCS1 + k_21*iSOCS1PRLRJ2a + k23*iSOCS1PRLRJ2a - kdeg*deg_ratio*iPRLRJ2a + kint_bound*PRLRJ2a - krec_bound*iPRLRJ2a;

dydt(15,1) = k5*PRLRJ2a*S5Ac - k_5*PRLRJ2aS5Ac - k6*PRLRJ2aS5Ac - k9*SHP2*PRLRJ2aS5Ac + k_9*PRLRJ2aS5AcSHP2 - kdeg*deg_ratio*PRLRJ2aS5Ac - kint_bound*PRLRJ2aS5Ac + krec_bound*iPRLRJ2aS5Ac;
dydt(16,1) = k5*iPRLRJ2a*S5Ac - k_5*iPRLRJ2aS5Ac - k6_internal*iPRLRJ2aS5Ac - k9*SHP2*iPRLRJ2aS5Ac + k_9*iPRLRJ2aS5AcSHP2 - kdeg*deg_ratio*iPRLRJ2aS5Ac + kint_bound*PRLRJ2aS5Ac - krec_bound*iPRLRJ2aS5Ac;

dydt(17,1) = k5*PRLRJ2a*S5Bc - k_5*PRLRJ2aS5Bc - k6*PRLRJ2aS5Bc - k9*SHP2*PRLRJ2aS5Bc + k_9*PRLRJ2aS5BcSHP2 - kdeg*deg_ratio*PRLRJ2aS5Bc - kint_bound*PRLRJ2aS5Bc + krec_bound*iPRLRJ2aS5Bc;
dydt(18,1) = k5*iPRLRJ2a*S5Bc - k_5*iPRLRJ2aS5Bc - k6_internal*iPRLRJ2aS5Bc - k9*SHP2*iPRLRJ2aS5Bc + k_9*iPRLRJ2aS5BcSHP2 - kdeg*deg_ratio*iPRLRJ2aS5Bc + kint_bound*PRLRJ2aS5Bc - krec_bound*iPRLRJ2aS5Bc;

dydt(19,1) = k6*PRLRJ2aS5Ac + k6_internal*iPRLRJ2aS5Ac - 2*k8A*pS5Ac*pS5Ac + 2*k_8A*pS5AcpS5Ac - k8AB*pS5Ac*pS5Bc + k_8AB*pS5AcpS5Bc - k11*PPX*pS5Ac + k_11*PPXpS5Ac - k13*pS5Ac*S5Ac + k_13*pS5AcS5Ac - k13*pS5Ac*S5Bc + k_13*pS5AcS5Bc;
dydt(20,1) = k6*PRLRJ2aS5Bc + k6_internal*iPRLRJ2aS5Bc - 2*k8B*pS5Bc*pS5Bc + 2*k_8B*pS5BcpS5Bc - k8AB*pS5Ac*pS5Bc + k_8AB*pS5AcpS5Bc - k11*PPX*pS5Bc + k_11*PPXpS5Bc - k13*pS5Bc*S5Bc + k_13*pS5BcS5Bc - k13*pS5Bc*S5Ac + k_13*S5AcpS5Bc;

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

