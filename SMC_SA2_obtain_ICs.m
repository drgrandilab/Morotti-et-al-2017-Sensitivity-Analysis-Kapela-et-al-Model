%% Sensitivity Analyisis - part 2/3
% 
% S. Morotti, M. Nieves-Cintrón, M.A. Nystoriak, M.F. Navedo, E. Grandi.
% Predominant contribution of L-type CaV1.2 channel stimulation to impaired
% intracellular calcium and cerebral artery vasoconstriction in diabetic 
% hyperglycemia. Channels. 2017. doi: 10.1080/19336950.2017.1293220.
% 
% Please, cite the above paper when using this code.

%% main obtain ICs
close all
clear all
clc

%% Initial conditions
load yfin_icaMbar_rest

y0n = yfinal;
N_state_vars = length(y0n);

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load SA_par_matrix_1000_s0p1 % sigma 0.1

[N_trials N_par]=size(all_parameters);

%% Input parameters
% INDEX:  1      2        3         4            5         6
%       {'none','V_hold','V_clamp','NE-NO_stim','Ko_ramp','RMP_Ca_ss','HG-ramp'};
prot_index = 1;
p_SA = ones(1,15);
p = [prot_index p_SA];

%% Define global variables
global tStep tArray_store
global NE_store NO_store K_o_store I_stim_store
global I_CaL_store I_KCa_store I_K_store I_PMCA_store I_NaK_store I_NCX_store
global INa_NSC_store ICa_NSC_store IK_NSC_store I_rel_store
global I_CaL_H_store

tStep = 1; tArray_store = zeros(1,1e6); NE_store = zeros(1,1e6); NO_store = zeros(1,1e6);
K_o_store = zeros(1,1e6); I_stim_store = zeros(1,1e6); I_CaL_store = zeros(1,1e6);
I_KCa_store = zeros(1,1e6); I_K_store = zeros(1,1e6); I_PMCA_store = zeros(1,1e6);
I_NaK_store = zeros(1,1e6); I_NCX_store = zeros(1,1e6); INa_NSC_store = zeros(1,1e6);
ICa_NSC_store = zeros(1,1e6); IK_NSC_store = zeros(1,1e6); I_rel_store = zeros(1,1e6);
I_CaL_H_store = zeros(1,1e6); 

%% Run cycle
duration = 480*1e3;
tspan = [0 duration];
options = odeset('RelTol',1e-5,'MaxStep',1);

all_ICs = zeros(N_trials,N_state_vars);

tic
parfor ii=1:N_trials,
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    
    p_SA = all_parameters(ii,:); % 15 parameters
    p = [prot_index p_SA];
    [t,y] = ode15s(@SMC_eccODEfile,tspan,y0n,options,p);
    all_ICs(ii,:) = y(end,:);
end
toc

all_ICs;
% columns: N state variables
% rows: N trials

%% Saving
%save SA_ICs_matrix_1000_s0p1 all_ICs % Control