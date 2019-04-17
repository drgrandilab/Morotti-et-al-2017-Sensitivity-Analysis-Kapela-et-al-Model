%% Kapela et al - Smooth Muscle Cell model

% J Theor Biol. 2008 Jul 21;253(2):238-60
% A mathematical model of Ca2+ dynamics in rat mesenteric smooth muscle cell:
% agonist and NO stimulation.
% Kapela A, Bezerianos A, Tsoukias NM.
% Dept. of Biomedical Engineering, Florida International University, Miami, FL, USA.
% 
% Matlab implementation by Stefano Morotti, University of California, Davis, Ca, USA.
% Main file.

close all
clear all
clc

%% Initial conditions
load yfin_icaMbar_rest      % ICs normal glucose
%load yfin_icaMbar_rest_GLU % ICs high glucose

y0n = yfinal;

%% Input parameters
% INDEX:  1      2        3         4            5         6          7
%       {'none','V_hold','V_clamp','NE-NO_stim','Ko_ramp','RMP_Ca_ss','HG-ramp'};
prot_index = 7; % see SMC_eccODEfile.m for protocol properties
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

%% Run simulation
tic
tspan = [0 240*1e3]; % [ms]
options = odeset('RelTol',1e-5,'MaxStep',1);
[t,y] = ode15s(@SMC_eccODEfile,tspan,y0n,options,p);
yfinal = y(end,:)';
toc

%% Save final conditions
%save yfin_icaMbar_rest yfinal
%save yfin_icaMbar_rest_GLU yfinal

%% Rename output variables
tArray = tArray_store(1:tStep); NE = NE_store(1:tStep); NO = NO_store(1:tStep); 
K_o = K_o_store(1:tStep); I_stim = I_stim_store(1:tStep); I_CaL = I_CaL_store(1:tStep); 
I_KCa = I_KCa_store(1:tStep); I_K = I_K_store(1:tStep); I_PMCA = I_PMCA_store(1:tStep);
I_NaK = I_NaK_store(1:tStep); I_NCX = I_NCX_store(1:tStep); 
INa_NSC = INa_NSC_store(1:tStep); ICa_NSC = ICa_NSC_store(1:tStep);
IK_NSC = IK_NSC_store(1:tStep); I_rel = I_rel_store(1:tStep); 
I_CaL_H = I_CaL_H_store(1:tStep);

%% Plot results
figure(1), set(gcf,'color','w')
subplot(2,1,1), hold on, set(gca,'box','off','tickdir','out','fontsize',12)
plot(1/1000*t,y(:,27)), ylabel('Em (mV)')
subplot(2,1,2), hold on, set(gca,'box','off','tickdir','out','fontsize',12)
plot(1/1000*t,1e6*y(:,23)), ylabel('[Ca2+]i (nM)'), xlabel('Time (s)')