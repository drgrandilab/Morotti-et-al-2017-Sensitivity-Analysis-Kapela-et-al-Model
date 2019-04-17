function dydt = SMC_eccODEfile(t,y,p)

% J Theor Biol. 2008 Jul 21;253(2):238-60
% A mathematical model of Ca2+ dynamics in rat mesenteric smooth muscle cell:
% agonist and NO stimulation.
% Kapela A, Bezerianos A, Tsoukias NM.
% Dept. of Biomedical Engineering, Florida International University, Miami, FL, USA.
% 
% Matlab implementation by Stefano Morotti, University of California, Davis, Ca, USA.
% Ordinary differential equations file.

% To modify external stimulations, go to section 'Stimulation protocols'. 
% To enable agonist-induced oscillations, set k_leak = 5, otherwise k_leak = 1.
% To remove desensitization to NE, set k_pG = 0.
% To introduce steady-state HG-induced modifications, set GLU to 1.

% UPDATE (June 19, 2018): delta_G set to a constant value (see line 315).

ydot = zeros(size(y));

%% Input parameters
prot_index = p(1); % protocol
p_SA = p(2:end); % sensitivity analysis

% INDEX:  1      2        3         4            5         6          7
%       {'none','V_hold','V_clamp','NE-NO_stim','Ko_ramp','RMP_Ca_ss','HG-ramp'};
switch prot_index
    case 1,
        protocol = 'none'; % no stimulation
	case 2,
        protocol = 'V_hold'; % Voltage clamped at fixed value 
	case 3,
        protocol = 'V_clamp'; % Voltage step
    case 4,
        protocol = 'NE-NO_stim'; % fig 6 Kapela et al. paper
    case 5,
        protocol = 'Ko_ramp'; % fig 5cd Kapela et al. paper
    case 6,
        protocol = 'RMP_Ca_ss'; % Sequence of voltage steps (-50, -40, -30 mV)
    case 7
        protocol = 'HG-ramp'; % Gradual increase in [glucose] (180 s)
end

%% D-glucose flag
GLU = 0; % 0=10mM (normal), 1=20mM (high) 

%% State variables
V_cGMP = y(1);
cGMP = y(2);
RS_G = y(3);
RS_PG = y(4);
G = y(5);
IP3	= y(6);
PIP2 = y(7);
DAG = y(8);
d_L = y(9);
f_L = y(10);
p_f = y(11);
p_s = y(12);
p_K = y(13);
q_1 = y(14);
q_2	= y(15);
h_IP3 = y(16);
P_SOC = y(17);
R_01 = y(18);
R_10 = y(19);
R_11 = y(20);
Ca_u = y(21);
Ca_r = y(22);
Ca_i = y(23);
Na_i = y(24);
K_i	= y(25);
Cl_i = y(26);
V_m = y(27);

%% Model Parameters
R = 8341.0;         % [mJmol-1K-1] gas constant
temp = 293.0;       % [K] absolute temperature
F = 96487.0;		% [Cmol-1] Faraday's constant
RT_F = R*temp/F;
z_K = 1;	    	% K ion valence
z_Na = 1;			% Na ion valence
z_Ca = 2;			% Ca ion valence
z_Cl = -1;          % Cl ion valence
N_Av = 6.022e23;    % Avogadro's constant

FoRT = F/R/temp;
Qpow = (temp-310)/10;

% Extracellular concentrations
Ca_o = 2.0;         % [mM] extracellular calcium
Na_o = 140.0;  		% [mM] extracellular sodium
Cl_o = 129.0;       % [mM] extracellular chloride
K_o = 5;            % [mM] extracellular potassium                 (PROTOCOL)
		
% Intracellular concentrations
NE = 0;             % [mM] norepinephrine concentration            (PROTOCOL)
NO = 0;             % [mM] exogenous nitric oxide concentration    (PROTOCOL)
vol_i = 1;          % [pl] volume of intracellular space
		
% Membrane capacitance and area 
C_m = 25;           % [pF] membrane capacitance 
A_m = C_m/1e6;      % [cm2] membrane area	
				
% Voltage dependent calcium current
LTCC_factor = 1; 
P_CaL = 1.88e-5*p_SA(1)*LTCC_factor; % [cm/s] max. conductance of L-type Ca channel
% Markov model: see equations
		
% Delayed rectifier current
g_K = 1.35*p_SA(3);         % [nS] max. conductance 
V_1_2 = -11.0;      % [mV] half-activation voltage
k = 15.0;           % [mV] slope factorof the Boltzmann function

% Alpha-adrenoceptor-activated nonselective cation channel NSC
Pmin_NSC = 0.4344;  % min. open probability 
PNa_NSC = (5.11e-7)*p_SA(11); % [cm/s]
PCa_NSC = (5.11e-7)*4.54*p_SA(11); % [cm/s]
PK_NSC = (5.11e-7)*1.06*p_SA(11); % [cm/s]	
d_NSCmin = 0.0244;  % 
K_NSC = 3.0e-3;     % [mM] half activation of NSC by IP3
	
% KATP current
g_KATP = 0.067*p_SA(4);     % [nS] max. background K current conductance

% Inward rectifier current
G_K_i = 0;          % inward rectifier constant
n_K_i = 0.5;        % inward rectifier constant

% Calcium-activated potassium current
P_BKCa = 3.9e-13*p_SA(2);   % [cm3/s] single channel permeability Mistry & Garland 1998
N_BKCa = 6.6e6;     % [1/cm2] channel density
tau_pf = 0.84;      % [ms] fast activation time constant, mean open times from Benham et al 1986  
tau_ps = 35.9;      % [ms] slow activation time constant
dV_1_2_KCa_NO = 46.3; % [mV] max. shift of the activation curve by direct NO 
dV_1_2_KCa_cGMP = 76; % [mV] max. shift of the activation curve by cGMP

% Calcium-induced calcium release
k_leak = 1;         % 1 - default, 5 - oscillation (during NE or K_o stim.)
Kr1 = 2500.0;       % [mM-2ms-1] activation rate constant
Kr2 = 1.05;         % [mM-1ms-1] inactivation rate constant
K_r1 = 0.0076;      % [ms-1] unbinding rate constant from activation
K_r2 = 0.084;       % [ms-1] unbinding rate constant from inactivation
I_relbar = 1*p_SA(14); % scaling factor
I_upbar = 3.34*(k_leak+1)*p_SA(15); % [pA] max. SR uptaking current
K_mup =  0.001;     % 0.08 % [mM] Michaelis constant of SR calcium pump
vol_u = 0.07;       % [pl] volume of uptake compartment
tau_tr = 1000.0;	% [ms] time constant of the internal diffusion
vol_r = 0.007;      % [pl] volume of release compartment
tau_rel= 0.0333;    % [ms]  
R_leak = 1.07e-5*(k_leak)*p_SA(13); % equal to R_10^2 during concentration clamp

% Ca buffering and cytosolic material balance
K_d = 2.6e-4;       % [mM] average binding affinity of CM for Ca
S_CMbar = 0.1;      % [mM] total concentration of CM sites for Ca
K_dB = 5.298e-4;    % [mM] average binding affinity of other buffers
B_Fbar = 0.1;       % [mM] total concentration of other buffer cites
vol_Ca = 0.7;       % [pl] intracellular volume available to free Ca
CSQNbar = 15;       % [mM] total calsequestrin concentration in release compartment
K_CSQN = 0.8;       % [mM] dissociation constant 
		
% Pump and exchanger currents
I_PMCAbar = 5.37*p_SA(9); % [pA] max. current through Ca pump
K_mPMCA = 170e-6;   % [mM] Michaelis constant of sarcolemmal Ca pump 
I_NaKbar = 2.3083*p_SA(7); % 1.91 % [pA/pF] max. Na/K current
K_mK = 1.6;         % [mM] K half saturation constant for Na/K pump
K_mNa = 22;         % [mM] Na half saturation constant for Na/K pump
Q_10_NaK = 1.87;    % the temperature coefficient Q10
gamma = 0.45;       %
g_NCX = 0.000487*p_SA(8); % [pA]
d_NCX = 0.0003;     %
L_cotr = 1.79e-8;   % [n(mol^2)/(J*s*(cm)^2)] = 1.79D-11 (mmol)^2/(J*s*(cm)^2) 

% Stretch-sensitive current: I_M
I_MCa = 0;          % [pA] Ca component
I_MNa = 0;          % [pA] Na component
I_MK = 0;           % [pA] K component
I_M = 0;            % [pA]

% Store operated non-selective cation channel
g_SOCCa = 0.0083*p_SA(10);   % [nS]
g_SOCNa = 0.0575*p_SA(10);   % [nS]
H_SOC = 1;          % Hill coefficient
K_SOC = 0.0001;     % [mM] half activation concentration
tau_SOC = 100;      % [ms] time constant for SOC activation 

% Chloride current
g_Cl = 0.23*p_SA(5); % [nS] max chloride conductance
g_NaKCl = 1*p_SA(6); % scaling factor

% IP3 receptor (Fink et al 1999)
I_IP3bar = 2880e-6*p_SA(12); % [1/ms]
K_IP3 = 0.12e-3;    % [mM] Bennett_JTheorBiol_2005
K_actIP3 = 0.17e-3; % [mM] Fink_BiophysJ_1999
K_inhIP3 = 0.1e-3;  % [mM]
k_onIP3 = 1.4;      % [1/(mM*ms)]

% Norepinephrine receptor (Bennett et al 2005)
R_TG = 2e4;         % total no. of receptors
K_1G = 0.01;        % [mM] unphosphorylated receptor dissociation constant
K_2G = 0.2;         % [mM] phosphorylated receptor dissociation constant
k_rG = 1.75e-7;     % [1/ms] receptor recycling rate
k_pG = 0.1e-3;      % [1/ms] receptor phosphorylation rate
k_eG = 6e-6;        % [1/ms] receptor endocytosis rate
ksi_G = 0.85;       % fraction of mobile receptors
G_TG = 1e5;         % total no. of G-protein molecules
k_degG = 1.25e-3;   % [1/ms] IP3 degradation rate
k_aG = 0.17e-3;     % [1/ms] G-protein activation rate
k_dG = 1.5e-3;      % [1/ms] G-protein deactivation rate
PIP2_T = 5e7;       % total number of PIP2 molecules
r_rG = 0.015e-3;    % [1/ms] PIP2 replenishment rate
K_cG = 0.4e-3;      % [mM] dissociation constant for Ca binding to PLC
alpha_G = 2.781e-8; % [1/ms] effective signal gain parameter
vol_IP3 = vol_i;    % [pL] volume for IP3 originally 0.5 pL
gamma_G = N_Av*vol_IP3*1e-15; % conversion of no. of PIP2 to mM concentration NO PAR

% cGMP formation
k1sGC = 2e3;        % [1/mM/ms]
k_1sGC = 15e-3;     % [1/ms]
k2sGC = 0.64e-5;    % [1/ms]
k_2sGC = 0.1e-6;    % [1/ms]
k3sGC = 4.2;        % [1/mM/ms]
kDsGC = 0.4e-3;     % [1/ms] 0.4e-3 corresponds to KmsGC = 93 nM
kDact_deactsGC = 0.1e-3; % [1/ms]
V_cGMPmax = 0.1*1.26e-6; % [mM/ms]
B5sGC = k2sGC/k3sGC; % 
A0sGC = ((k_1sGC+k2sGC)*kDsGC+k_1sGC*k_2sGC)/(k1sGC*k3sGC); % 
A1sGC = ((k1sGC+k3sGC)*kDsGC+(k2sGC+k_2sGC)*k1sGC)/(k1sGC*k3sGC); % 
kpde_cGMP = 0.0695e-3; % [1/ms]

% Stimulation current
I_stim = 0;         % [pA]                                          (PROTOCOL)

%% Simulation protocol
% Default vaues: NE = 0; NO = 0; K_o = 5; I_stim = 0.

switch protocol
    case 'none',
        I_stim = 0;
    case 'NE_stim', % Norepinephrine stimulation
        if t>=10e3, NE = 1e-3; end      % [mM] beginning of stimulation 
        if t>=70e3, NE = 0; end         % [mM] end of stimulation 
    case 'NO_stim', % Nitric oxide stimulation
        if t>=30e3, NO = 1e-3; end      % [mM] beginning of stimulation 
        if t>=70e3, NO = 0; end         % [mM] end of stimulation 
    case 'Ko_stim', % Extracellular potassium stimulation
        if t>=10e3, K_o = 30; end      % [mM] beginning of stimulation 
        if t>=70e3, K_o = 5; end       % [mM] end of stimulation
	case 'I_stim', % Current stimulation
        if t>=10e3, I_stim = 5; end    % [pA] current injection 
        if t>=40e3, I_stim = -5; end   % [pA] current injection
        if t>=70e3, I_stim = 0; end    % [pA] end of current injection
    case 'NE-NO_stim', % Fig. 6 stimulation
        if t>=10e3, NE = 3e-3; end      % [mM] beginning of NE stimulation 
        if t>=70e3, NO = 1e-3; end      % [mM] beginning of NO stimulation
    case 'Ko_ramp', % Fig. 5cd stimulation
        K_o = 5 + 5 * floor(t/200e3);
    case 'V_clamp', % Voltage step
        if t>=12.5e3 && t<17.5e3, 
            V_clamp = 5;
        else
            V_clamp = -59.333;
        end
        R_clamp = 0.01;
        I_stim = (V_clamp-V_m)/R_clamp*C_m;
	case 'V_hold', % Voltage clamp
        V_clamp = -80; % -59.333;
        R_clamp = 0.01;
        I_stim = (V_clamp-V_m)/R_clamp*C_m;
    case 'RMP_Ca_ss', % Voltage steps
        duration = 300e3;
        if t <= duration,
            I_stim = 0;
        else
            if t > duration && t <= 2*duration, V_clamp = -50; end
            if t > 2*duration && t <= 3*duration, V_clamp = -40; end
            if t > 3*duration, V_clamp = -30; end
            R_clamp = 0.01;
            I_stim = (V_clamp-V_m)/R_clamp*C_m;
        end
    case 'HG-ramp', % Gradual increase in [glucose] (180 s)
        I_stim = 0;
end

%% Glucose-dependent parameters
if GLU == 1,
    % HG-dependent steady-state modifications
    P_CaL = (1+1)*P_CaL; % ICa
    g_K = (1-0.63)*g_K; % IK
    P_BKCa = (1-0.58)*P_BKCa; % IK
end

if prot_index == 7,
    % HG-administration in 180 s (ramp) - fig 1 Morotti et al. paper
    P_CaL = (1+1*(t)/180e3)*1.88e-5*(t<180e3)+(1+1)*1.88e-5*(t>=180e3); % ICa
    g_K = (1-0.63*(t)/180e3)*1.35*(t<180e3)+(1-0.63)*1.35*(t>=180e3); % IK
    P_BKCa = (1-0.58*(t)/180e3)*3.9e-13*(t<180e3)+(1-0.58)*3.9e-13*(t>=180e3); % IK
end

%% cGMP formation
V_cGMPbar = V_cGMPmax*(B5sGC*NO+NO^2)/(A0sGC+A1sGC*NO+NO^2);
if ((V_cGMPbar-V_cGMP)>=0),
    tausGC = 1/(k3sGC*NO+kDact_deactsGC);
else
    tausGC = 1/(kDact_deactsGC+k_2sGC);
end
V_cGMP_dot = (V_cGMPbar-V_cGMP)/tausGC;
cGMP_dot = V_cGMP-kpde_cGMP*cGMP*cGMP/(1e-3+cGMP);
			
%% Norepinephrine receptor
RS_G_dot = (k_rG*ksi_G*R_TG-(k_rG+k_pG*NE/(K_1G+NE))*RS_G-k_rG*RS_PG);
RS_PG_dot = NE*(k_pG*RS_G/(K_1G+NE)-k_eG*RS_PG/(K_2G+NE));
rho_rG = NE*RS_G/(ksi_G*R_TG*(K_1G+NE));
%delta_G = k_dG*G/(k_aG*(G_TG-G));
delta_G = 0.001235; % from Lemon et al 2003 https://doi.org/10.1016/S0022-5193(03)00079-1
G_dot = k_aG*(delta_G+rho_rG)*(G_TG-G)-k_dG*G;
r_hG = alpha_G*Ca_i/(K_cG+Ca_i)*G;
IP3_dot = r_hG/gamma_G*PIP2-k_degG*IP3;
PIP2_dot = -(r_hG+r_rG)*PIP2-r_rG*gamma_G*IP3+r_rG*PIP2_T;
DAG_dot = r_hG/gamma_G*PIP2-k_degG*DAG;
	
%% Reversal potentials
E_Ca = RT_F/z_Ca * log(Ca_o/Ca_i); 
E_Na = RT_F/z_Na * log(Na_o/Na_i);
E_K = RT_F/z_K * log(K_o/K_i);
E_Cl = RT_F/z_Cl * log(Cl_o/Cl_i);

%% Voltage-dependent calcium current: I_CaL
% Original GHK formulation
tau_d_L = 2.5*exp(-((V_m+40)/30)^2)+1.15;
d_Lbar = 1/(1+exp(-(V_m)/8.3));
d_L_dot = (d_Lbar-d_L)/tau_d_L;
f_Lbar = 1/(1+exp((V_m+42.0)/9.1));
tau_f_L = 65*exp(-((V_m+35)/25)^2)+45;
f_L_dot = (f_Lbar-f_L)/tau_f_L;
if (V_m == 0),
    I_CaL = d_L*f_L*P_CaL*A_m*1e6*z_Ca*F*(Ca_i-Ca_o);
else
    I_CaL = d_L*f_L*P_CaL*A_m*1e6*V_m*((z_Ca*F)^2)/(R*temp)*(Ca_o-Ca_i*exp(V_m*z_Ca/(RT_F)))/(1-exp(V_m*z_Ca/(RT_F)));
end
I_Catot = I_CaL;

%% Calcium-activated potassium current: I_KCa (NO- & cGMP-dependent))
if (V_m == 0),
    i1_KCa = 1e6*P_BKCa*F*(K_i-K_o);
else
    i1_KCa = 1e6*P_BKCa*V_m*F/RT_F*(K_o-K_i*exp(V_m/RT_F))/(1-exp(V_m/RT_F));
end
R_NO = NO/(NO+200e-6);
R_cGMP = cGMP^2/(cGMP^2+(1.5e-3)^2);
V_1_2_KCa = -41.7*log10(Ca_i)-128.2-dV_1_2_KCa_NO*R_NO-dV_1_2_KCa_cGMP*R_cGMP; % NO-dependent activation from Yang 2005
p_obar = 1/(1+exp(-(V_m-V_1_2_KCa)/18.25));	
p_f_dot = (p_obar-p_f)/tau_pf;
p_s_dot = (p_obar-p_s)/tau_ps;
P_KCa = 0.17*p_f+0.83*p_s;
I_KCa = 1*A_m*N_BKCa*i1_KCa*P_KCa;

%% Delayed rectifier potassium current: I_K
p_Kbar = 1/(1+exp(-(V_m-V_1_2)/k));
tau_p_K = 61.49*exp(-0.0268*V_m);
p_K_dot = (p_Kbar-p_K)/tau_p_K;
q_bar = 1/(1+exp((V_m+40)/14));
q_1_dot = (q_bar-q_1)/371;
q_2_dot = (q_bar-q_2)/2884;
I_K = 1*g_K*p_K*(0.45*q_1+0.55*q_2)*(V_m-E_K);
	
%% Inward rectifier potassium current: I_K_i
g_maxK_i = G_K_i*(K_o^n_K_i); % not used, G_K_i set to 0
I_K_i = g_maxK_i*(V_m-E_K)/(1+exp((V_m-E_K)/28.89));
		
%% Background potassium current: I_KATP
I_KATP = g_KATP*(V_m-E_K);

%% IP3 receptor calcium current: I_IP3 (NE-dependent)
NE_factor = 1;
h_IP3_dot = k_onIP3*(K_inhIP3-(Ca_i+K_inhIP3)*h_IP3);
I_IP3 = I_IP3bar*((NE_factor*(IP3/(IP3+K_IP3))*Ca_i/(Ca_i+K_actIP3)*h_IP3)^3)*(Ca_u-Ca_i)*z_Ca*F*vol_Ca;

%% Sarcolemmal calcium pump current: I_PMCA
I_PMCA = I_PMCAbar*Ca_i/(Ca_i+K_mPMCA);

%% Sodium-potassium pump: I_NaK
I_NaK = (Q_10_NaK^((temp-309.15)/10))*C_m*I_NaKbar*((K_o^1.1)/((K_o^1.1)+(K_mK^1.1))*(Na_i^1.7)/((Na_i^1.7)+(K_mNa^1.7)))*(V_m+150)/(V_m+200);

%% Sodium-calcium exchanger current: I_NCX (cGMP-dependent)
Fi_F = exp(gamma*V_m*F/(R*temp));
Fi_R = exp((gamma-1)*V_m*F/(R*temp));
I_NCX = 1*(1+0.55*cGMP/(cGMP+(45e-3)))*g_NCX*((Na_i^3)*Ca_o*Fi_F-(Na_o^3)*Ca_i*Fi_R)/(1+d_NCX*((Na_o^3)*Ca_i+(Na_i^3)*Ca_o));

%% Calcium-dependent chloride current: I_Cl (cGMP-dependent)
alpha_Cl = (cGMP^3.3)/((cGMP^3.3)+(6.4e-3)^3.3);
P_Cl = (Ca_i^2)/(Ca_i^2+(365e-6)^2)*0.0132+(Ca_i^2)/(Ca_i^2+(400e-6*(1-alpha_Cl*0.9))^2)*alpha_Cl;
I_Cl = P_Cl*g_Cl*C_m*(V_m-E_Cl);

%% Sodium-potassium-chloride cotransport current: I_NaKCl_Cl (cGMP-dependent)
I_NaKCl_Cl = g_NaKCl*(1+7/2*cGMP/(cGMP+6.4e-3))*(-A_m*L_cotr*R*temp*z_Cl*F*log(Na_o/Na_i*K_o/K_i*(Cl_o/Cl_i)^2));

%% Alpha-adrenoceptor-activated nonselective cation channel NSC: I_NSC (NE-dependent)
Po_NSC = Pmin_NSC+(1-Pmin_NSC)/(1+exp(-(V_m-47.12)/24.24));
if (V_m == 0),
    INa_NSC = 1*(NE_factor*DAG/(DAG+K_NSC)+d_NSCmin)*Po_NSC*PNa_NSC*A_m*1e6*F*(Na_i-Na_o);
    ICa_NSC = 1*(0*NE_factor*DAG/(DAG+K_NSC)+d_NSCmin)*Po_NSC*PCa_NSC*A_m*1e6*z_Ca*F*(Ca_i-Ca_o);
    IK_NSC = 1*(NE_factor*DAG/(DAG+K_NSC)+d_NSCmin)*Po_NSC*PK_NSC*A_m*1e6*F*(K_i-K_o);
else
    INa_NSC = 1*(NE_factor*DAG/(DAG+K_NSC)+d_NSCmin)*Po_NSC*PNa_NSC*A_m*1e6*V_m*(F^2)/(R*temp)*(Na_o-Na_i*exp(V_m/RT_F))/(1-exp(V_m/RT_F));
    ICa_NSC = 1*(0*NE_factor*DAG/(DAG+K_NSC)+d_NSCmin)*Po_NSC*PCa_NSC*A_m*1e6*V_m*((z_Ca*F)^2)/(R*temp)*(Ca_o-Ca_i*exp(V_m*z_Ca/RT_F))/(1-exp(V_m*z_Ca/RT_F));
    IK_NSC = 1*(NE_factor*DAG/(DAG+K_NSC)+d_NSCmin)*Po_NSC*PK_NSC*A_m*1e6*V_m*(F^2)/(R*temp)*(K_o-K_i*exp(V_m/RT_F))/(1-exp(V_m/RT_F));
end
I_NSC = ICa_NSC+INa_NSC+IK_NSC;

%% Store operated non-selective cation current: I_SOC
P_SOCbar = 1/(1+(Ca_u/K_SOC)^(H_SOC));
P_SOC_dot = (P_SOCbar-P_SOC)/tau_SOC;
I_SOCCa = 1*P_SOC*g_SOCCa*(V_m-E_Ca);
I_SOCNa = 1*P_SOC*g_SOCNa*(V_m-E_Na);
I_SOC = I_SOCCa+I_SOCNa;

%% Calcium-induced calcium release (and leak): I_rel 
R_00 = 1-R_01-R_10-R_11;
R_10_dot = Kr1*(Ca_i^2)*R_00-(K_r1+Kr2*Ca_i)*R_10+K_r2*R_11;
R_11_dot = Kr2*Ca_i*R_10-(K_r1+K_r2)*R_11+Kr1*(Ca_i^2)*R_01;
R_01_dot = Kr2*Ca_i*R_00+K_r1*R_11-(K_r2+Kr1*(Ca_i^2))*R_01;
I_rel = (I_relbar*R_10^2+R_leak)*(Ca_r-Ca_i)*(2*F*vol_r)/tau_rel;

%% SERCA uptaking current: I_up
I_up = I_upbar*Ca_i/(Ca_i+K_mup);

%% Calcium diffusion between SERCA compartmtents: I_tr
I_tr = (Ca_u-Ca_r)*(2*F*vol_u)/tau_tr;

%% SERCA calcium concentrations
Ca_u_dot = (I_up-I_tr-I_IP3)/(2*F*vol_u);
Ca_r_dot = (I_tr-I_rel)/(2*F*vol_r)/(1+CSQNbar*K_CSQN/((K_CSQN+Ca_r)^2));

% %% Ca buffering and cytosolic material balance
% S_CM = S_CMbar*K_d/(K_d+Ca_i); % concentration of free calmodulin sites
% CaCM = S_CMbar-S_CM; % concentration of Ca-calmodulin complexes % not used

%% Intracellular concentrations
I_Catotm = I_SOCCa+I_Catot-2*I_NCX+I_PMCA+ICa_NSC+I_MCa;
Ca_i_dot = -(I_Catotm+I_up-I_rel-I_IP3)/(2*F*vol_Ca)/(1+S_CMbar*K_d/((K_d+Ca_i)^2)+B_Fbar*K_dB/((K_dB+Ca_i)^2));
 
I_Natotm = -0.5*I_NaKCl_Cl+I_SOCNa+3*I_NaK+3*I_NCX+INa_NSC+I_MNa;
Na_i_dot = -(I_Natotm)/(F*vol_i);

I_Ktotm = -0.5*I_NaKCl_Cl+I_K+I_KCa+I_K_i+IK_NSC+I_KATP-2*I_NaK+I_MK;
K_i_dot = -(I_Ktotm)/(F*vol_i);

I_Cltotm = I_NaKCl_Cl+I_Cl;
Cl_i_dot = -(I_Cltotm)/(z_Cl*F*vol_i);

%% Transmembrane potential
V_m_dot = -1/C_m*(-I_stim+I_Cl+I_SOC+I_Catot+I_K+I_KCa+I_K_i+I_M+I_NCX+I_NaK+I_PMCA+I_NSC+I_KATP);

%% Output
ydot(1) = V_cGMP_dot;
ydot(2) = cGMP_dot;
ydot(3) = RS_G_dot;
ydot(4) = RS_PG_dot;
ydot(5) = G_dot;
ydot(6) = IP3_dot;
ydot(7) = PIP2_dot;
ydot(8) = DAG_dot;
ydot(9) = d_L_dot;
ydot(10) = f_L_dot;
ydot(11) = p_f_dot;
ydot(12) = p_s_dot;
ydot(13) = p_K_dot;
ydot(14) = q_1_dot;
ydot(15) = q_2_dot;
ydot(16) = h_IP3_dot;
ydot(17) = P_SOC_dot;
ydot(18) = R_01_dot;
ydot(19) = R_10_dot;
ydot(20) = R_11_dot;
ydot(21) = Ca_u_dot;
ydot(22) = Ca_r_dot;
ydot(23) = Ca_i_dot;
ydot(24) = Na_i_dot;
ydot(25) = K_i_dot;
ydot(26) = Cl_i_dot;
ydot(27) = V_m_dot;

dydt = ydot;

%% Other outputs (saved as global variables)
global tStep tArray_store
global NE_store NO_store K_o_store I_stim_store
global I_CaL_store I_KCa_store I_K_store I_PMCA_store I_NaK_store I_NCX_store
global INa_NSC_store ICa_NSC_store IK_NSC_store I_rel_store I_CaL_H_store

if t > tArray_store(tStep), % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;
end
tArray_store(tStep) = t;

NE_store(tStep) = NE; NO_store(tStep) = NO; K_o_store(tStep) = K_o;
I_stim_store(tStep) = I_stim; I_CaL_store(tStep) = I_Catot;
I_KCa_store(tStep) = I_KCa; I_K_store(tStep) = I_K;
I_PMCA_store(tStep) = I_PMCA; I_NaK_store(tStep) = I_NaK;
I_NCX_store(tStep) = I_NCX; INa_NSC_store(tStep) = INa_NSC;
ICa_NSC_store(tStep) = ICa_NSC; I_rel_store(tStep) = I_rel;
IK_NSC_store(tStep) = IK_NSC; I_CaL_H_store(tStep) = I_CaL;