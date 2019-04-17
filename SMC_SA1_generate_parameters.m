%% Sensitivity Analyisis - part 1/3
% 
% S. Morotti, M. Nieves-Cintrón, M.A. Nystoriak, M.F. Navedo, E. Grandi.
% Predominant contribution of L-type CaV1.2 channel stimulation to impaired
% intracellular calcium and cerebral artery vasoconstriction in diabetic 
% hyperglycemia. Channels. 2017. doi: 10.1080/19336950.2017.1293220.
% 
% Please, cite the above paper when using this code.

%% main generates random parameters
close all
clear all
clc

%% Parameters
% 1) ICaL
% 2) IBKCa
% 3) IK
% 4) IKleak
% 5) IClCa
% 6) INaKCl
% 7) INaK
% 8) INCX
% 9) IPMCA
% 10) ISOC
% 11) INSC
% 12) IIP3
% 13) ISRleak
% 14) ISRrel
% 15) ISRup

parameter_names = {'ICaL','IBKCa','IK','IKleak','IClCa'...
    'INaKCl','INaK','INCX','IPMCA','ISOC'...
    'INSC','IIP3','ISRleak','ISRrel','ISRup'} ;

n_parameters = length(parameter_names);
baseline_parameters = ones(1,n_parameters);

%% Random variations
variations = 1000; % number of trials

sigmaG = 0.1*ones(1,n_parameters); % standard deviation for each parameter

all_parameters = zeros(variations,n_parameters);
for ii = 1:n_parameters,
    scaling = exp(sigmaG(ii)*randn(1,variations)) ;
    newparams = baseline_parameters(ii)*scaling ;
    all_parameters(:,ii) = newparams ;
end

all_parameters;
% columns: N parameters
% rows: N trials

%% Saving
%save SA_par_matrix_1000_s0p1 all_parameters parameter_names