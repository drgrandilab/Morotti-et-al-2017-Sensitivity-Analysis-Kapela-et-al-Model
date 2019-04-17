%% Sensitivity Analyisis - part 3/3
% 
% S. Morotti, M. Nieves-Cintrón, M.A. Nystoriak, M.F. Navedo, E. Grandi.
% Predominant contribution of L-type CaV1.2 channel stimulation to impaired
% intracellular calcium and cerebral artery vasoconstriction in diabetic 
% hyperglycemia. Channels. 2017. doi: 10.1080/19336950.2017.1293220.
% 
% Please, cite the above paper when using this code.

%% main PLS
close all
clear all
clc

color = [0 114 189]/255;

disp('Sensitivity Analysis')
disp('----------------------------------------------------------------')

%% Load parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
% and array 'parameter_names'
% parameter_names = {'ICaL','IBKCa','IK','IKleak','IClCa'...
%     'INaKCl','INaK','INCX','IPMCA','ISOC'...
%     'INSC','IIP3','ISRleak','ISRrel','ISRup'} 

load SA_par_matrix_1000_s0p1

[N_trials N_pars] = size(all_parameters);

%% Load outputs (taken from initial conditions matrix)
% load matrix all_ICs (columns: N state variables, rows: N trials)
load SA_ICs_matrix_1000_s0p1

%              -Em            [Ca2+]i
all_outputs = [-all_ICs(:,27) all_ICs(:,23)];
% the PLS routine works only with positive numbers as inputs!

output_names = {'|Em|','[Ca2+]i'};
output_units = {'mV','mM'};
N_outputs = 2;

%% Outputs baseline model
output_baseline = [-59.4 68]; % mV, nM

disp('Outputs in the baseline model:')
disp(['Em = ',num2str(output_baseline(1),4),' mV'])
disp(['[Ca2+]i = ',num2str(output_baseline(2),4),' nM'])
disp('----------------------------------------------------------------')

%% Outputs - population level
all_outputs_mean = mean(all_outputs);
all_outputs_std_dev = std(all_outputs);

disp('Average outputs at the population level (mean +/- std dev):')
disp(['Em = -',num2str(all_outputs_mean(1),4),' +/- ',num2str(all_outputs_std_dev(1),4),' mV'])
disp(['[Ca2+]i = ',num2str(1e6*all_outputs_mean(2),4),' +/- ',num2str(1e6*all_outputs_std_dev(2),4),' nM'])
disp('----------------------------------------------------------------')

%% Call the PLS routine
good_count = N_trials;
X = log(all_parameters);
Y = log(all_outputs);
% N_pars, N_outputs, N_trials, good_count

% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = PLS_nipals(X,Y,rank(X));

% Regression Coefficients
disp('Regression Coefficients (Em, [Ca2+]i):')
Reg_Coeff = Bpls
disp('----------------------------------------------------------------')

% Calculate agreement of values predicted by regression (Yhat = Bpls*X) 
% with original outputs (Y)
SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
R2each = SSYR./SSYT;

%% Plot Regression Coefficients
figure,set(gcf,'color','w','Position',[100,100,800,400])
for subdex = 1:N_outputs,
	subplot(1,2,subdex)
	% Plot data points
    plot(exp(Y(:,subdex)),exp(Yhat(:,subdex)),'Marker','o','LineStyle','none','Color',color);
    xlabel(['Actual ', output_names{subdex}])
    ylabel(['Predicted ', output_names{subdex}])
    title(['R^2 = ',num2str(R2each(subdex),4)])
    set(gca,'box','off','tickdir','out','fontsize',10)
    % Plot identity line
    g_ylim = get(gca,'ylim') ;
    g_xlim = get(gca,'xlim') ;
    minpoint = min([g_ylim(1),g_xlim(1)]);
    maxpoint = max([g_ylim(2),g_xlim(2)]);
    hold on
    plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
end

figure,set(gcf,'color','w','Position',[100,100,800,400])
for subdex2 = 1:N_outputs,
    subplot(1,2,subdex2)
	bar(Bpls(:,subdex2),'FaceColor',color)
	title(output_names(subdex2))
    set(gca,'box','off','tickdir','out','fontsize',10)
    set(gca,'XTick',1:N_pars)
    set(gca,'XTickLabel',parameter_names)
    set(gca,'XLim',[0 N_pars+1])
    rotateXLabels( gca(), 90)
end

%% Plot Perturbations
scale = (1-0.50:0.01:1+0.50);
xlim1 = 0.7; xlim2 = 1.3;
ylim1 = 0.9; ylim2 = 1.1;

% Em
var_index = 1;
var_value = output_baseline(1); % mV

B1 = Bpls(1,var_index); % LTCC
mod1 = var_value.*scale.^B1;
B2 = Bpls(2,var_index); % BKCa
mod2 = var_value.*scale.^B2;
B3 = Bpls(3,var_index); % KV
mod3 = var_value.*scale.^B3;

figure; set(gcf,'color','w')
subplot(1,2,1),hold on, plot(scale,mod1,scale,mod2,scale,mod3)
set(gca,'box','off','tickdir','out','fontsize',10)
xlabel('Scale factor (-)');
ylabel('Em (mV)');
legend('LTCC','BKCa','KV2.1')
xlim([xlim1 xlim2])
ylim1adj = var_value*ylim1; ylim2adj = var_value*ylim2;
ylim([ylim2adj ylim1adj])
    
% [Ca2+]i
var_index = 2;
var_value = output_baseline(2); % nM

B1 = Bpls(1,var_index); % LTCC
mod1 = var_value.*scale.^B1;
B2 = Bpls(2,var_index); % BKCa
mod2 = var_value.*scale.^B2;
B3 = Bpls(3,var_index); % KV
mod3 = var_value.*scale.^B3;

%figure; set(gcf,'color','w')
set(gcf,'color','w','Position',[100,100,800,400])
subplot(1,2,2),hold on, plot(scale,mod1,scale,mod2,scale,mod3)
set(gca,'box','off','tickdir','out','fontsize',10)
xlabel('Scale factor (-)');
ylabel('[Ca2+]i (nM)');
%legend('LTCC','BKCa','KV2.1')
xlim([xlim1 xlim2])
ylim1adj = var_value*ylim1; ylim2adj = var_value*ylim2;
ylim([ylim1adj ylim2adj])