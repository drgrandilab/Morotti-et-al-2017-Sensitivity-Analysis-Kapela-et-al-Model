Morotti et al. code for linear regression analysis of populations of electrophysiological models.

This package contains the Matlab code used to perform sensitivity analysis on a population of electrophysiological
models with the method of linear logistic regression.
Here, the Kapela et al. model of rat mesenteric smooth muscle cell (Journal of Theoretical Biology 253 (2008)
238– 260) is used to create a family of 1000 model variants by perturbing model parameters, and to investigate
how model parameters influence membrane potential and calcium concentration.


________________________________________________________________________________________________________________
Contents:

readme.txt		 	this file

	.m files - model components

SMC_masterCompute.m		loads initial conditions and runs the simulations with the baseline model
SMC_eccODEfile.m		contains ODEs for the Kapela et al. model

SMC_SA1_generate_parameters.m	generates random perturbations required for building a population of models
SMC_SA2_obtain_ICs.m		generates steady-state conditions for each model in the population
SMC_SA3_sensitivity_analysis.m	performs linear regression analysis and plots the results		

PLS_nipals.m			function used for linear regression analysis
rotateXLabels.m			function used for plotting the results

	.mat files - initial conditions (ICs)

yfin_icaMbar_rest.mat	 	control model (steady-state, GLU = 0 in SMC_eccODEfile.m)
yfin_icaMbar_rest_GLU.mat	HG model (steady-state, GLU = 1 in SMC_eccODEfile.m)

	.mat files - others

SA_par_matrix_1000_s0p1.mat	matrix with parameter perturbations in each element of the population
SA_ICs_matrix_1000_s0p1.mat	matrix with ICs for each element of the population
________________________________________________________________________________________________________________


Reference:
S. Morotti, M. Nieves-Cintrón, M.A. Nystoriak, M.F. Navedo, E. Grandi.
Predominant contribution of L-type CaV1.2 channel stimulation to impaired intracellular calcium and cerebral
artery vasoconstriction in diabetic hyperglycemia.
Channels (Austin). 2017 Jul 4;11(4):340-346. doi: https://doi.org/10.1080/19336950.2017.1293220

Please, cite the above paper when using these codes.
