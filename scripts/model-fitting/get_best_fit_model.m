function get_best_fit_model(putative, data_rates, params, AN, ...
	model_params, fit_params_all, modelpath)
% get_best_fit_model identifies and saves the best-fitting IC model 
% parameters for a neuron
%
%   Inputs:
%     putative       - String identifier for the neuron (used in filenames)
%     data_rates     - Experimental spike rates for all stimulus conditions
%     params         - Cell array of stimulus parameter structures
%     AN             - Auditory nerve model outputs (cell array)
%     model_params   - Structure of lateral inhibition model parameters
%     fit_params_all - Matrix of fitted IC model parameters (rows: CF ranges)
%     modelpath      - Directory path where the best model will be saved
%
%   Outputs:
%     None (saves best-fitting AN/IC model and parameters as a .mat file in modelpath)
%
%   Author: J. Fritzinger


%% Evaluate model fits  

num_paramCF = size(AN, 1);
mse_all = NaN(num_paramCF, 1);
for iparamCF = 1:num_paramCF

	% Run model with fit parameters
	AN_sub = AN(iparamCF,:);
	fit_params = fit_params_all(iparamCF,:);
	CS_params = [fit_params(1:2) 0.001];
	BMFs = fit_params(3:5);
	nstim = size(params, 2);
	model_outputs = cell(nstim, 1);
	for istim = 1:nstim
		param = params{istim};
		an_sout = squeeze(AN_sub{istim}.an_sout);
		an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
		an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
		model_outputs{istim} = modelLateralSFIE_BMF(param, model_params, ...
			an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
			'BMFs', BMFs);
	end

	% Calculate MSE 
	mse = minimize_IC_model_fit(data_rates, AN_sub, params, model_params,...
		fit_params);
	mse_all(iparamCF) = mse;

end

% Get best parameters
[~, min_ind] = min(mse_all);


%% Run model with best parameters
 
paramCF = min_ind;
AN_sub = AN(paramCF,:);
fit_params = fit_params_all(paramCF,:);
CS_params = [fit_params(1:2) 0.001];
BMFs = fit_params(3:5);
nstim = size(params, 2);
IC_best = cell(nstim, 1);
for istim = 1:nstim
	param = params{istim};
	an_sout = squeeze(AN_sub{istim}.an_sout);
	an_sout_lo = squeeze(AN_sub{istim}.an_sout_lo);
	an_sout_hi = squeeze(AN_sub{istim}.an_sout_hi);
	IC_best{istim} = modelLateralSFIE_BMF(param, model_params, ...
		an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
		'BMFs', BMFs);
end
AN_best = AN(paramCF, :);

%% Save model and best parameters

filename = sprintf('%s_BestModel.mat', putative);
save(fullfile(modelpath, filename),...
	'AN_best', "IC_best", "model_params", "params", "fit_params")

