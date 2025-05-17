function lateral_model = run_model(params, AN, iparamCF, paramS, paramD)
% run_model is a function that is used in save_model_parameters.m and runs
% the lateral inhibition / broad inhibition model with specific parameter
% changes for figure 11. 
%
% J. Fritzinger, updated 5/17/25


% Model parameters
model_params.config_type = 'BS inhibited by off-CF BS';
model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.num_CFs = 1;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 1; % how many times to run the AN model
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.BMF = 100;

model_params.type = 'Lateral Model';
lm_params = [paramS paramS paramD];
num_stim = size(params,2);
for ist = 1:num_stim
	lateral_model{ist} = modelLateralSFIE(params{ist}, ...
		model_params, AN{iparamCF, ist}.an_sout, AN{iparamCF, ist}.an_sout_lo, ...
		AN{iparamCF, ist}.an_sout_hi,'CS_params', lm_params);
end
end