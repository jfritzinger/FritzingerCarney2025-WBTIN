function run_AN_model(params, CF, paramCF, putative, savepath)
% run_AN_model  Simulate auditory nerve (AN) responses for multiple CF
% ranges and stimuli
%
%   Inputs:
%     params    - Cell array of stimulus parameter structures
%     CF        - Characteristic frequency of the neuron (Hz)
%     paramCF   - Vector of CF offsets (in octaves) to simulate lateral interactions
%     putative  - String identifier for the neuron (used in output filename)
%     savepath  - Directory path where results will be saved
%
%   Outputs:
%     None (results are saved as a .mat file in savepath)
%
%   The function configures model parameters (species, fiber type, etc.),
%   runs the model for each combination of CF offset and stimulus,
%   and saves the resulting AN structure and parameters.
%
%   Author: J. Fritzinger

%% Model parameters

model_params.range = 2; % 1 = population model, 2 = single cell model
model_params.species = 1; % 1 = cat, 2 = human
model_params.num_CFs = 1;
model_params.nAN_fibers_per_CF = 10;
model_params.cohc = 1; % (0-1 where 1 is normal)
model_params.cihc = 1; % (0-1 where 1 is normal)
model_params.nrep = 10; % how many times to run the AN model
model_params.fiberType = 3; % AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR)
model_params.implnt = 1; % 0 = approximate model, 1=exact powerlaw implementation(See Zilany etal., 2009)
model_params.noiseType = 1; % 0 = fixed fGn, 1 = variable fGn) - this is the 'noise' associated with spont. activity of AN fibers - see Zilany et al., 2009. "0" lets you "freeze" it.
model_params.which_IC = 1; % 2 = ModFilt; 1 = SFIE model
model_params.onsetWin = 0.020; % exclusion of onset response, e.g. to omit 1st 50 ms, use 0.050
model_params.BMF = 100;
model_params.type = 'Lateral Model';
model_params.config_type = 'BS inhibited by off-CF BS';

% Run AN Model with Multiple CFs
timerVal = tic;
num_stim = size(params, 2);
model_params.paramCF = paramCF;
num_paramCF = length(paramCF);
AN = cell(1, num_stim);
for iparamCF = 1:num_paramCF
	timerVal2 = tic;

	% Run AN Model
	model_params.lateral_CF = [CF*2^(-1*paramCF(iparamCF)), CF, CF*2^paramCF(iparamCF)];
	model_params.CFs = model_params.lateral_CF;
	model_params.CF_range = model_params.CFs(2);

	for ist = 1:num_stim
		timerVal3 = tic;

		% Used for manuscript
		AN{iparamCF, ist} = modelLateralAN(params{ist}, model_params);
		AN{iparamCF, ist}.CF_span = paramCF(iparamCF);

		% Test with efferents
		%efferent = 0;
		%AN{ist} = modelLateralEfferentAN(params{ist}, model_params, efferent);

		disp(['Stim ' num2str(ist) ' took ' num2str(toc(timerVal3)/60) ' minutes'])
	end
	disp(['AN model took ' num2str(toc(timerVal2)/60) ' minutes'])
end

% Save AN model results
filename = sprintf('%s_AN.mat', putative);
save(fullfile(savepath, filename), 'params', 'AN', ...
	'model_params', '-v7.3')
disp(['All AN models took ' num2str(toc(timerVal)/60) ' minutes'])

end