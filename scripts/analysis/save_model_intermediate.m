%% save_model_intermediate
% Saves broad inhibition output and intermediate cell outputs for plotting
% in figure 10. 
%
% J. Fritzinger, updated 5/17/25
clear 

%% Generate SAM noise

param.type = 'MTFN';
param.Fs = 100000;
param.mnrep = 30;
param.ramp_dur = 0.05;
param.noise_state = 0;
param.noise_band = [100, 10000];
param.dur = 1;
param.reptim = 1.5;
param.fms = [2, 600, 2];
param.mdepths = [0,0,1];
param.binmode = 2;
param.No = 30;
param.spl = 30;
param.raised_sine = 1;
param.onsetWin = 25;
param = generate_MTF(param);
param.num_stim = size(param.stim, 1);

%% Run model

IC = cell(2,1);
for iexample = 1:2
	switch iexample
		case 1
			range = 0.62;
			CS_params = [1 0.28 0];
			BMFs = [87 164 111];
			CF = 1000;
		case 2
			range = 0.5;
			CS_params = [0.19 0.00 0.001];
			BMFs = [300	38	10];
			CF = 3482;
	end

	% Model parameters
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
	model_params.config_type = 'BS inhibited by off-CF BS';

	% Run AN Model
	timerVal = tic;
	model_params.lateral_CF = [CF*2^(-1*range), CF, CF*2^range];
	model_params.CFs = model_params.lateral_CF;
	model_params.CF_range = model_params.CFs(2);

	% Used for manuscript
	AN = modelLateralAN(param, model_params);
	disp(['AN model took ' num2str(toc(timerVal)) ' seconds'])

	% Run IC model
	an_sout = squeeze(AN.an_sout);
	an_sout_lo = squeeze(AN.an_sout_lo);
	an_sout_hi = squeeze(AN.an_sout_hi);
	IC{iexample} = modelLateralSFIE_BMF(param, model_params, ...
		an_sout, an_sout_lo, an_sout_hi, 'CS_params', CS_params,...
		'BMFs', BMFs);
	IC{iexample}.CFs = [CF*2^(-1*range), CF, CF*2^range];
	IC{iexample}.CS_params = CS_params;
	IC{iexample}.BMFs = BMFs;

% Save model responses
[base, datapath, ~, ~] = get_paths();
filename = ['Model_Component_Example_' num2str(iexample) '.mat'];
save(fullfile(datapath, filename), 'param', 'IC', ...
	'model_params', '-v7.3')

end