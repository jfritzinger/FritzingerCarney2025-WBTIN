%% save_AN_models
% This script saves AN model output for broad inhibition model for model
% neurons with CF = 1000, 3000, and 5000 Hz. Used in
% save_model_parameters.m script that runs broad inhibition results and
% those outputs are used in figure 11. 
%
% J. Fritzinger

%% Set up pathing and naming

% Set up paths
[base, datapath, ~, ~] = get_paths();

% Add paths
addpath(fullfile(base, 'scripts', 'helper-functions'), '-end')
addpath(fullfile(base, 'scripts', 'stim-generation'), '-end')
addpath(fullfile(base, 'scripts', 'model-lat-inh'), '-end');

% Choose CF 
putative = 'CF_ 1000'; % Options: 1000, 3000, or 5000 for the manuscript
switch putative
	case 'CF_1000'
		CF = 1000;
	case 'CF_3000'
		CF = 3000;
	case 'CF_5000'
		CF = 5000;
end

% Set up save paths
savepath = fullfile(base, 'model-outputs', putative);
mkdir(savepath)

%% Create stimuli 

% MTF parameters
params{1}.type = 'MTFN';
params{1}.Fs = 100000;
params{1}.mnrep = 10;
params{1}.ramp_dur = 0.05;
params{1}.noise_state = 0;
params{1}.noise_band = [100, 10000];
params{1}.dur = 1;
params{1}.reptim = 1.5;
params{1}.fms = [2, 600, 3];
params{1}.mdepths = [0,0,1];
params{1}.binmode = 2;
params{1}.No = 30;
params{1}.spl = 30;
params{1}.raised_sine = 1;
params{1}.onsetWin = 25;
params{1} = generate_MTF(params{1});
params{1}.num_stim = size(params{1}.stim, 1);

% TIN parameters
params{2}.type = 'TIN';
params{2}.Fs = 100000;
params{2}.CF = CF;
params{2}.version = 5;
params{2}.dur = 0.3;
params{2}.reptim = 0.6;
params{2}.mnrep = 5;
params{2}.ramp_dur= 0.01;
params{2}.No = [3, 23, 43];
params{2}.SNR = [30 40];
params{2}.binmode = 2; % binaural
params{2}.choice = 2; % Frozen noise (1), warm noise (2), all warm (3)
params{2}.condition = 0; % NoSo (0) or NoSpi (1)
params{2}.bw_choice = 1; % NB (1) or WB(2) or 100 Hz(3) or match with sliding WB-TIN (4)
params{2} = generate_TIN(params{2});
params{2}.num_stim = size(params{2}.stim, 1);

% WB-TIN parameters
No = [3, 23, 43];
for ii = 1:3
	params{2+ii}.Fs = 100000;
	params{2+ii}.mnrep = 10;
	params{2+ii}.physio = 0;
	params{2+ii}.SNR = [30 40];
	params{2+ii}.CF = CF;
	params{2+ii}.No = No(ii);
	params{2+ii}.fpeak_mid = CF;
	params{2+ii}.stp_otc = 15;
	params{2+ii}.reptim = 0.6;
	params{2+ii}.dur = 0.3;
	params{2+ii}.ramp_dur = 0.02;
	params{2+ii}.range = 3;
	params{2+ii}.bandwidth = 4;
	params{2+ii}.version = 5;
	params{2+ii} = generate_WBTIN(params{2+ii});
	params{2+ii}.num_stim = size(params{2+ii}.stim, 1);
end

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
paramCF = 0.25:0.25:1.5;

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

		disp(['Stim ' num2str(ist) ' took ' num2str(toc(timerVal3)/60) ' minutes'])
	end
	disp(['AN model took ' num2str(toc(timerVal2)/60) ' minutes'])
end

% Save AN model results
filename = sprintf('%s_AN.mat', putative);
save(fullfile(savepath, filename), 'params', 'AN', ...
	'model_params', '-v7.3')
disp(['All AN models took ' num2str(toc(timerVal)/60) ' minutes'])
