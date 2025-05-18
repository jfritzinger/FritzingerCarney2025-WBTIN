%% basecamp.m
% This script automates the process of fitting the broad inhibition model
% to a specific neural response or runs the AN model for a specific CF.
% It performs the following steps:
%   1. Sets up required paths and loads neuron-specific data.
%   2. Generates auditory stimuli (MTF, TIN, WB-TIN)
%   3. Runs the Auditory Nerve (AN) model to simulate neural responses
%		to the stimuli.
%   4. Fits an Inferior Colliculus (IC) model to the simulated responses
%		using optimization routines.
%   5. Evaluates model fits by comparing simulated and recorded neural
%		data, and generates summary plots.
%
% To use: Set the desired neuron in the 'Choose neuron' section. For
% running a model with a specific CF, not a putative unit, only the AN
% model will run.
% Output: Model results and evaluation figures are saved to the
% specified output directory.
%
% Author: J. Fritzinger
clear

%% Choose neuron

% Set up paths
[base, datapath, ~, ~] = get_paths();

% Add paths
addpath(fullfile(base, 'scripts', 'helper-functions'), '-end')
addpath(fullfile(base, 'scripts', 'stim-generation'), '-end')
addpath(fullfile(base, 'scripts', 'UR_EAR_2022a'), '-end');
addpath(fullfile(base, 'scripts', 'UR_EAR_2022a', 'source'), '-end');
addpath(fullfile(base, 'scripts', 'model-fitting'), '-end');
addpath(fullfile(base, 'scripts', 'model-lat-inh'), '-end');

% Could put any neuron here
putative = 'R29_TT4_P2_N16';
switch putative
	case 'R29_TT4_P2_N16'
		CF = 2639;
	case 'R29_TT4_P2_N2'
		CF = 3482;
end

% Set up save paths
savepath = fullfile(base, 'model-outputs', putative);
mkdir(savepath)

%% Create stimuli

% Load in data / params
load(fullfile(datapath, 'Neural_Data', [putative '.mat']), 'data');

% MTFN parameters
params{1} = data{3,2}; % Gets binaural WB-TIN stimuli
params{1}.Fs = 100000;
params{1}.dur = 0.7;
params{1}.mnrep = 5;
params{1} = generate_MTF(params{1});
params{1}.num_stim = size(params{1}.stim, 1);

% TIN parameters
[~,closestIndex] = min(abs(data{9,2}.fpeaks-CF));
params{2}.type = 'TIN';
params{2}.Fs = 100000;
params{2}.CF = data{9,2}.fpeaks(closestIndex);
params{2}.version = 5;
params{2}.dur = 0.3;
params{2}.reptim = 0.6;
params{2}.mnrep = 5;
params{2}.ramp_dur= 0.01;
params{2}.No = 23;
params{2}.SNR = [30 40];
params{2}.binmode = 2; % binaural
params{2}.choice = 2; % Frozen noise (1), warm noise (2), all warm (3)
params{2}.condition = 0; % NoSo (0) or NoSpi (1)
params{2}.bw_choice = 1; % NB (1) or WB(2) or 100 Hz(3) or match with sliding WB-TIN (4)
params{2} = generate_TIN(params{2});
params{2}.num_stim = size(params{2}.stim, 1);

% WB-TIN parameters
params{3} = data{7,2}; % Gets binaural WB-TIN stimuli
params{3}.Fs = 100000;
params{3}.mnrep = 5;
params{3}.physio = 0;
params{3}.SNR = 40;
params{3}.CF = CF;
params{3} = generate_WBTIN(params{3});
params{3}.num_stim = size(params{3}.stim, 1);

%% Run AN model

paramCF = 0.125:0.125:1.5;
run_AN_model(params, CF, paramCF, putative, savepath);


%% Run IC fitting procedure

run_IC_model_fmincon(savepath, putative, CF)

%% Evaluate

% Load in neural data
[~, datapath, ~, ~] = get_paths();
load(fullfile(datapath, 'Neural_Data', [putative '.mat']), 'data');
data_rates = analyze_data(data, CF); % Analyze data and put in correct form

% Load in AN
filename = sprintf('%s_AN.mat', putative);
load(fullfile(savepath, filename), 'params', 'AN', 'model_params')

% Load in IC parameter values
filename = sprintf('%s_IC.mat', putative);
load(fullfile(savepath, filename), 'fit_params_all')

% Evaluate fits
evaluate_fits(CF, data_rates, params, AN, model_params, fit_params_all)
get_best_fit_model(putative, data_rates, params, AN, model_params, fit_params_all, savepath)
pdf_evaluate_fits(putative, CF, data_rates, params, AN, model_params, fit_params_all, savepath)
