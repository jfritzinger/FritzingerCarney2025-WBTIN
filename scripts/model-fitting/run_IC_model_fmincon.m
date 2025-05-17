function run_IC_model_fmincon(savepath, putative, CF)
% run_IC_model_fmincon  Fit IC model parameters to neural and AN data 
% using fmincon optimization
%
%   Inputs:
%     savepath   - Directory path where AN model results are stored and where
%                  IC fitting results will be saved
%     putative   - String identifier for the neuron (used in filenames)
%     CF         - Characteristic frequency of the neuron (Hz)
%
%   Outputs:
%     None (saves optimized IC model parameters as a .mat file in savepath)
%
%   The function:
%     - Loads neural data and AN model outputs
%     - For each AN CF range, fits IC model parameters (strengths and BMFs)
%       using constrained nonlinear optimization (fmincon)
%     - Saves all fitted parameter sets to disk for further evaluation
%
%   Author: J. Fritzinger

%% Load in data

% Load in neural data
[~, datapath, ~, ~] = get_paths();
load(fullfile(datapath, 'Neural_Data', [putative '.mat']), 'data');
data_rates = analyze_data(data, CF); % Analyze data and put in correct form

% Load in AN model 
filename = sprintf('%s_AN.mat', putative);
load(fullfile(savepath, filename), 'params', 'AN', 'model_params')


%% For each AN response, fit IC model to data
timerVal2 = tic;

num_paramCF = size(AN, 1);
fit_params_all = zeros(num_paramCF,5);
for iparamCF = 1:num_paramCF
	AN_sub = AN(iparamCF,:);
	best_fval = Inf;
	for ipass = 1
		timerVal = tic;

		S1_init = 0.1 + (0.7 - 0.1) * rand(1);
		S2_init = 0.1 + (0.7 - 0.1) * rand(1);
		BMF1_init = 25 + (250 - 25) * rand(1);
		BMF2_init = 25 + (250 - 25) * rand(1);
		BMF3_init = 25 + (250 - 25) * rand(1);

		%x0 = [0.5 0.5 100 100 100]; % Slow, Shigh, D, BMFlow, BMFon, BMFhigh
		x0 = [S1_init S2_init BMF1_init BMF2_init BMF3_init];
		lb = [0   0   10  10  10 ];
		ub = [1   1   300 300 300];

		fun = @(x) minimize_IC_model_fit(...
			data_rates, AN_sub, params, model_params, x);
		options = optimoptions('fmincon','Display','off',...
			'Algorithm','sqp');
		[fit_params, fval] = fmincon(...
			fun, x0, [], [], [], [], lb, ub, [], options);
		

		if fval < best_fval
			best_x = fit_params;
			best_fval = fval;
		end

		fprintf('Fitting took %0.02f minutes\n', toc(timerVal)/60);
	end
	fit_params_all(iparamCF, :) = best_x;
end

% Save IC model fit results
filename = sprintf('%s_IC.mat', putative);
save(fullfile(savepath, filename), 'fit_params_all')
fprintf('Fitting all models took %0.02f minutes\n', toc(timerVal2)/60);

end
