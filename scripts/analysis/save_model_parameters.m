%% save_model_parameters.m
% This struct saves the broad inhibition model output for figure 11, which
% details how the WB-TIN and MTF responses change with changing model
% parameters. 
%
% J. Fritzinger
clear

%%
timerVal = tic;

% Load in data
[~, datapath] = get_paths();
CF = 3000; 
filename = sprintf('CF_%d_AN.mat', CF);
load(fullfile(datapath, sprintf('CF_%d', CF), filename), 'params', 'AN')

%% Run model - strength parameter
paramS = [0.1 0.3 0.5];
paramD = 0;
iparamCF = 4;

model = cell(3, 1);
for imodel = 1:3
	model{imodel} = run_model(params, AN, iparamCF, paramS(imodel), paramD);
end
save(fullfile(datapath, 'ModelParams1.mat'),...
	"model", "params", "paramS", "paramD", "iparamCF", "CF")
elapsedTime = toc(timerVal)/60;
disp(['Stim took ' num2str(elapsedTime) ' minutes'])

%% Run model - CF parameter
paramS = 0.4;
paramD = 0;
iparamCF = [1 3 5];

model = cell(3, 1);
for imodel = 1:3
	model{imodel} = run_model(params, AN, iparamCF(imodel), paramS, paramD);
end
save(fullfile(datapath, 'ModelParams2.mat'),...
	"model", "params", "paramS", "paramD", "iparamCF", "CF")
elapsedTime = toc(timerVal)/60;
disp(['Stim took ' num2str(elapsedTime) ' minutes'])

%% Run model - delay parameter
% Not running this, delay doesn't have any effect on WB-TIN response

% paramS = 0.3;
% paramD = [0 0.001 0.002];
% iparamCF = 1;
% for imodel = 1:3
% 	model{imodel} = run_model(params, AN, iparamCF, paramS, paramD(imodel));
% end
% save(fullfile(datapath, 'ModelParams3.mat'),...
% 	"model", "params", "paramS", "paramD", "iparamCF", "CF")
% elapsedTime = toc(timerVal)/60;
% disp(['Stim took ' num2str(elapsedTime) ' minutes'])

%% Run model - different CFs
paramS = 0.4;
paramD = 0;
iparamCF = 4;
name = 'CF';
CF = [1000 3000 5000]; % 1kHz may be better

model = cell(3, 1);
for imodel = 1:3
	filename = sprintf('%d_AN.mat', CF(imodel));
	load(fullfile(filepath, sprintf('CF_%d', CF(imodel)), filename), 'params', 'AN')
	model{imodel} = run_model(params, AN, iparamCF, paramS, paramD);
	if imodel == 1
		params1 = params;
	elseif imodel == 2
		params2 = params;
	else
		params3 = params;
	end
	clear params
	clear AN
end
save(fullfile(datapath, 'ModelParams4.mat'),...
	"model", "paramS", "paramD", "iparamCF", "CF","params1", ...
	"params2", "params3")
elapsedTime = toc(timerVal)/60;
disp(['Stim took ' num2str(elapsedTime) ' minutes'])