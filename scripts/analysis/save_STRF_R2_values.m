%% save_STRF_R2_values.m
% This script loads in each WB-TIN response and STRF and creates a WB-TIN
% prediction based on the STRF. It saves a .mat file called 'STRFModel_23'
% that is used in figure 6. Takes a long time to run. 
%
% J. Fritzinger
clear 

%% Load in data

[~, datapath] = get_paths();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), ...
	'PreserveVariableNames',true);
num_units = size(sessions, 1);


%% Find sessions of interest

WB_noise_con(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_3);
WB_noise_con(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_23);
WB_noise_con(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_43);
WB_noise(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_3);
WB_noise(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_23);
WB_noise(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_43);
STRF_bin = cellfun(@(s) contains(s, 'R'), sessions.STRF);
STRF_con = cellfun(@(s) contains(s, 'R'), sessions.STRF_con);

bin = WB_noise(:,2) & STRF_bin;
contra = WB_noise_con(:,2) & STRF_con;

%% Analysis
binmodes = {'Contra', 'Binaural'};
Nos = {'3', '23', '43'};
ind = 1;
temp = struct();
for ibin = 2

	if ibin == 1
		index = find(contra);
	else
		index = find(bin);
	end

	for isesh = 1:length(index)

		s_ind = index(isesh);
		putative_neuron = sessions.Putative_Units{s_ind};
		CF = sessions.CF(s_ind);
		load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');

		% Get data for each stimulus
		params_STRF = data{4,ibin};
		param_WB = data(7,ibin); % Gets binaural WB-TIN stimuli, 7 is 23 dB, 8 is 43 dB

		% General analysis
		data_STRF = analyzeSTRF(params_STRF);
		data_WB = analyzeWBTIN(param_WB, []);
		param_WB = param_WB{1};
		data_WB = data_WB{1};

		% Calculate model response
		[R2, avModel, stdModel, ratio, max_all] = modelWBTINSTRF(param_WB,...
			data_STRF, data_WB);

		% Display progress
		fprintf('%s Done, %.2f through %s \n', putative_neuron,...
			isesh/length(index), binmodes{ibin})

		% Add to struct
		temp(isesh).putative = putative_neuron;
		temp(isesh).R2 = R2;
		temp(isesh).avModel = avModel;
		temp(isesh).stdModel = stdModel;
		temp(isesh).ratio = ratio;
		temp(isesh).max_all = max_all;
		temp(isesh).WB_No = param_WB.No;
		temp(isesh).SNRs = param_WB.SNR;

		if ibin == 1
			STRFmodel_con = temp(isesh);
			filename = [putative_neuron '_STRF_Contra_23.mat'];
			save(fullfile(datapath, filename), "STRFmodel_con")
		else
			STRFmodel = temp(isesh);
			filename = [putative_neuron '_STRF_Bin_23.mat'];
			save(fullfile(datapath, filename), "STRFmodel")
		end
	end

	%% Save Data

	if ibin == 1
		STRF_contra = temp;
		save(fullfile(datapath, 'STRFModel_23_contra.mat'), "STRF_contra")
	else
		STRF = temp;
		save(fullfile(datapath, 'STRFModel_23.mat'), "STRF")
	end
	clear temp

end

%% Load in each and save as one matrix - used if error occurs, because
% this script takes a long time to run

% binmodes = {'Contra', 'Binaural'};
% Nos = {'3', '23', '43'};
% ind = 1;
% for ibin = 2
% 	if ibin == 1
% 		index = find(contra);
% 	else
% 		index = find(bin);
% 	end
% 
% 	for isesh = 1:length(index)
% 		s_ind = index(isesh);
% 
% 		% Load in session
% 		putative_neuron = sessions.Putative_Units{s_ind};
% 		CF = sessions.CF(s_ind);
% 		if ibin == 1
% 			filename = [putative_neuron '_STRF_Contra_23.mat'];
% 			load(fullfile(datapath, filename), "STRFmodel_con")
% 		else
% 			filename = [putative_neuron '_STRF_Bin_23.mat'];
% 			load(fullfile(datapath, filename), "STRFmodel")
% 		end
% 		temp(isesh) = STRFmodel; % Add to struct
% 	end
% 
%     % Save Data
% 	if ibin == 1
% 		STRF_contra = temp;
%         save(fullfile(datapath, 'STRFModel_23_contra.mat'), "STRF_contra")
% 	else
% 		STRF = temp;
%         save(fullfile(datapath, 'STRFModel_23.mat'), "STRF")
% 	end
% end
