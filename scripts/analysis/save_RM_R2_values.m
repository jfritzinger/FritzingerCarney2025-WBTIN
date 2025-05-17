%% save_RM_R2_values.m 
% This script loads in all neurons and caluculates the variance explained
% (R2) of the WB-TIN and response map. Saves a spreadsheet called
% 'RM_comparison' that is used in a supplemental figure. 
%
% J. Fritzinger
clear

%% Load in spreadsheet
timerVal = tic;

[base, datapath, ~, ~] = get_paths();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_units = size(sessions, 1);


%% Find sessions of interest

WB_noise_con(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_3);
WB_noise_con(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_23);
WB_noise_con(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_43);
WB_noise(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_3);
WB_noise(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_23);
WB_noise(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_43);
RM_bin = cellfun(@(s) contains(s, 'R'), sessions.type_RM);
RM_con = cellfun(@(s) contains(s, 'R'), sessions.type_RM_con);
bin = WB_noise(:,3) & RM_bin;
contra = WB_noise_con(:,3) & RM_con;

%% Set up table

filename = 'RM_Comparison.xlsx';

% Create Table
varNames = ["Putative", "binmode", "RM", "R2_RM", "BF_spl", "iBF", "SNR", ...
	"WB_No", "WB_No_overall"];
varTypes = ["string", "double", "string", "double", "double", "double",...
	"double", "double", "double"];
table_size = [3000 length(varNames)];
analysisTable = table('Size',table_size,'VariableTypes',varTypes, ...
	'VariableNames',varNames);


%% Analysis
binmodes = {'Contra', 'Binaural'};
Nos = {'3', '23', '43'};
iind = 1;
for ibin = 1:2

	if ibin == 1
		index = find(contra);
	else
		index = find(bin);
	end
	
	for isesh = 1:length(index)
		s_ind = index(isesh);

		% Load in session
		putative_neuron = sessions.Putative_Units{s_ind};
		CF = sessions.CF(s_ind);
		RM_type = sessions.RM_Type(s_ind);
		load(fullfile(datapath,'Neural_Data', [putative_neuron '.mat']), 'data');

		% Get data for each stimulus
		params_RM = data{2,ibin};
		params_WB = data(6:8,ibin); % Gets binaural WB-TIN stimuli
		emptyCells = cellfun(@isempty,params_WB);
		params_WB(emptyCells) = [];

		% General analysis
		data_RM = analyzeRM(params_RM);
		datas_WB = analyzeWBTIN(params_WB, []);

		for iNo = 1:length(datas_WB)
			param_WB = params_WB{iNo};
			data_WB = datas_WB{iNo};

			% Calculate model response
			SNRs = param_WB.SNR;
			snr_choice = 1:length(SNRs);
			[R2, BF_spl, iBF] = calculateRMR2(param_WB, data_WB, data_RM, snr_choice);

			for isnr = 1:length(SNRs)
				analysisTable.Putative{iind} = putative_neuron;
				analysisTable.binmode(iind) = ibin;
				analysisTable.R2_RM(iind) = R2(isnr);
				analysisTable.BF_spl(iind) = BF_spl;
				analysisTable.iBF(iind) = iBF;
				analysisTable.SNR(iind) = SNRs(isnr);
				analysisTable.RM{iind} = RM_type{1};
				analysisTable.WB_No(iind) = param_WB.No;
				analysisTable.WB_No_overall(iind) = param_WB.No_overall;
				iind = iind +1;
			end
		end

		% Display progress
		fprintf('%s Done, %.2f through %s \n', putative_neuron, isesh/length(index), binmodes{ibin})
	end
end

%% Save Table 

writetable(analysisTable, fullfile(datapath, filename));
elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])
