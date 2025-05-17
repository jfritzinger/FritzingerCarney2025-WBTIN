%% save_predictable_variance
%
% This script loads in the putative neurons spreadsheet and calculates the
% predictable variance for each WB-TIN presentation for each neuron. 
%
% J. Fritzinger, updated 5/17/25

%% Load in spreadsheet

[base, datapath, savepath, ppi] = get_paths();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);
num_units = size(sessions, 1);

%% Analyze predictable variance

WB_noise_con(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_3);
WB_noise_con(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_23);
WB_noise_con(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_43);
WB_noise(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_3);
WB_noise(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_23);
WB_noise(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_43);

timerVal = tic;
binmodes = {'Contra', 'Binaural'};
Nos = {'3', '23', '43'};
Vp_all = cell(2, 3);
for ibin = 2
	for iNo = 2

		if ibin == 1
			index = find(WB_noise_con(:,iNo));
		else
			index = find(WB_noise(:,iNo));
		end

		num_sesh = length(index);
		Vp_section = NaN(num_sesh,3);
		parfor isesh = 1:num_sesh
			s_ind = index(isesh);

			% Load in session
			putative_neuron = sessions.Putative_Units{s_ind};
			d = load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']));

			% Get data for each stimulus
			params_WB = d.data(5+iNo,ibin); % Gets binaural WB-TIN stimuli

			% General analysis
			data_WB = analyzeWBTIN(params_WB, []);

			% Save data into struct that will be plotted
			SNR_vals = [20 30 40];
			if isempty(data_WB{1}.fpeaks)
				data_WB{1}.V_p = NaN(1,3);
			elseif length(data_WB{1}.V_p) == 2
				which_exists = ismember(SNR_vals, data_WB{1}.SNRs);
				temp = NaN(1,3);
				temp(which_exists) = data_WB{1}.V_p;
				data_WB{1}.V_p = temp;
			end
			Vp_section(isesh,:) = data_WB{1}.V_p;
		end
		Vp_all{ibin, iNo} = Vp_section;
	end
end

Vp_combined = [];
for ii = 1:6
	Vp_combined = [Vp_combined; Vp_all{ii}];
end

% Save 
save(fullfile(datapath, 'FigS1_Vp.mat'), 'Vp_combined')
elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])