%% save_dog_spreadsheet
%
% Script that loads in a WBTIN population and saves a population data table
% that has all of the DoG parameters in it
%
% J. Fritzinger, updated 5/17/25
clear

%% Load in data
timerVal = tic;

[base, datapath, savepath, ~] = get_paths();
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
contra = WB_noise_con(:,1) | WB_noise_con(:,2) | WB_noise_con(:,3);
bin = WB_noise(:,1) | WB_noise(:,2) | WB_noise(:,3);

%% Get all WBTIN data & add DoG fit to all
filename = 'DoGAnalysis.xlsx';

% Create Table
varNames = ["Putative", "MTF_shape", "DSID", ...
	"No", "NoiseBW", "SNR", "Binmode", "CF", "RM", "g_exc", "g_inh", ...
	"g_ratio", "sigma_exc", "sigma_inh", "sigma_ratio", "CF_exc", ...
	"CF_inh", "offset_diff", "R2", "R2_er", "SNR_er", "V_p", "ratio1"];
varTypes = ["string", "double", "double", "double", "double","double",...
	"double", "double","string", "double","double","double","double",...
	"double", "double", "double", "double","double", "double", "double", ...
	"double", "double", "double"];
table_size = [3000 length(varNames)];
analysisTable = table('Size',table_size,'VariableTypes',varTypes, 'VariableNames',varNames);

iind = 1;
for ibin = 1:2
	if ibin == 1
		index = find(contra);
	else
		index = find(bin);
	end
	num_neurons = length(index);

	for isesh = 1:num_neurons

		% Load in session
		putative_neuron = sessions.Putative_Units{index(isesh)};
		CF = sessions.CF(index(isesh));
		RM = sessions.RM_Type{index(isesh)};
		MTF = sessions.MTF{index(isesh)};
		load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');
		data_all = data;

		if strcmp(MTF, 'BE')
			MTF_type = 1;
		elseif strcmp(MTF, 'BS')
			MTF_type = 2;
		elseif contains(MTF, 'H')
			MTF_type= 3;
		else
			MTF_type = 4;
		end

		% Get data for each stimulus
		params_WB = data_all(6:8,ibin); % Gets binaural WB-TIN stimuli
		emptyCells = cellfun(@isempty,params_WB);
		params_WB = params_WB(~emptyCells);
		data_WB = analyzeWBTIN(params_WB, CF);

		% Loop through the number of WB-TIN datasets for a session
		num_ds = length(data_WB);
		for ids = 1:num_ds

			param = params_WB{ids};
			data = data_WB{ids};

			% If there is no WBTIN data (spikes didn't match or some other error)
			if isempty(data) || isempty(data.fpeaks)
				continue
			end

			% Create DoG Fit
			type = 'Log';
			mse_type = 2;
			[DOGparams, R2, noise_alone_norm_all, max_rate_all] =...
				fitDifferenceofGaussians(data, mse_type);

			% Save DoG parameters for every SNR
			num_SNRs = length(param.SNR);
			for isnr = 1:num_SNRs
				noise_alone_norm = noise_alone_norm_all(isnr);
				max_rate = max_rate_all(isnr);

				analysisTable.Putative{iind} = putative_neuron;
				analysisTable.DSID(iind) = param.dsid;
				analysisTable.No(iind) = param.No;
				analysisTable.NoiseBW(iind) = param.bandwidth;
				analysisTable.Binmode(iind) = param.binmode;
				analysisTable.CF(iind) = CF;
				analysisTable.RM{iind} = RM;
				analysisTable.SNR(iind) = param.SNR(isnr);
				analysisTable.MTF_shape(iind) = MTF_type;

				% Save found parameters % manuscript
				analysisTable.g_exc(iind) = DOGparams(isnr, 1);
				analysisTable.g_inh(iind) = DOGparams(isnr, 2);
				analysisTable.g_ratio(iind) = DOGparams(isnr, 2)/DOGparams(isnr, 1);
				analysisTable.sigma_exc(iind) = DOGparams(isnr, 3);
				analysisTable.sigma_inh(iind) = DOGparams(isnr, 4);
				analysisTable.sigma_ratio(iind) = DOGparams(isnr, 4)/DOGparams(isnr, 3); % inh/exc
				analysisTable.CF_exc(iind) = DOGparams(isnr, 5);
				analysisTable.CF_inh(iind) = DOGparams(isnr, 6);
				analysisTable.offset_diff(iind) = DOGparams(isnr, 5)-DOGparams(isnr, 6);

				% Calculate Coefficient of Determination & Ratios
				fpeaks_dog = log10(data.fpeaks(2:end,1));
				[dog_rate, ~, ~] = createDoG(DOGparams(isnr, :), fpeaks_dog);
				dog_rate = (dog_rate+noise_alone_norm).*max_rate;

				[hat_r2er, ~] = r2er_n2m(dog_rate', data.rate_matrix(:,2:end));
				SNR_er = snr_er_est(data.rate_matrix(:,2:end));
				r_v_ratio = R2(isnr) / data.V_p(isnr);

				analysisTable.R2(iind) = R2(isnr);
				analysisTable.R2_er(iind) = hat_r2er;
				analysisTable.SNR_er(iind) = SNR_er;
				analysisTable.V_p(iind) = data.V_p(isnr);
				analysisTable.ratio1(iind) = r_v_ratio;
				iind = iind +1;

			end
		end
	end
end

%% Save analysisTable

writetable(analysisTable, fullfile(datapath, filename));

elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])
