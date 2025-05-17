%% Fig5_LME
% J. Fritzinger, updated 5/20/24
%
% This script runs the linear mixed effects model on the WB-TIN population
% data to get population-level trends
clear

%% Load in data
timerVal = tic;

[base, datapath, savepath, ppi] = get_paths();

% Load in spreadsheet
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), 'PreserveVariableNames',true);

%% Find sessions of interest

WB_noise_con(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_3);
WB_noise_con(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_23);
WB_noise_con(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_con_43);
WB_noise(:,1) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_3);
WB_noise(:,2) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_23);
WB_noise(:,3) = cellfun(@(s) contains(s, 'R'), sessions.SPEC_slide_WB_noise_43);

contra = WB_noise_con(:,1) | WB_noise_con(:,2) | WB_noise_con(:,3);
bin = WB_noise(:,1) | WB_noise(:,2) | WB_noise(:,3);

%MTF_types = sessions.MTF;
MTF_types = sessions.MTF_at_100;
RM_types = sessions.RM_Type;


%% Create table

varNames = ["Putative", "MTF", "CF", "CF_Group", "RM", "No", "SNR", ...
	"Binmode", "Tone_Freq", "Rate"];
varTypes = ["string", "string", "double", "string", "string", "double", "double", ...
	"string", "string", "double"];
table_size = [21000 length(varNames)];
analysisTable = table('Size',table_size,'VariableTypes',varTypes, 'VariableNames',varNames);

%% Two photon plot of all WB-TIN

binmodes = {'Contra', 'Binaural'};

SNRs = [20, 30, 40];
num_SNRs = length(SNRs);

Nos = [3, 23, 43];
num_Nos = length(Nos);

CF_groups = {'Low', 'Medium', 'High'}; % Low < 2kHz, Medium 2-4kHz, High > 4kHz

index_bin = find(bin);
WB_mat = NaN(num_Nos, 3, length(index_bin),10000); % 3 No x 3 SNR x # neurons 
index_con = find(contra);
WB_mat_contra = NaN(num_Nos, 3, length(index_con),10000);
iind = 1;
for ibin = 1:2

	if ibin == 1
		index = index_con;
	else
		index = index_bin;
	end
	num_sesh = length(index);
	V_p = NaN(num_Nos, num_sesh);

	% Sort by CF & get RM/MTF type
	CFs = sessions.CF(index);
	[~, order_ind] = sort(CFs);
	if ibin == 1
		CFs_ordered_contra = CFs(order_ind);
		MTFs_contra = MTF_types(index(order_ind));
		RMs_contra = RM_types(index(order_ind));
	else
		CFs_ordered = CFs(order_ind);
		MTFs = MTF_types(index(order_ind));
		RMs = RM_types(index(order_ind));
	end
	
	for iNo = 1:3 %1:num_Nos
		for isesh = 1:num_sesh
			s_ind = index(order_ind(isesh));

			% Load in session
			putative_neuron = sessions.Putative_Units{s_ind};
			load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');
			CF = CFs(order_ind(isesh));
			MTF_type = MTF_types{index(order_ind(isesh))};
			RM_type = RM_types{index(order_ind(isesh))};

			if CF < 2000
				igroup = 1;
			elseif CF >= 2000 && CF < 4000
				igroup = 2;
			else
				igroup = 3;
			end

			% Get data for each stimulus
			if isempty(data{5+iNo,ibin})

			else
				params_WB = data(5+iNo,ibin); % Gets binaural WB-TIN stimuli

				% General analysis
				data_WB = analyzeWBTIN(params_WB, []);
				data_WB = data_WB{1};
				noise_alone = mean(data_WB.rate(1,:));

				for iSNR = 1:3
					if length(data_WB.SNRs) == 3
						WB = data_WB.rate(:,iSNR); % 30 dB SNR
						V_p(iNo, isesh) = data_WB.V_p(iSNR); % 30 dB SNR
					elseif length(data_WB.SNRs) == 2
						if iSNR == 1
							continue
						else
							WB = data_WB.rate(:,iSNR-1); % 30 dB SNR
							V_p(iNo, isesh) = data_WB.V_p(iSNR-1); % 30 dB SNR
						end
					end
					fpeaks = data_WB.fpeaks(:,end);
					fpeaks_re_CF = log2(fpeaks/CF);

					% Align by CF (approximately)
					f = linspace(-3.5, 3.5, 10000);
					[~, f_ind(1)] = min(abs(fpeaks_re_CF(2)-f));
					[~, f_ind(2)] = min(abs(fpeaks_re_CF(end)-f)); % find indices
					f_interp = linspace(f(f_ind(1)),f(f_ind(2)), f_ind(2)-f_ind(1));

					% Interpolate & get z-score
					r_interp = interp1(fpeaks_re_CF(2:end), WB(2:end),f_interp, 'spline');
					z_rate = (r_interp - noise_alone) / std(r_interp);
					if ibin == 2
						WB_mat(iNo, iSNR, isesh, f_ind(1):f_ind(2)-1) = z_rate;
					else
						WB_mat_contra(iNo, iSNR, isesh, f_ind(1):f_ind(2)-1) = z_rate;
					end

					% Create giant table of everything
					tone_freq_names = -1:0.25:1;
					names = strtrim(cellstr(num2str(tone_freq_names'))');
					tone_freqs = 1:-0.25:-1;
					for itone = 1:9
						[~,ind] = min(abs(f_interp+tone_freqs(itone)));
						analysisTable.Putative{iind} = putative_neuron;
						analysisTable.MTF{iind} = MTF_type;
						analysisTable.CF(iind) = CF;
						analysisTable.CF_Group{iind} = CF_groups{igroup};
						analysisTable.RM{iind} = RM_type;
						analysisTable.No(iind) = Nos(iNo);
						analysisTable.SNR(iind) = SNRs(iSNR);
						analysisTable.Binmode{iind} = binmodes{ibin};
						analysisTable.Tone_Freq(iind) = names{itone};
						analysisTable.Rate(iind) = z_rate(ind);
						iind = iind + 1;
					end
				end
			end
			fprintf('%s, No = %d, %d percent done \n', binmodes{ibin},Nos(iNo), round(isesh/num_sesh*100))

		end
	end
end

%% Get rid of RM types with really small sample sizes 
% 
% types = unique(analysisTable.RM);
% for ii = 1:length(types)
% 	num_type(ii) = sum(strcmp(types{ii}, analysisTable.RM));
% end
% 
% % Problems are with inhibitory, O, off, on, unusual, on/off
% for ii = 2:7
% 	ind_inhibitory = strcmp(types{ii}, analysisTable.RM);
% 	analysisTable(ind_inhibitory,:) = [];
% end

%% Save spreadsheet 

filename = 'Population_Data.xlsx';
writetable(analysisTable, fullfile(datapath, filename));

elapsedTime = toc(timerVal)/60;
disp(['This took ' num2str(elapsedTime) ' minutes'])

%% Save data 

save(fullfile(datapath, 'Population_Data.mat'), "WB_mat_contra", "WB_mat", "f",...
	 "CFs_ordered", "CFs_ordered_contra", "MTFs_contra", "RMs_contra", ...
	"RMs", "MTFs")
