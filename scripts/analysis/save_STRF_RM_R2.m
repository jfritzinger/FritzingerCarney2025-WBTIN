%% save_STRF_RM_R2.m


%% Load in data

[~, datapath, ~, ~] = get_paths();
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), ...
	'PreserveVariableNames',true);
num_units = size(sessions, 1);

load(fullfile(datapath, 'STRFModel_23.mat'), 'STRF')
%load(fullfile(datapath, 'STRFModel_23_contra.mat'), 'STRF_contra')
%load(fullfile(datapath, "Fig12_RM_contra.mat"), "RM_contra")

%% Sessions of interest

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
STRF_units = {STRF.putative}';

binmodes = {'Contra', 'Binaural'};
Nos = {'3', '23', '43'};
ind = 1;
for ibin = 2

	if ibin == 1
		index = find(contra);
		R2_STRF_contra = NaN(length(index), 3);
	else
		index = find(bin);
		R2_STRF_bin = NaN(length(index), 3);
	end

	RM_R2 = NaN(length(index), 3);
	num_sesh = length(index);
	%MTF_type = zeros(num_sesh,1);
	for isesh = 1:num_sesh
		s_ind = index(isesh);

		% Load in session
		putative_neuron = sessions.Putative_Units{s_ind};
		CF = sessions.CF(s_ind);
		RM_type_one = sessions.RM_Type(s_ind);

		if strcmp(RM_type_one{1}, 'V') || strcmp(RM_type_one{1}, 'I') ...
				|| strcmp(RM_type_one{1}, 'On')

			RM_type(isesh) = RM_type_one;
			load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');

			% Get data for each stimulus
			params_RM = data{2, 2};
			params_STRF = data{4, ibin};

			if isempty(params_STRF)
				disp(putative_neuron)
			end

			params_WB = data(7,ibin); % Gets binaural WB-TIN stimuli
			params_MTF = data{3, ibin};

			if isempty(params_MTF)&& ibin == 1
				params_MTF = data{3,2};
			elseif isempty(params_MTF)&& ibin == 2
				params_MTF = data{3,1};
			end

			if isempty(params_RM)
				params_RM = data{2,1};
			end

			% General analysis
			data_RM = analyzeRM(params_RM);
			data_WB = analyzeWBTIN(params_WB, []);
			params_WB = params_WB{1};
			data_WB = data_WB{1};

			% RM Correlation
			if isempty(data_WB.fpeaks)
				RM_R2(isesh,:) = NaN;
			else

				% Find BF nearest to WBTIN overall level
				lbw = params_WB.fpeak_mid * 2^(-params_WB.bandwidth/2);
				hbw = params_WB.fpeak_mid * 2^(params_WB.bandwidth/2);
				overall_spl = params_WB.No + 10*log10(hbw-lbw);
				[~,iBF] = min(abs(data_RM.SPL-overall_spl));
				BF_spl = data_RM.SPL(iBF);

				num_snr = length(params_WB.SNR);
				for isnr = num_snr

					% Analysis
					f_interp = logspace(log10(data_WB.fpeaks(2,1)), log10(data_WB.fpeaks(end,1)), 50);
					r_WBTIN = data_WB.rate(2:end,isnr);
					r_interp_WBTIN = interp1(data_WB.fpeaks(2:end,1),r_WBTIN,f_interp,'pchip'); % interpolates rate

					r_RM = data_RM.rates(:,iBF);
					f_mid = logspace(log10(data_RM.freqs(1)), log10(data_RM.freqs(end)), 700);
					r_mid_RM = interp1(data_RM.freqs,r_RM,f_mid,'pchip'); % interpolates rate
					[~, starting] = min(abs(f_mid-f_interp(1)));
					[~, ending] = min(abs(f_mid-f_interp(end)));
					r_interp_RM = interp1(f_mid(starting:ending),r_mid_RM(starting:ending),f_interp,'pchip'); % interpolates rate

					R_int = corrcoef(r_interp_RM,r_interp_WBTIN);
					r2_RM = R_int(1, 2).^2;

					if num_snr == 3
						RM_R2(isesh, isnr) = r2_RM;
					elseif num_snr == 2
						RM_R2(isesh, isnr+1) = r2_RM;
					end
				end
			end

			% STRF data
			if ibin == 2
				R2_temp = STRF(isesh).R2;
				if length(R2_temp)==3
					R2_STRF_bin(isesh, :) = R2_temp;
				elseif length(R2_temp)==2
					R2_STRF_bin(isesh, 2:3) = R2_temp;
				end
			else
				R2_temp = STRF_contra(isesh).R2;
				if length(R2_temp)==3
					R2_STRF_contra(isesh, :) = R2_temp;
				elseif length(R2_temp)==2
					R2_STRF_contra(isesh, 2:3) = R2_temp;
				end
			end
		end
	end

	% Add to matrix for analysis
	if ibin == 1
		R2_RM_contra = RM_R2;
		RM_type_contra = RM_type;
	else
		R2_RM_bin = RM_R2;
		RM_type_bin = RM_type;
	end
	clear MTF_type
end

%% Save

filename = 'STRF_RM_R2.mat';
save(fullfile(datapath, filename), "R2_RM_bin", "RM_type_bin", ...
	"R2_STRF_bin");