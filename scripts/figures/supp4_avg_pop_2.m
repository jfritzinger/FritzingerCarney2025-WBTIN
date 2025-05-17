function supp4_avg_pop_2(save_fig)
% supp_avg_pop_2 plots Figure S4.
%	This function plots V-type and MTF type population averages, combining
%	contra and binaural responses. 
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN, compareMultiple
%	 .mat files required: Population_Data.mat
%	 spreadsheets required: -
%
% Author: J. Fritzinger
% Created: 2025-02-18; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------

%% Load in data

% Load in data
[~, datapath, ~, ppi] = get_paths();
load(fullfile(datapath, 'Population_Data.mat'), 'WB_mat', ...
	'RMs', 'MTFs', 'WB_mat_contra', 'RMs_contra', 'MTFs_contra', 'f')

%% Put contra and binaural together!

WB_mat = cat(3, WB_mat, WB_mat_contra);
RMs = cat(1, RMs, RMs_contra);
MTFs = cat(1, MTFs, MTFs_contra);

names_MTF = {'BE', 'BS', 'F', 'H_B_E', 'H_B_S'};
names_RM = {'I-type', 'V-type'};
ind_BE = strcmp('BE', MTFs);
ind_BS = strcmp('BS', MTFs);
ind_F = strcmp('F', MTFs);
ind_HBE = strcmp('H_BE', MTFs);
ind_HBS = strcmp('H_BS', MTFs);
ind_I = strcmp('I', RMs);
ind_V = strcmp('V', RMs);
SNR = 3;

%% Set up figure

figure('Position',[100,100,5*ppi,3*ppi])
tiledlayout(1, 2,'Padding','compact', 'TileSpacing','compact')
legsize = 7;
titlesize = 9;
fontsize = 8;
labelsize = 13;
loc = 'northwest';

%% Create figure 

for itype = 1:2
	if itype == 1
		ind_temp = ind_I;
	else
		ind_temp = ind_V;
	end
	WB_temp1 = WB_mat(:, SNR, ind_BS&ind_temp, :);
	WB_temp2 = WB_mat(:, SNR, ind_BE&ind_temp, :);
	WB_temp3 = WB_mat(:, SNR, ind_F&ind_temp, :);
	WB_temp4 = WB_mat(:, SNR, ind_HBE&ind_temp, :);
	WB_temp5 = WB_mat(:, SNR, ind_HBS&ind_temp, :);

	nexttile
	compareMultiple(names_MTF, names_RM{itype}, f, legsize, loc, ...
		WB_temp1, WB_temp2, WB_temp3, WB_temp4, WB_temp5);
	if itype == 1
		ylabel('Noise ref. z-score')
	end
	xlabel('Tone Freq.w.r.t. CF (oct)')
	set(gca, 'fontsize', fontsize)
	if itype == 1
		title('I-type', 'fontsize', titlesize)
	else
		title('V-type', 'fontsize', titlesize)
	end
end

%% Set positions and labels

% Figure labels
height = 0.94;
left = [0.01 0.52];
labels = {'A', 'B'};
row = reshape(repmat(1:2, 2, 1), 4, 1);
col = repmat(1:2, 1, 2);
for ii = 1:2
	irow = row(ii);
	icol = col(ii);
	annotation('textbox',[left(icol) height(irow) 0.0259 0.066],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'FitBoxToText','off','EdgeColor','none');
end

%% Export

if save_fig == 1
	saveFigure('Supp4_avg_pop_2')
end
end
