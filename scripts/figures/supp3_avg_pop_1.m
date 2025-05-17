function supp3_avg_pop_1(save_fig)
% supp_avg_pop_1 plots Figure S3.
%	This function plots average response and +/- 1 SEM for contralateral
%	WB-TIN responses for 0, 20 and 40 dB SNR and 3, 23, 43 dB SPL. 
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN, 
%	 .mat files required: Population_Data.mat
%	 spreadsheets required: -
%
% Author: J. Fritzinger
% Created: 2025-02-18; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------


%% Load in data

[~, datapath, ~, ppi] = get_paths();
load(fullfile(datapath, 'Population_Data.mat'), 'WB_mat_contra', 'f')

%% Set up figure

figure('Position',[50,50,6*ppi,2*ppi])
tiledlayout(1, 3,'Padding','loose', 'TileSpacing','compact')
legsize = 7;
titlesize = 9;
fontsize = 8;
labelsize = 13;
CF_tests = -1:0.25:1;
colors = {'b', 'r', 'g'};

%% Plots 

SNRs = {'20', '30', '40'};
for iSNR = 1:3
	nexttile
	hold on
	for ii = 1:9
		xline(CF_tests(ii), '--', 'color', [0.5 0.5 0.5])
	end
	for iNo = 1:3
		MTF_ind = 1:size(WB_mat_contra,3);
		WB_MTF_mat = WB_mat_contra(:,iSNR, MTF_ind, :);
		WB_array_z = squeeze(WB_MTF_mat(iNo, 1, :,:));

		oct_mat = WB_array_z;
		num_samples = sum(~isnan(oct_mat),1, 'omitnan');
		oct_avg_all = mean(oct_mat, 1, 'omitnan');
		oct_std_all = std(oct_mat, 1, 'omitnan')./sqrt(num_samples);

		not_nan_ind = ~isnan(oct_avg_all);
		oct_avg = oct_avg_all(not_nan_ind);
		oct_std = oct_std_all(not_nan_ind);
		f_notnan = f(not_nan_ind);


		plot(f_notnan, oct_avg, 'color', 'k', 'LineWidth',1.5);
		patch([f_notnan flip(f_notnan)], [oct_avg-oct_std flip(oct_avg+oct_std)], ...
			colors{iNo}, 'FaceAlpha',0.25, 'EdgeColor','none')
		yline(0)
		xline(0)
		xlim([-1.5 1.5])
		ylim([-1 1.5])
		set(gca, 'FontSize', fontsize)
		title([SNRs{iSNR} ' dB SNR'], 'fontsize', titlesize)

		if iSNR == 2
			xlabel('Tone Frequency w.r.t. CF (octaves)')
		end
		if iSNR == 1
			ylabel('Noise-ref. z-score')
		end
	end
end
leg1 = {'', '', '', '', '', '', '', '', '', '', '3 dB N_0', '','',...
	'','23 dB N_0', '','','','43 dB N_0'};
leg = legend(leg1, 'fontsize', legsize);
leg.ItemTokenSize = [10,8];

%% Set positions and labels

% Annotate 'Contra';
annotation('textbox',[0.05 0.3 0.10 0.053],'String',{'Contralateral'},...
	'Rotation',90,'FontWeight','bold','FontSize',titlesize+2,'EdgeColor','none');

% Figure labels
height = [0.94 0.46];
left = [0.05 0.35 0.63];
labels = {'A', 'B', 'C'};
row = reshape(repmat(1:2, 3, 1), 6, 1);
col = repmat(1:3, 1, 2);
for ii = 1:3
	irow = row(ii);
	icol = col(ii);
	annotation('textbox',[left(icol) height(irow) 0.026 0.065],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'FitBoxToText','off','EdgeColor','none');
end

%% Export 

if save_fig == 1
	saveFigures('Supp3_avg_pop_1')
end
end