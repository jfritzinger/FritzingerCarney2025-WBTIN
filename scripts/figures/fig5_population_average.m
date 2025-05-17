function fig5_population_average(save_fig)
% population_average plots Figure 4. 
%	This function plots population average plots for binaural condition,
%	WB-TIN. It plots SNR = 20, 30, 40 dB and N0 = 3, 23, and 43 dB SPL. It
%	also plots average I- and V-type responses split into low, medium, and
%	high CFs. 
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN.m, compareMultiple.m
%	 .mat files required: 'Population_Data.mat'
%	 spreadsheets required: -
%
% Author: J. Fritzinger
% Created: 2025-01-22; Last revision: 2025-03-06
%
% -------------------------------------------------------------------------


%% Load in data
[~, datapath, ~, ppi] = get_paths();

% Load in data
name = 'Population_Data.mat';
load(fullfile(datapath, name), 'CFs_ordered', 'RMs', 'WB_mat', ...
	'WB_mat_contra', 'RMs_contra', 'CFs_ordered_contra', 'f');

%% Set up figure

loc = 'northeastoutside';
h = gobjects(5, 1);
position1 = [0.852,0.757,0.13354,0.1388];

figure('Position',[50,50,4.567*ppi,2.3*ppi]);
legsize = 5;
titlesize = 8;
fontsize = 7;
linewidth = 1;
labelsize = 13;

%% Plot 

CF_tests = -1:0.25:1;
colors = {'b', 'r', 'g'};
for iSNR = 1:3
	h(iSNR) = subplot(2, 3, iSNR);
	hold on
	for ii = 1:9
		xline(CF_tests(ii), '--', 'color', [0.5 0.5 0.5])
	end
	for iNo = 1:3
		MTF_ind = 1:size(WB_mat,3);
		WB_MTF_mat = WB_mat(:,iSNR, MTF_ind, :);
		WB_array_z = squeeze(WB_MTF_mat(iNo, 1, :,:));

		oct_mat = WB_array_z;
		num_samples = sum(~isnan(oct_mat),1, 'omitnan');
		oct_avg_all = mean(oct_mat, 1, 'omitnan');
		oct_std_all = std(oct_mat, 1, 'omitnan')./sqrt(num_samples);

		not_nan_ind = ~isnan(oct_avg_all);
		oct_avg = oct_avg_all(not_nan_ind);
		oct_std = oct_std_all(not_nan_ind);
		f_notnan = f(not_nan_ind);

		plot(f_notnan, oct_avg, 'color', 'k', 'LineWidth',linewidth);
		patch([f_notnan flip(f_notnan)], [oct_avg-oct_std flip(oct_avg+oct_std)], ...
			colors{iNo}, 'FaceAlpha',0.25, 'EdgeColor','none')
	end

	% Label plot
	yline(0)
	xline(0)
	xlim([-1.5 1.5])
	ylim([-1 1.5])
	set(gca, 'FontSize', fontsize)
	grid on
	if iSNR == 1
		ylabel('Noise-ref. z-score')
		title('20 dB SNR')
	elseif iSNR == 2
		xlabel('Tone Freq. w.r.t. CF (Oct.)')
		title('30 dB SNR')
		yticklabels([])
	else
		title('40 dB SNR')
		yticklabels([])
		leg2 = {'', '', '', '', '', '', '', '', '','', '3 dB, n=207', '',...
			'23 dB, n=207','','43 dB, n=207'};
		leg = legend(leg2, 'fontsize', legsize, 'Position',position1);
		leg.ItemTokenSize = [5,5];
	end
end





%% Second figure: RMs on different plots 

% Put contra and binaural together! 
WB_mat = cat(3, WB_mat, WB_mat_contra);
CFs_ordered = cat(1, CFs_ordered, CFs_ordered_contra);
RMs = cat(1, RMs, RMs_contra);

names_RM = {'I-type', 'V-type'};
names_CF = {'Low CF', 'Med CF', 'High CF'};
ind_I = strcmp('I', RMs);
ind_V = strcmp('V', RMs);
ind_low = CFs_ordered<2000;
ind_med = CFs_ordered>=2000&CFs_ordered<4000;
ind_high = CFs_ordered>=4000;
SNR = 3;

% CF groups plotted in 2 RM plots 
for itype = 1:2
	if itype == 1
		ind_temp = ind_I;
	else
		ind_temp = ind_V;
	end

	WB_temp1 = WB_mat(:, SNR, ind_low&ind_temp, :);
	WB_temp2 = WB_mat(:, SNR, ind_med&ind_temp, :);
	WB_temp3 = WB_mat(:, SNR, ind_high&ind_temp, :);

	%nexttile
	h(3+itype) = subplot(2, 3, itype+3);
	compareMultiple(names_CF, names_RM{itype}, f, legsize, loc, ...
		WB_temp1, WB_temp2, WB_temp3);
	set(gca, 'FontSize', fontsize)

	title(names_RM{itype}, 'fontsize', titlesize)
	if itype == 1
		ylabel('Noise ref. z-score')
	end
	xlabel('Tone Freq. w.r.t. CF (Oct.)')
end

%% Set positions and labels

% Set positions
height = 0.31;
width = 0.25;
left = [linspace(0.13, 0.67, 3) 0.12 0.58];
bottom = [repmat(0.62, 3, 1) repmat(0.12,3, 1) ];
for ii = 1:5
	set(h(ii), 'position', [left(ii) bottom(ii) width height])
end

% Figure labels
height = [0.97 0.45];
left = [0.0 0.0 0.5];
labels = {'A', 'B', 'C', 'D' 'E'};
row = [1 2 2];
for ii = 1:3
	irow = row(ii);
	annotation('textbox',[left(ii) height(irow) 0.0259 0.066],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'FitBoxToText','off','EdgeColor','none');
end

% Annotate 'Diotic' and 'Contra'
annotation('textbox',[0.05 0.15 0.15 0.07],'String','RM Type',...
	'Rotation',90,'FontWeight','bold','FontSize',titlesize,...
	'EdgeColor','none');
annotation('textbox',[0.04 0.68 0.1 0.053],'String','Diotic',...
	'Rotation',90,'FontWeight','bold','FontSize',titlesize,...
	'EdgeColor','none');

%% Export

if save_fig == 1
	saveFigures('Fig5_population_average')
end
end


