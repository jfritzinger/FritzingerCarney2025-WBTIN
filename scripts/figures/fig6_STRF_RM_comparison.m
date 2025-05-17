function fig6_STRF_RM_comparison(save_fig)
% STRF_RM_comparison plots Figure 5.
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN.m, analyzeRM, analyzeWBTIN
%	 .mat files required: STRFModel_23.mat, STRF_RM_R2.mat, Neural_Data
%	 spreadsheets required: Data_Table.xlsx
%
% Author: J. Fritzinger
% Created: 2025-01-30; Last revision: 2025-03-06
%
% -------------------------------------------------------------------------


%% Load in data

% Load in spreadsheet
[~, datapath, ~, ppi] = get_paths();
sessions = readtable(fullfile(datapath, 'Data_Table.xlsx'), ...
	'PreserveVariableNames',true);

% Load in data
load(fullfile(datapath, 'STRFModel_23.mat'), 'STRF')
load(fullfile(datapath, 'STRF_RM_R2.mat'), "R2_RM_bin", "RM_type_bin", ...
	"R2_STRF_bin");

%% Set up figure

figure('Position',[50,50,3.346*ppi,5.5*ppi]);
legsize = 7;
fontsize = 7;
linewidth = 1;
labelsize = 13;
RM_colors = {'k', 'k', '#A49BD0', '#5E50A9', '#20116B'};
WB_colors = {'#3F985C', '#3F985C', '#3F985C'};

%% RM analysis

h = gobjects(9, 1);
for ineuron = 1:2
	if ineuron == 1 % V-type
		putative_neuron = 'R25_TT2_P8_N14';
	else % I-type
		putative_neuron = 'R27_TT4_P7_N16';
	end

	% Load in session
	load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');
	indx = find(strcmp(sessions.Putative_Units, putative_neuron));
	CF = sessions.CF(indx);

	% Get data for each stimulus
	params_RM = data{2, 2};
	params_WB = data(7,2); % Gets binaural WB-TIN stimuli

	% General analysis
	data_RM = analyzeRM(params_RM);
	datas_WB = analyzeWBTIN(params_WB, []);

	% Find peak rates
	if ineuron == 1
		max_rate = 30;
	else
		max_rate = 41;
	end

	iNo = 1;
	param_WB = params_WB{iNo};
	data_WB = datas_WB{iNo};

	% Find BF nearest to WBTIN overall level
	lbw = param_WB.fpeak_mid * 2^(-param_WB.bandwidth/2);
	hbw = param_WB.fpeak_mid * 2^(param_WB.bandwidth/2);
	overall_spl = param_WB.No + 10*log10(hbw-lbw);
	[~,iBF] = min(abs(data_RM.SPL-overall_spl));
	BF_spl = data_RM.SPL(iBF);

	%
	freqs_interp = logspace(log10(data_WB.fpeaks(2,1)), ...
		log10(data_WB.fpeaks(end,1)), 50);
	rate_WBTIN = data_WB.rate(2:end,end);
	rate_interp_WBTIN = interp1(data_WB.fpeaks(2:end,1),rate_WBTIN,...
		freqs_interp,'pchip'); % interpolates rate

	rate_RM = data_RM.rates(:,iBF);
	freqs_mid = logspace(log10(data_RM.freqs(1)), log10(data_RM.freqs(end)), 700);
	rate_mid_RM = interp1(data_RM.freqs,rate_RM,freqs_mid,'pchip'); 
	[~, starting] = min(abs(freqs_mid-freqs_interp(1)));
	[~, ending] = min(abs(freqs_mid-freqs_interp(end)));
	rate_interp_RM = interp1(freqs_mid(starting:ending),...
		rate_mid_RM(starting:ending),freqs_interp,'pchip'); % interpolates rate

	R_int = corrcoef(rate_interp_RM,rate_interp_WBTIN);
	r2 = R_int(1, 2).^2;

	% Plot
	% Plot RM/WBTIN
	plot_ind = [1 2];
	h(plot_ind(ineuron)) = subplot(4, 3, plot_ind(ineuron));
	hold on
	line([data_WB.fpeaks(2) data_WB.fpeaks(end)], [1 1]* data_WB.rate(1,end),...
		'Color','#034E1C', 'linewidth', linewidth);
	line([data_RM.freqs(1) data_RM.freqs(end)], [1 1]*data_RM.spont,...
		'Color', '#20116B', 'linewidth', linewidth);
	area(data_RM.freqs, data_RM.rates(:,iBF), 'EdgeColor',RM_colors{iBF}, ...
		'FaceColor',RM_colors{iBF}, 'FaceAlpha',0.5);
	errorbar(data_WB.fpeaks(:,1),data_WB.rate(:,end),data_WB.rate_std(:,end)...
		/(sqrt(param_WB.nrep)), 'Color',WB_colors{iNo}, 'linewidth', linewidth);
	xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',linewidth)

	set(gca,'FontSize',fontsize)
	grid on
	set(gca, 'XScale', 'log');
	xlim([data_WB.fpeaks(2) data_WB.fpeaks(end)])
	ylim([0 max_rate+5])

	xL=xlim;
	yL=ylim;
	message = ['R^2 = ' num2str(round(r2, 2))];
	text(1.05*xL(1),1*yL(2),message,'HorizontalAlignment','left',...
		'VerticalAlignment','top', 'FontSize',legsize)
	xticklabels([])
	hLegend = legend('', '', [num2str(BF_spl) ' dB SPL'], [...
		num2str(round(param_WB.No_overall)) ' dB SPL'], 'fontsize',...
		legsize, 'box', 'off');
	hLegend.ItemTokenSize = [6,6];
	ylabel('Avg. Rate (sp/s)')
	xlabel('Tone Freq. (kHz)')
	xticklabels(xticks/1000)
end



%% Plot

for ineuron = 1:2
	switch ineuron
		case 1 % Good
			putative_neuron = 'R24_TT2_P13_N2';
		case 2 % Bad
			putative_neuron = 'R29_TT2_P3_N2';
	end

	% Load in session
	load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');
	load(fullfile(datapath, 'Neural_Data', [putative_neuron '_STRF.mat']), 'data_STRF')
	indx = find(strcmp(sessions.Putative_Units, putative_neuron));
	CF = sessions.CF(indx);

	% Get data for each stimulus
	params_WB = data(7,2); % Gets binaural WB-TIN stimuli

	% General analysis
	data_WB = analyzeWBTIN(params_WB, []);
	params_WB = params_WB{1};
	data_WB = data_WB{1};

	% Get STRF model results
	model_STRF = STRF(strcmp({STRF.putative}, putative_neuron));

	% Plot STRF
	plot_ind = [3 5];
	h(plot_ind(ineuron)) = subplot(4, 3, plot_ind(ineuron));
	STRF_mat = data_STRF.H2ex_strf-data_STRF.H2in_strf;
	imagesc(data_STRF.t.*1000, data_STRF.f./1000, STRF_mat, data_STRF.clims_strf);
	set(gca,'Ydir','normal','XLim',data_STRF.tlims.*1000,'YLim',...
		[params_WB.fpeaks(2) params_WB.fpeaks(end)]./1000)
	hold on
	yline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',2)
	colormap(redblue)
	grid on
	set(gca,'FontSize',fontsize)
	xlabel('Time (ms)');
	ylabel('Freq. (kHz)')
	box off

	% Plot real response
	plot_ind = [4 6];
	h(plot_ind(ineuron)) = subplot(4, 3, plot_ind(ineuron));
	hold on
	noise_alone = mean(data_WB.rate(1,:));
	num_SNRs = length(params_WB.SNR);
	px = [params_WB.fpeaks(2) params_WB.fpeaks(end) params_WB.fpeaks(end) ...
		params_WB.fpeaks(2)];
	rate_pos = noise_alone+mean(data_WB.rate_std(1,1:num_SNRs))/...
		(sqrt(params_WB.nrep));
	rate_neg = noise_alone-mean(data_WB.rate_std(1,1:num_SNRs))/...
		(sqrt(params_WB.nrep));
	py = [rate_neg rate_neg rate_pos rate_pos];
	patch(px, py, [0.9 0.9 0.9], 'EdgeColor', [0.8 0.8 0.8])
	line([data_WB.fpeaks(2) data_WB.fpeaks(end)], [1 1]*noise_alone,...
		'Color','k', 'linewidth', linewidth);
	errorbar(data_WB.fpeaks(:,end),data_WB.rate(:,end),...
		data_WB.rate_std(:,end)/(sqrt(params_WB.nrep)), 'LineWidth',...
		linewidth, 'Color', WB_colors{3});

	% Plot STRF prediction
	RM_colors = {'#8856a7', '#A49BD0', '#5E50A9', '#20116B'};
	hold on
	noise_alone_m = mean(model_STRF(1).avModel(1,:));
	line([params_WB.fpeaks(2) params_WB.fpeaks(end)], [1 1]*noise_alone,...
		'Color','k', 'LineWidth',linewidth);
	ratio =  max(data_WB.rate(:,end)-noise_alone)/...
		max(model_STRF(1).avModel(:,end)-noise_alone_m);
	plot(params_WB.fpeaks,((model_STRF(1).avModel(:,end)-...
		noise_alone_m).*ratio)+noise_alone, 'LineWidth',linewidth, ...
		'Color',RM_colors{1});
	xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',linewidth)
	hold off
	xlabel('Tone Freq. (kHz)')
	ylabel('Avg. Rate (sp/s)')
	set(gca, 'XScale', 'log');
	grid on
	xlim([params_WB.fpeaks(2) params_WB.fpeaks(end)]);
	set(gca,'FontSize',fontsize)
	ylim([0 model_STRF(1).max_all+7])
	if ineuron == 2
		hLegend = legend(char('', 'Noise Alone', 'WB-TIN', '', 'STRF Model'), ...
			'location', 'southwest', 'Fontsize', legsize-2, 'box', 'off');
		hLegend.ItemTokenSize = [6,6];
	end
	xticklabels(xticks/1000)
	fprintf('R^2 = %0.2f\n', model_STRF(1).R2(end))
end

%% Scatter Plot

isnr = 3;
r2_STRF = R2_STRF_bin(:,isnr);
r2_RM = R2_RM_bin(:,isnr);
RM_type = RM_type_bin;

empty_units = cellfun('isempty', RM_type);
RM_type = RM_type(~empty_units);
r2_STRF = r2_STRF(~empty_units);
r2_RM = r2_RM(~empty_units);

h(7) = subplot(4, 3, 7); % [1 4 7 10 13 16]
h(8) = subplot(4, 3, 8);
h(9) = subplot(4, 3, 9);

% Scatter plot
a = 7;
hold(h(a), "on")
for itype = 1:3
	if itype == 1
		subtype_i = cellfun(@(a) strcmp(a, 'V'), RM_type);
	elseif itype == 2
		subtype_i = cellfun(@(a) strcmp(a, 'I'), RM_type);
	else
		subtype_i = cellfun(@(a) strcmp(a, 'On'), RM_type);
	end
	r2_STRF_subset = r2_STRF(subtype_i);
	r2_RM_subset = r2_RM(subtype_i);
	scatter(h(a), r2_STRF_subset,r2_RM_subset, 18, 'filled', 'MarkerEdgeColor','k')
end
plot(h(a),[0 median(r2_STRF, 'omitnan')], [median(r2_RM, 'omitnan') ...
	median(r2_RM, 'omitnan')],'color', [0.4 0.4 0.4], 'LineWidth',2)
plot(h(a),[median(r2_STRF, 'omitnan') median(r2_STRF, 'omitnan')], ...
	[0 median(r2_RM, 'omitnan')],'color', [0.4 0.4 0.4], 'LineWidth',2)
set(h(a), 'FontSize', fontsize, 'YLabel', [], 'XLabel', [])
xlim(h(a),[0 1])
ylim(h(a),[0 1])
xticks(h(a), 0:0.2:1)
yticks(h(a), 0:0.2:1)
xticklabels(h(a), [])
yticklabels(h(a), [])
grid on
hLegend = legend(h(a), 'V-type', 'I-type', 'Onset', '','');
hLegend.ItemTokenSize = [4,4];

% Create distribution plot for the X-axis (horizontal)
b = 8;
edges = 0:0.05:1;
hold(h(b), "on")
xline(h(b), 0, 'k')
histogram(h(b), r2_STRF,edges, 'Orientation', 'vertical');
xlabel(h(b),'STRF Model R^2')
set(h(b), 'fontsize', fontsize)
xlim(h(b), [0 1])
hold(h(b), 'on')
xline(h(b), median(r2_STRF, 'omitnan'),'color', [0.4 0.4 0.4], 'LineWidth',2)
grid on
yticks(h(b), [0 40])

% Create distribution plot for the Y-axis (vertical) below the scatter plot
c = 9;
hold(h(c), "on")
yline(h(c), 0, 'k')
histogram(h(c), r2_RM, edges,  'Orientation', 'horizontal');
set(h(c), 'XDir', 'reverse');
ylabel(h(c),'RM R^2')
set(h(c), 'fontsize', fontsize)
ylim(h(c), [0 1])
xticks(h(c), [0 20])
hold(h(c), 'on')
yline(h(c), median(r2_RM, 'omitnan'),'color', [0.4 0.4 0.4], 'LineWidth',2)
grid on


%% Position plots
set(gcf, 'Position',[50,50,3.346*ppi,6*ppi])

left = [0.13 0.63];
bottom = [0.835 0.625 0.42 0.05];
width = 0.34;
height = 0.13;

% [left, bottom, width, height]
set(h(1), 'Position', [left(1) bottom(1) width height]);
set(h(2), 'Position', [left(2) bottom(1) width height]);

set(h(3), 'Position', [left(1) bottom(2) width height]);
set(h(4), 'Position', [left(2) bottom(2) width height]);

set(h(5), 'Position', [left(1) bottom(3) width height]);
set(h(6), 'Position', [left(2) bottom(3) width height]);

set(h(7), 'Position', [left(1)+0.15 bottom(4)+0.06 0.69 0.24]);
set(h(8), 'Position', [left(1)+0.15 bottom(4)      0.69 0.04]);
set(h(9), 'Position', [left(1)+0.02 bottom(4)+0.06 0.08 0.24]);

% Add labels
labels = {'A', 'B', 'C', 'D'};
annotation('textbox',[0.01 0.955 0.077 0.054],'String',labels{1},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.74 0.077 0.054],'String',labels{2},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.535 0.077 0.054],'String',labels{3},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.01 0.335 0.077 0.054],'String',labels{4},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Export

if save_fig == 1
	saveFigures('Fig6_STRF_RM_comparison')
end
end