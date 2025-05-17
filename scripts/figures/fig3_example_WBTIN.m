function fig3_example_WBTIN(save_fig)
% example_WBTIN plots Figure 3.
%	This function plots an example neuron response to response map,
%	modulation transfer function, STRF, on-CF NB and WB TIN, and three
%	levels of WB-TIN
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN, analyzeRM, analyzeMTF, analyzeWBTIN,
%	 analyzeSTRF, analyzeNBTIN
%	 .mat files required: Neural_Data
%	 spreadsheets required: Data_Table.xlsx
%
% Author: J. Fritzinger
% Created: 2024-05-29; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------

%% Example Session

% R29_TT4_P2_N11
putative_neuron = 'R29_TT4_P2_N11';
y_ticks = [0 20 40 60];
y_range = [0 64];
MTF_range = [0 70];
x_label = [1000 2000 3000 4000 6000 8000];

%% Load in spreadsheet

[~, datapath, ~, ppi] = get_paths();

% Load in spreadsheet & data
spreadsheet_name = 'Data_Table.xlsx';
sessions = readtable(fullfile(datapath, spreadsheet_name), ...
	'PreserveVariableNames',true);
load(fullfile(datapath, 'Neural_Data', [putative_neuron '.mat']), 'data');

% Find example in spreadsheet
s_ind = strcmp(sessions.Putative_Units, putative_neuron);
CF = sessions.CF(s_ind);

% Get data for each stimulus
params_RM = data{2,2}; % Gets binaural response map
params_MTF = data{3,2}; % Gets binaural MTFN
params_WB = data(6:8,2); % Gets binaural WB-TIN stimuli
params_STRF = data{4,2};
params_NB = data(9:11,2);
plot_range = [params_WB{1}.fpeaks(2) params_WB{1}.fpeaks(end)];


%% Figure Parameters

figure('Position', [50,50,6.929*ppi,3.2*ppi]);
CF_color = [0.3 0.3 0.3];
legsize = 6;
titlesize = 8;
fontsize = 7;
spont_color = [0.4 0.4 0.4];
linewidth = 1;
labelsize = 13;

%% Plot RM

% Analysis
data_RM = analyzeRM(params_RM);

% Plot
h = gobjects(10, 1);
for ind = 1:3
	h(ind) = subplot(3, 4, ind);
	hold on

	if ind == 1
		area(data_RM.freqs,data_RM.rates(:,3),'EdgeColor', '#A49BD0',...
			'FaceColor', '#A49BD0','LineWidth',linewidth) % 33 dB
		spl = '33 dB SPL';
	elseif ind == 2
		area(data_RM.freqs,data_RM.rates(:,4),'EdgeColor', '#5E50A9',...
			'FaceColor', '#5E50A9','LineWidth',linewidth) % 53 dB
		spl = '53 dB SPL';
	elseif ind == 3
		area(data_RM.freqs,data_RM.rates(:,5),'EdgeColor', '#20116B',...
			'FaceColor', '#20116B','LineWidth',linewidth) % 73 dB
		spl = '73 dB SPL';
	end
	plot(data_RM.freqs([1 end]),[1 1]*data_RM.spont,'-','LineWidth',...
		linewidth, 'Color',spont_color)
	xline(CF, '--', 'Color', CF_color,'LineWidth',linewidth);
	text(0.05, 0.90, spl, 'Units', 'normalized', ...
		'VerticalAlignment', 'top', 'FontSize',legsize)
	grid on
	hold off
	ylim([0 max(data_RM.rates, [], 'all')+5])
	yticks([0 40])
	xlim(plot_range)
	xticks(x_label)
	xticklabels(x_label./1000)
	set(gca, 'XScale', 'log','fontsize',fontsize);
	if ind == 1
		xlabel('Tone Freq. (kHz)')
		hLegend = legend('', '', 'CF', 'Location','northeast ', ...
			'Box','off');
		hLegend.ItemTokenSize = [12,6];
	elseif ind == 2
		ylabel('Avg. Rate (sp/s)')
		xticklabels([])
	elseif ind == 3
		title('Response Map', 'FontSize', titlesize)
		xticklabels([])
	end
end


%% Plot MTF

% Analysis
data_MTF = analyzeMTF(params_MTF);

% Plot
h(4) = subplot(3, 4, 4);
hold on
line([1 data_MTF.fms(end)], [1 1]*data_MTF.rate(1),'Color',spont_color,...
	'LineWidth',linewidth);
errorbar(data_MTF.fms,data_MTF.rate, ...
	data_MTF.rate_std/sqrt(params_MTF.nrep),'.k', 'LineWidth',...
	linewidth, 'CapSize',4);
line(data_MTF.fms,data_MTF.rate,'Color','k', 'Marker','.', ...
	'MarkerSize',5, 'MarkerFaceColor','w', 'LineWidth', linewidth);
hold off
set(gca, 'XScale', 'log');
xlim([data_MTF.fms(1) data_MTF.fms(end)])
xticks([1 2 5 10 20 50 100 200 500])
xlabel('Modulation Freq. (Hz)')
ylim(MTF_range)
set(gca,'fontsize',fontsize)
grid on
title('MTF', 'FontSize', titlesize)
ylabel('Avg. Rate (sp/s)')
hLegend = legend('Unmodulated', 'Location','southwest', ...
	'Box','off');
hLegend.ItemTokenSize = [6,6];

%% Plot Binaural WB-TIN

% Analysis
data_WB = analyzeWBTIN(params_WB, CF);

% Plot
colors = {'#82BB95', '#3F985C', '#034E1C'};
num_WB =  length(data_WB);
for ind = 1:num_WB

	h(5+ind) = subplot(3,4, 4+ind);
	hold on
	line([data_WB{ind}.fpeaks(2) data_WB{ind}.fpeaks(end)], ...
		[1 1]*mean(data_WB{ind}.rate(1)),'Color',spont_color, ...
		'LineWidth', linewidth);
	errorbar(data_WB{ind}.fpeaks(:,1),data_WB{ind}.rate(:,1),...
		data_WB{ind}.rate_std(:,1)/(sqrt(params_WB{ind}.nrep)), ...
		'LineWidth', linewidth,'Color', colors{2}, 'CapSize',4);
	errorbar(data_WB{ind}.fpeaks(:,2),data_WB{ind}.rate(:,2),...
		data_WB{ind}.rate_std(:,2)/(sqrt(params_WB{ind}.nrep)), ...
		'LineWidth', linewidth,'Color', colors{3}, 'CapSize',4);
	xline(CF, '--', 'Color', CF_color, 'LineWidth', linewidth);
	hold off

	for ii = 1:length(params_WB{ind}.SNR)
		leg(ii,:) = [num2str(round(params_WB{ind}.SNR(ii),1)) 'dB SNR'];
	end

	xlim(plot_range);
	ylim(y_range)
	set(gca, 'XScale', 'log');
	xticks(x_label)
	xticklabels(x_label./1000)
	yticks(y_ticks)
	xlabel('Tone Freq. (kHz)')
	ylabel('Avg Rate (sp/s)')
	if ind == 1
		hLegend = legend(char('Noise Alone', leg), 'location','northwest',...
			'FontSize',legsize, 'Box','off');
		hLegend.ItemTokenSize = [6,6];
	end
	set(gca,'fontsize',fontsize)
	grid on
	title(['N_0 = ' num2str(params_WB{ind}.No) 'dB SPL'], 'FontSize',...
		titlesize)
end

%% Plot STRF

% Analysis
data_STRF = analyzeSTRF(params_STRF);

% Plot
h(5) = subplot(3,4, 8);
imagesc(data_STRF.t*1000, data_STRF.f./1000, ...
	data_STRF.H2ex_strf-data_STRF.H2in_strf, data_STRF.clims_strf);
set(gca,'Ydir','normal','XLim',data_STRF.tlims*1000, ...
	'YLim',[0 plot_range(2)]./1000)
colormap(redblue)
xlabel('Latency (ms)');
ylabel('Freq. (kHz)')
set(gca,'fontsize',fontsize)
title('STRF', 'FontSize', titlesize)
box off

%% Plot On-CF NB-TIN

% Analysis
data_NB = analyzeNBTIN(params_NB, CF);
num_NB = length(data_NB);
num_SNRs = length(params_NB{1}.SNR);
rate_onCF = zeros(num_NB, num_SNRs);
rate_std_onCF = zeros(num_NB, num_SNRs);
label = cell(num_NB, 1);
for ind = 1:num_NB

	NBTIN = data_NB{ind};
	param = params_NB{ind};
	
	% Find frequency in fpeaks nearest to CF
	[~,closestIndex] = min(abs(NBTIN.fpeaks(:, 1)-CF));

	% Create a 3x3 array by repeating the row 3 times
	SNRs = categorical(param.SNRNo);
	rate_onCF(ind, :) = NBTIN.rate(closestIndex,:);
	rate_std_onCF(ind,:) = NBTIN.rate_std(closestIndex,:);

	ttest_r = NBTIN.ttest_sig;
	if ttest_r < 0.05
		sig = '*';
	else
		sig = '';
	end

	label{ind} = sprintf('N_o = %d, p = %.2f%s', param.No, ttest_r, sig);
end

% Plot on-CF results
SNRs = repmat(SNRs, 3, 1);
h(9) = subplot(3,4, 9);
errorbar(SNRs', rate_onCF', (rate_std_onCF/sqrt(param.nrep))', ...
	'Linewidth', linewidth);
hLegend = legend(label, 'location', 'northeast', 'FontSize',legsize, ...
	'Box','off');
hLegend.ItemTokenSize = [6,6];
xlabel('SNR (dB SPL)')
ylabel('Avg. Rate (sp/s)')
set(gca, 'fontsize', fontsize)
grid on
title('NB-TIN','FontSize',titlesize)
ylim([0 120])
yticklabels(0:20:60)
box off

%% Plot On-CF WB-TIN

% Calculate legend with t-test results
for ind = 1:length(data_WB)
	WBTIN = data_WB{ind};
	param = params_WB{ind};
	if WBTIN.ttest_sig < 0.05
		sig = '*';
	else
		sig = '';
	end
	label{ind} = sprintf('N_o = %d, p = %.2f%s', param.No, ...
		WBTIN.ttest_sig, sig);
	WBTIN_rate_onCF(ind,:) = WBTIN.rate_onCF';
end

% Plot on-CF results
h(10) = subplot(3,4, 10);
errorbar(WBTIN.SNR_onCF', WBTIN_rate_onCF, ...
	(WBTIN.rate_std_onCF/sqrt(param.nrep))', 'Linewidth', linewidth);
hLegend = legend(label, 'location', 'northwest', 'FontSize',legsize, ...
	'Box','off');
hLegend.ItemTokenSize = [6,6];
xlabel('SNR (dB SPL)')
grid on
set(gca, 'fontsize', fontsize)
title('WB-TIN', 'FontSize',titlesize)
ylim([0 100])
yticklabels(0:20:60)
ylabel('Avg. Rate (sp/s)')
box off


%% Position All Plots

width = 0.17;
height = 0.32;
left = linspace(0.065, 0.815, 4);
bottom = [0.09, 0.62];

set(h(1), 'Position', [left(1),bottom(2),width,height/3]) % [bottom, left, width, height]
set(h(2), 'Position', [left(1),bottom(2)+height/3,width,height/3])
set(h(3), 'Position', [left(1),bottom(2)+height/3*2,width,height/3])

set(h(4), 'Position', [left(2),bottom(2),width,height])
set(h(5), 'Position', [left(3),bottom(2),width,height])
set(h(9), 'Position', [left(4),bottom(2),width,height])

set(h(6), 'Position', [left(1),bottom(1),width,height])
set(h(7), 'Position', [left(2),bottom(1),width,height])
set(h(8), 'Position', [left(3),bottom(1),width,height])
set(h(10), 'Position', [left(4),bottom(1),width,height])

%% Annotations

labels = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
labelleft= repmat(linspace(0.01, 0.76, 4), 1, 2);
labelbottom = [repmat(0.97,1, 4) repmat(0.44, 1, 4)];
for ii = 1:8
	annotation('textbox',[labelleft(ii) labelbottom(ii) 0.071 0.044],...
		'String',labels{ii},'FontWeight','bold','FontSize',labelsize,...
		'EdgeColor','none');
end

%% Export

if save_fig == 1
	saveFigure('Fig3_example_WBTIN')
end
end