function fig7_dog_analysis(save_fig)
% dog_analysis plots Figure 7.
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN, fitDifferenceofGaussians, createDoG
%	 .mat files required: Neural_Data
%	 spreadsheets required: DoG_Analysis.xlsx, Data_Table.xlsx
%
% Author: J. Fritzinger
% Created: 2025-01-10; Last revision: 2025-05-12
%
% -------------------------------------------------------------------------

%% Load in data

[~, datapath, ~, ppi] = get_paths();
analysisTable = readtable(fullfile(datapath, 'DoG_Analysis.xlsx'), ...
	'PreserveVariableNames',true);
sessions = readtable(fullfile(datapath, 'Data_Table.xlsx'),...
	'PreserveVariableNames',true);

%% Set up figure

figure('Position',[50,50,6.929*ppi,4*ppi]);
scattersize = 6;
legsize = 6;
titlesize = 8;
fontsize = 7;
linewidth = 1;
labelsize = 13;

%% Seperate data into groups

suprathreshold_ind = analysisTable.SNR==40; % 40 dB SNR
No23 = analysisTable.No==20 | analysisTable.No==23; % 3, 23, and 43 No
good_fit = analysisTable.R2>=0.4; % Good fit
good_unit = analysisTable.V_p>=0.4; % Good unit
binaural = analysisTable.Binmode==2; % Binmode
noise_bw = analysisTable.NoiseBW==4 | analysisTable.NoiseBW==3; % BW = 4 octaves

%% Example unit 

% Examples
putative_neuron = 'R29_TT4_P2_N10';

% Find example in spreadsheet
load(fullfile(datapath,'Neural_Data', [putative_neuron '.mat']), 'data');
s_ind = strcmp(sessions.Putative_Units, putative_neuron);
CF = sessions.CF(s_ind);

% Get data for each stimulus
type = 'Log';
params_WB = data(8,2); % Gets binaural WB-TIN stimuli, 43 dB SPL 

% Analysis
data_WB = analyzeWBTIN(params_WB, []);
data_WB = data_WB{1};
[DOGparams, ~, ~, ~] = fitDifferenceofGaussians(data_WB, 2);
CF_exc = DOGparams(end, 5); % DoG Parameters
CF_inh = DOGparams(end, 6); % DoG Parameters

% WBTIN
fpeaks = data_WB.fpeaks;
fpeaks = log10(fpeaks);
noise_alone = data_WB.rate(1, end); %mean(data_WB.rate(1, :)); %

% Noise-alone analysis
px = [fpeaks(2) fpeaks(end) fpeaks(end) fpeaks(2)];
rate_pos = mean(data_WB.rate(1,end))+mean(data_WB.rate_std(1,end))/...
	(sqrt(params_WB{1}.nrep));
rate_neg = mean(data_WB.rate(1,end))-mean(data_WB.rate_std(1,end))/...
	(sqrt(params_WB{1}.nrep));
py = [rate_neg rate_neg rate_pos rate_pos];

% Plot data
h(15) = subplot(3, 7, 15);
hold on
patch(px, py, [0.9 0.9 0.9], 'EdgeColor', [0.9 0.9 0.9])
line([fpeaks(2) fpeaks(end)], [1 1]* noise_alone,'Color',[0.4 0.4 0.4],...
	'LineWidth', linewidth); % noise alone
xline(log10(CF), '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth); % CF line
errorbar(fpeaks(:,end),data_WB.rate(:, end),data_WB.rate_std(:,end)/...
	(sqrt(params_WB{1}.nrep)),'.', 'LineWidth', linewidth,'Color', '#3F985C');

% Plot DoG
fpeaks_dog = logspace(log10(fpeaks(2,1)), log10(fpeaks(end,1)), 100);
[dog, gauss_exc, gauss_inh] = createDoG(DOGparams(end,:), fpeaks_dog);
plot(fpeaks_dog, gauss_exc+noise_alone, 'LineWidth',...
	linewidth, 'Color', 'r')
plot(fpeaks_dog, -gauss_inh+noise_alone, ...
	'LineWidth',linewidth, 'Color', 'b')
plot(fpeaks_dog, dog+noise_alone, 'Color', ...
	[0.8500 0.3250 0.0980], 'LineWidth',linewidth, 'color', 'k');

% Labels
xlabel('Tone Freq. (kHz)')
ylabel('Avg. Rate (sp/s)')
xlim([fpeaks(2, 1) fpeaks(end,1)])
set(gca, 'XScale', 'log');
grid on
xticks_lin = [1000, 2000, 3000, 4000, 5000, 7000];
xticks(log10(xticks_lin))
xticklabels(xticks_lin./1000)
set(gca,'FontSize',fontsize)
title('Example Fit', 'FontSize',titlesize)

%% Plot good R^2 fits 
 
% Analysis
data_good_ind = suprathreshold_ind & No23 & noise_bw & good_unit;
data_bad_ind = suprathreshold_ind & No23 & noise_bw & ~good_unit;
R2_values = analysisTable.R2(data_good_ind);
R2_values_bad = analysisTable.R2(data_bad_ind);

% Plot 
h(1) = subplot(3, 7, 1);
edges = linspace(0, 1, 21);
hold on
histogram(R2_values_bad, edges,'FaceColor', "#D95319")
histogram(R2_values, edges, 'FaceColor',"#0072BD");
xlabel('R^2')
set(gca, 'FontSize', fontsize)
title('Fit R^2 Values', 'fontsize', titlesize)
xticks(0:0.5:1)
xlim([0 1])
ylim([0 50])
ylabel('# Fits')
hLegend = legend('V_p < 0.4','V_p > 0.4', 'location', 'northwest',...
	'fontsize', legsize, 'box', 'off');
hLegend.ItemTokenSize = [6,6];

%% Strengths

toAnalyze = suprathreshold_ind & No23 & good_fit & binaural & noise_bw & good_unit;
ind = find(toAnalyze);
sig = analysisTable.g_ratio(ind);

% Histogram of g_ratio
h(2) = subplot(3, 7, 2);
edges = linspace(0, 1.6, 15);
histogram(sig, edges)
hold on
xline(mean(sig), 'r', 'LineWidth', linewidth)
xline(median(sig), ':r', 'LineWidth', linewidth)
xlim([0 1.7])
ylim([0 50])
ylabel('# Neurons')
xlabel('Str. Ratio (Inh/Exc)')
set(gca, 'FontSize', fontsize)
grid on
title('Strength Ratio', 'FontSize',titlesize)
hLegend = legend('', 'Mean', 'Median', 'fontsize', legsize, ...
	'location', 'northwest', 'box', 'off');
hLegend.ItemTokenSize = [10,14];

%% Plot Sigmas

No_names = {'3', '23', '43'};
sigma_exc = zeros(3, 2); sigma_inh = zeros(3,2);
for itype = 1:3
	iNo = 2;

	toAnalyze = suprathreshold_ind & No23 & good_fit & binaural & noise_bw & good_unit;
	ind = find(toAnalyze);
	CFs = analysisTable.CF(ind);
	if itype == 1
		sig = 10.^analysisTable.CF_exc(ind).*(10.^analysisTable.sigma_exc(ind)-1);
		sigma_exc(iNo,1:2) = [mean(sig./CFs) median(sig./CFs)];
	elseif itype == 2
		sig = 10.^analysisTable.CF_inh(ind).*(10.^analysisTable.sigma_inh(ind)-1);
		sigma_inh(iNo,1:2) = [mean(sig./CFs) median(sig./CFs)];
	else
		sig_exc = 10.^analysisTable.CF_exc(ind).*(10.^analysisTable.sigma_exc(ind)-1);
		sig_inh = 10.^analysisTable.CF_inh(ind).*(10.^analysisTable.sigma_inh(ind)-1);
		sig = sig_inh./sig_exc;
		fprintf('No=%s: %d/%d neurons above 1\n', No_names{iNo}, sum(sig>=1), length(sig))
	end

	% Fit linear regression line
	mdl = fitlm(log10(CFs), log10(sig));
	x = 0.3:0.5:10000;
	p(1) = mdl.Coefficients.Estimate(2,1);
	p(2) = mdl.Coefficients.Estimate(1,1);
	p(3) = mdl.Coefficients.pValue(2);
	if itype == 1
		mdlfit2 = 10^p(2) .* x.^p(1);
	end
	mdlfit = 10^p(2) .* x.^p(1);

	% Plot scatter and fit lines 
	h(2+itype) = subplot(3, 7, 2+itype);
	hold on
	if itype == 1
		if p(3) < 0.05
			plot(x/1000, mdlfit, 'r');
		else
			plot(x/1000, mdlfit, '--r');
		end
		plot([min(CFs) max(CFs)]/1000, [min(CFs) max(CFs)], 'k', 'linewidth',1)
		scatter(CFs/1000, sig,scattersize, 'filled', 'r', 'MarkerEdgeColor', 'k')
	elseif itype == 2
		plot(x/1000, mdlfit2, 'r');
		if p(3) < 0.05
			plot(x/1000, mdlfit, 'b');
		else
			plot(x/1000, mdlfit, '--b');
		end
		plot([min(CFs) max(CFs)]/1000, [min(CFs) max(CFs)], 'k', 'linewidth',1)
		scatter(CFs/1000, sig,scattersize, 'filled', 'b', 'MarkerEdgeColor', 'k')
	else
		scatter(CFs/1000, sig,scattersize, 'filled', 'MarkerEdgeColor', ...
			'k', 'MarkerFaceColor', '#009E73')
		yline(1, 'k')
	end

	% Figure labels
	axislabels = [100 200 500 1000 2000 5000 10000];
	xlim([0.3 11])
	set(gca, 'FontSize', fontsize, 'yscale', 'log', 'xscale', 'log')
	xticks(axislabels./1000)
	grid on
	if itype == 1
		title('Bandwidth (\sigma)', 'FontSize',titlesize)
		ylabel('\sigma_e_x_c (kHz)')
	elseif itype == 2
		ylabel('\sigma_i_n_h (kHz)')
	elseif itype == 3
		ylim([0.1 10])
		yticks([0.1 0.2 0.5 1 2 5])
		xlabel('CF (kHz)')
		ylabel('\sigma_i_n_h/\sigma_e_x_c')
	end

	if itype == 1 || itype == 2
		xticklabels([])
		ylim([40 15000])
		yticks(axislabels)
		yticklabels(axislabels./1000)
		hLegend = legend(sprintf('p = %.3f', p(3)),...
			'Location', 'NorthWest', 'fontsize', legsize, 'box', 'off');
		hLegend.ItemTokenSize = [6,6];
		msg = sprintf('y = 10^{%0.2f} \\times x^{%0.2f}', p(2), p(1));
		text(0.22, 0.15, msg, 'Units', 'normalized', ...
			'VerticalAlignment', 'top', 'FontSize',legsize+1)
	end
end


%% CF Scatter plots

for itype = 1:3
	toAnalyze = suprathreshold_ind & No23 & good_fit & binaural & noise_bw;
	ind = find(toAnalyze);
	CFs = analysisTable.CF(ind)/1000;

	if itype == 1
		CF_exc = 10.^analysisTable.CF_exc(ind)/1000;
	elseif itype == 2
		CF_inh = 10.^analysisTable.CF_inh(ind)/1000;
	else
		CF_exc = 10.^analysisTable.CF_exc(ind);
		CF_inh = 10.^analysisTable.CF_inh(ind);
		CF_offset = log2(CF_exc) - log2(CF_inh);
	end

	h(5+itype) = subplot(3, 7, 5+itype);
	hold on
	if itype == 1 || itype == 2
		plot([min(CFs) max(CFs)], [min(CFs) max(CFs)], 'k', 'LineWidth',1)
		scatter(CFs, CF_exc, scattersize, 'filled', 'r', ...
			'MarkerEdgeColor', 'k');
	elseif itype == 2
		plot([min(CFs) max(CFs)], [min(CFs) max(CFs)], 'k', 'LineWidth',1)
		scatter(CFs, CF_inh, scattersize, 'filled', 'b',...
			'MarkerEdgeColor', 'k');
	elseif itype == 3
		yline(0, 'k')
		scatter(CFs, CF_offset,scattersize, 'filled', ...
			'MarkerEdgeColor', 'k', 'MarkerFaceColor', '#009E73');
	end

	set(gca, 'FontSize', fontsize)
	xticks(axislabels./1000)
	set(gca, 'xscale', 'log')
	xlim([0.3 11])

	if itype == 1
		title('Center Freq. (f)', 'FontSize',titlesize)
		yticks(axislabels./1000)
		ylim([0.3 17])
		set(gca, 'yscale', 'log')
		ylabel('{\it f}_{exc} (kHz)')
		xticklabels([])
	elseif itype == 2
		yticks(axislabels./1000)
		ylim([0.3 17])
		set(gca, 'yscale', 'log')
		ylabel('{\it f}_{inh} (kHz)')
		xticklabels([])
	elseif itype == 3
		ylim([-2 2])
		xlabel('CF (kHz)')
		ylabel('{\it f}_{exc} - {\it f}_{inh} (kHz)')
	end
	grid on
end

% Can you convert to octaves w.r.t. the excitatory center frequency?
%differences = (CF_exc - CF_inh)./analysisTable.CF(ind);
%diff_wrt_CF = mean(differences);
%exc_CF = log2(CF_exc./analysisTable.CF(ind));
%inh_CF = log2(CF_inh./analysisTable.CF(ind));
%diff_wrt_CF = mean(exc_CF - inh_CF);


%% Scatter plot 
scattersize = 18;

is_interest = No23 & suprathreshold_ind & binaural & good_fit;
indices = find(is_interest);

colors_str = [253,204,138; 252,141,89; 215,48,31]./256;

sigma_ratio = log(analysisTable.sigma_inh(indices) ./ analysisTable.sigma_exc(indices));
g_ratio = log(analysisTable.g_ratio(indices));

CFs = analysisTable.CF(indices);
edges = [250 2000 4000 15000];
[~,~,bin] = histcounts(CFs, edges);
colors = zeros(length(CFs), 3);
for ii = 1:length(CFs)
	colors(ii,:) = colors_str(bin(ii),:);
end

h(9) = subplot(3, 7, 9);
scatter(sigma_ratio, g_ratio, scattersize, 'filled', 'MarkerEdgeColor','k')
hold on
xline(0)
yline(0)
xlabel('Log Bandwidth Ratio (\sigma_i_n_h/\sigma_e_x_c)')
ylabel('Log Strength Ratio (g_i_n_h/g_e_x_c)')
set(gca, 'fontsize', fontsize)
xlim([-2.5 4.5])
ylim([-2.5 2])
sigma_ex = zeros(5, 1);
g_ex = zeros(5, 1);
for iexample = 1:5
	switch iexample
		case 1 % Q1
			putative_neuron = 'R25_TT2_P8_N16'; % : R25 TT3 P9 N3
		case 2 % Q2
			putative_neuron = 'R24_TT2_P12_N2';
		case 3 % Q2
			putative_neuron = 'R25_TT3_P8_N3';
		case 4 % Q3
			putative_neuron = 'R27_TT3_P1_N5';
		case 5 % Q4
			putative_neuron = 'R27_TT2_P8_N1';
	end
	
	index = find(strcmp(putative_neuron, analysisTable.Putative) & No23 & ...
		binaural & suprathreshold_ind);
	sigma_ex(iexample) = log(analysisTable.sigma_inh(index) ./ analysisTable.sigma_exc(index));
	g_ex(iexample) = log(analysisTable.g_ratio(index));

end
scatter(sigma_ex, g_ex, scattersize,'filled', 'MarkerFaceColor', [253,141,60]./256, ...
	'MarkerEdgeColor','k')

% Count # in each quadrant 
quad(1) = sum(sigma_ratio>0 & g_ratio>0);
quad(2) = sum(sigma_ratio<0 & g_ratio>0);
quad(3) = sum(sigma_ratio<0 & g_ratio<0);
quad(4) = sum(sigma_ratio>0 & g_ratio<0);
for iquad = 1:4
	fprintf('Quadrant %d: n=%d\n',iquad, quad(1))
end 

%% Plot examples

for iexample = 1:5

	% Load in one example
	switch iexample
		case 1 % Q1
			putative_neuron = 'R25_TT2_P8_N16'; % : R25 TT3 P9 N3
		case 2 % Q2
			putative_neuron = 'R25_TT3_P8_N3';
		case 3 % Q2
			putative_neuron = 'R24_TT2_P12_N2';
		case 4 % Q3
			putative_neuron = 'R27_TT3_P1_N5';
		case 5 % Q4
			putative_neuron = 'R27_TT2_P8_N1';
	end
	
	% Find example in spreadsheet
	load(fullfile(datapath,'Neural_Data', [putative_neuron '.mat']), 'data');
	s_ind = strcmp(sessions.Putative_Units, putative_neuron);
	CF = sessions.CF(s_ind);

	% Get data for each stimulus
	type = 'Log';
	params_WB = data(7,2); % Gets binaural WB-TIN stimuli, 23 dB 
	data_WB = analyzeWBTIN(params_WB, []);
	data_WB = data_WB{1};
	[DOGparams, ~, ~, ~] = fitDifferenceofGaussians(data_WB, 2);

	% WBTIN
	fpeaks = data_WB.fpeaks;
	fpeaks = log10(fpeaks);
	rate = data_WB.rate;
	rate_std = data_WB.rate_std;
	param = params_WB{1};
	noise_alone = mean(rate(1,end));

	% Noise-alone analysis
	px = [fpeaks(2) fpeaks(end) fpeaks(end) fpeaks(2)];
	rate_pos = mean(rate(1,end))+mean(rate_std(1,end))/(sqrt(param.nrep));
	rate_neg = mean(rate(1,end))-mean(rate_std(1,end))/(sqrt(param.nrep));
	py = [rate_neg rate_neg rate_pos rate_pos];

	% Plot data
	h(9+iexample) = subplot(3, 7, 9+iexample);
	box on
	hold on
	patch(px, py, [0.9 0.9 0.9], 'EdgeColor', [0.9 0.9 0.9])
	line([fpeaks(2) fpeaks(end)], [1 1]* rate(1,end),'Color',...
		[0.4 0.4 0.4], 'LineWidth', linewidth); % noise alone
	xline(log10(CF), '--', 'Color', [0.4 0.4 0.4], 'LineWidth',linewidth);
	errorbar(fpeaks(:,end),rate(:, end),rate_std(:,end)/(sqrt(param.nrep)),'.', ...
		'LineWidth', linewidth,'Color', '#3F985C');

	% Plot DoG
	fpeaks_dog = logspace(log10(fpeaks(2,1)), log10(fpeaks(end,1)), 100);
	[dog, ~, ~] = createDoG(DOGparams(end,:), fpeaks_dog);
	plot(fpeaks_dog, dog+noise_alone, 'Color', ...
		[0.8500 0.3250 0.0980], 'LineWidth',linewidth, 'color', 'k');

	% Labels
	if ismember(iexample, [2, 4, 5])
		xlabel('Tone Freq. (kHz)')
	end
	if ismember(iexample, [3, 4])
		ylabel('Rate (sp/s)')
	else
		yticklabels([])
	end
	grid on
	xlim([fpeaks(2, 1) fpeaks(end,1)])
	set(gca, 'XScale', 'log');
	xticks_lin = [1000, 2000, 3000, 5000, 7000];
	xticks(log10(xticks_lin))
	xticklabels(xticks_lin./1000)
	box off
	set(gca, 'fontsize', fontsize)
end


%% Annotations / Rearrange figure 

left = [0.0600    0.2500    0.4500    0.6600];
bottom = linspace(0.1, 0.66, 3);
width = 0.13;
height = 0.27;

set(h(15), 'Position', [left(1) 0.66+height-0.18 0.11 0.18]);
set(h(1), 'Position', [left(1) 0.42 0.11 0.18]) % [left, bottom, width, height]
set(h(2), 'Position', [left(1) bottom(1) 0.11 0.18])
set(h(8), 'Position', [left(2) bottom(1) width height])
set(h(7), 'Position', [left(2) bottom(2) width height])
set(h(6), 'Position', [left(2) bottom(3) width height])
set(h(5), 'Position', [left(3) bottom(1) width height])
set(h(4), 'Position', [left(3) bottom(2) width height])
set(h(3), 'Position', [left(3) bottom(3) width height])
set(h(9), 'Position', [left(4) 0.3 0.32 0.4])

width = 0.09;
height = 0.12;
set(h(10), 'Position', [0.89 0.80 width height])
set(h(11), 'Position', [0.775 0.80 width height])
set(h(12), 'Position', [left(4) 0.80 width height])
set(h(13), 'Position', [left(4) 0.08 width height])
set(h(14), 'Position', [0.88 0.08 width height])

% Create textbox
annotation('textbox',[left(1)-0.06 0.98 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1)-0.06 0.63 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(1)-0.06 0.32 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2)-0.04 0.98 0.0826 0.0385],'String',{'D'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(3)-0.04 0.98 0.0826 0.0385],'String',{'E'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(4)-0.04 0.98 0.0826 0.0385],'String',{'F'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

% Create arrow
arrow_locs = [0.7334, 0.7294, 0.7569, 0.6145;...
	0.801 0.7595, 0.731 0.559;...
	0.9639 0.8476 0.7534 0.5833;...
	0.6613 0.7640 0.2118 0.4109;...
	0.9680 0.8336 0.2001 0.4548];
for iarrow = 1:5
	annotation('arrow',arrow_locs(iarrow, 1:2),arrow_locs(iarrow, 3:4),...
		'LineWidth',1, 'HeadLength',5, 'HeadWidth',5);
end

%% Export Figure

if save_fig == 1
	saveFigures('Fig7_dog_analysis')
end
end