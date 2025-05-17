function supp1_data_distribution(save_fig)
% supp_data_distribution plots Figure S1.
%	This function loads in the putative neurons spreadsheet and plots the MTF
%	distribution, BMF distribution, WMF distribution, hybrid BMF/WMF
%	distribution, CF distribution for all neurons, and predictable variance
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN
%	 .mat files required: FigS1_Vp.mat
%	 spreadsheets required: Data_Table.xlsx
%
% Author: J. Fritzinger
% Created: 2024-05-29; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------


%% Load in spreadsheet

[~, datapath, ~, ppi] = get_paths();
sessions = readtable(fullfile(datapath, 'Data_Table.xlsx'), ...
	'PreserveVariableNames',true);
load(fullfile(datapath, 'FigS1_Vp.mat'), 'Vp_combined')

%% Set up figure

figure('Position',[100,100,4*ppi,6*ppi])
fontsize = 7;
titlesize = 9;
labelsize = 12;

%% Get MTFs for each putative neuron

MTFs = sessions.MTF;
num_sesh = length(MTFs);
MTF_type = zeros(num_sesh,1);
for isesh = 1:num_sesh
	MTF_shape = MTFs{isesh};
	if contains(MTF_shape, 'H')
		MTF_type(isesh) = 3;
	elseif strcmp(MTF_shape, 'BE')
		MTF_type(isesh) = 1;
	elseif strcmp(MTF_shape, 'BS')
		MTF_type(isesh) = 2;
	else % Flat
		MTF_type(isesh) = 4;
	end
end
MTF_names = categorical({'BE','BS','Hybrid','Flat'});
MTF_names = reordercats(MTF_names,{'BE','BS','Hybrid','Flat'});
num_BE = sum(MTF_type==1);
num_BS = sum(MTF_type==2);
num_H = sum(MTF_type==3);
num_n = sum(MTF_type==4);
num_types = [num_BE num_BS num_H num_n];

% Plot
h(1) = subplot(4, 2, 1);
bar(MTF_names,num_types, 'black');
set(gca, 'FontSize', fontsize)
title('MTF Type', 'FontSize',titlesize)
grid on
ylim([0 110])
ylabel('Number of Neurons')


%% BMFs

% Get BMFs/WMFs
BE_MTFs = strcmp(sessions.MTF, 'BE');
BMFs = sessions.BMF(BE_MTFs);
BMFs(isnan(BMFs)) = [];

edges = [0.2 2 4 8 16 32 64 128 254 512 1028];
edges2 = zeros(10, 1);
for iedge = 1:10
	edges2(iedge) = sqrt(edges(iedge)*edges(iedge+1));
end

% BMFs
h(2) = subplot(4, 2, 3);
histogram(BMFs, edges2,'FaceColor', '#0072BD', 'EdgeColor','k')
hold on
xline(exp(median(log(BMFs(BMFs~=0)))), 'k', 'LineWidth',1.5)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'FontSize', fontsize)
title('BMF & WMF', 'fontsize', titlesize)
set(gca, 'XScale', 'log');

str = ['n=' num2str(length(BMFs))];
annotation('textbox',[0.14 0.61 0.096 0.053],'String',str,...
	'FontSize',10,'EdgeColor','none');

%% WMFs

BS_MTFs = strcmp(sessions.MTF, 'BS');
WMFs = sessions.WMF(BS_MTFs);
WMFs(isnan(WMFs)) = [];

h(3) = subplot(4, 2, 5);
histogram(WMFs, edges2,'FaceColor', '#D95319', 'EdgeColor','k')
hold on
ylabel('Number of Neurons')
xline(exp(median(log(WMFs(WMFs~=0)))), 'k', 'LineWidth',1.5)
set(gca, 'FontSize', fontsize)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'XScale', 'log');


str = ['n=' num2str(length(WMFs))];
annotation('textbox',[0.14 0.41 0.096 0.053],'String',str,...
	'FontSize',10,'EdgeColor','none');

%% Hybrids

H_MTFs = contains(sessions.MTF, 'H');
WMFs = sessions.WMF(H_MTFs);
WMFs(isnan(WMFs)) = [];
BMFs = sessions.BMF(H_MTFs);
BMFs(isnan(BMFs)) = [];

% Plot 
h(4) = subplot(4, 2, 7);
histogram(BMFs, edges2)
hold on
histogram(WMFs, edges2)
xlabel('BMF or WMF (Hz)')
hLegend = legend('BMF', 'WMF', 'Location','west');
hLegend.ItemTokenSize = [6,6];
set(gca, 'FontSize', fontsize)
xticks([2 4 8 16 32 64 128 254 512])
set(gca, 'XScale', 'log');

str = ['n=' num2str(length(BMFs))];
annotation('textbox',[0.14 0.21 0.096 0.053],'String',str,...
	'FontSize',10,'EdgeColor','none');

%% Get CFs for each putative neuron

% WBTIN Diotic
CFs = sessions.CF;
edges = [0 500 1000 2000 4000 8000 13000];
names = categorical({'<0.5', '0.5-1', '1-2', '2-4', '4-8', '8+'});
names = reordercats(names,{'<0.5', '0.5-1', '1-2', '2-4', '4-8', '8+'});
CF = CFs;
CF(CF==0) = [];
[N, ~] = histcounts(CF, edges);

% Plot
h(5) = subplot(4, 2, 2);
bar(names,N,'FaceColor', 'k', 'EdgeColor','k');
grid on
ylabel('Number of Neurons')
xlabel('CF (kHz)')
set(gca, 'FontSize', fontsize)
title('CF Distribution', 'fontsize', titlesize)
ylim([0 110])

%% Predictable Variance

% Plot
for iSNR = 1:3
	Vp = Vp_combined(:,iSNR);

	good_Vp = Vp>0.4;
	bad_Vp = Vp<=0.4;
	Vp_good = Vp(good_Vp);
	Vp_bad = Vp(bad_Vp);

	% Plot V_p
	h(5+iSNR) = subplot(4, 2, 2+2*iSNR); % 4, 6, 8
	edges = linspace(0, 1, 21);
	histogram(Vp_good, edges);
	hold on
	histogram(Vp_bad, edges)
	set(gca, 'FontSize', fontsize)
	xticks(0:0.2:1)
	if iSNR == 1
		xticklabels([])
		title('V_p for WB-TIN Data', 'fontsize', titlesize)
	elseif iSNR == 2
		xticklabels([])
		ylabel('Number of Fits')
	else
		xlabel('Proportion of V_p')
	end
	xlim([0 1])
	grid on
	hold on
	ylim([0 52])
	n_good = length(Vp_good);
	n_noisy = length(Vp_bad);
	leg = {['Good, n=' num2str(n_good)], ['Noisy, n=' num2str(n_noisy)]};
	hLegend = legend(leg, 'location', 'northwest');
	hLegend.ItemTokenSize = [6,6];
end

%% Rearrange Tiles 

dimensions = [0.334 0.2]; % width, height
left = [repmat(0.127, 1, 4) repmat(0.622, 1, 4)]; % 0.127 first column
bottom = repmat([0.765 0.465 0.265 0.065], 1, 2);
for ii = 1:8
	set(h(ii), 'Position', [left(ii) bottom(ii) dimensions])
end

labels = {'BE', 'BS', 'Hybrid'};
bottom = [0.559 0.344 0.129];
for ii = 1:3
	annotation('textbox',[0.07 bottom(ii) 0.104 0.043],'String',labels{ii},...
		'EdgeColor','none', 'Rotation',90);
end

labels = {'40 dB SNR', '30 dB SNR', '20 dB SNR'};
bottom = [0.1 0.3 0.5];
for ii = 1:3
	annotation('textbox',[0.55 bottom(ii) 0.3 0.043],'String',labels{ii},...
		'EdgeColor','none', 'Rotation',90);
end

labels = {'A', 'B', 'C', 'D'};
left = repmat([0.046 0.546], 1, 2);
bottom = [repmat(0.957, 1, 2) repmat(0.654, 1, 2)];
for ii = 1:4
	annotation('textbox',[left(ii) bottom(ii) 0.096 0.053],'String',labels{ii},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
end

%% Export figure

if save_fig == 1
	saveFigure('Supp1_data_distribution')
end
end