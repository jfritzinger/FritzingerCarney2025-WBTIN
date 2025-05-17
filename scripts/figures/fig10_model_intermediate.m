function fig10_model_intermediate(save_fig)
% model_intermediate plots Figure 10.
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required:
%	 .mat files required:
%	 spreadsheets required: 
%
% Author: J. Fritzinger
% Created: 2025-01-21; Last revision: 2025-03-07
%
% -------------------------------------------------------------------------


%% Analysis and Plotting
colors = {'#1b9e77', '#d95f02', '#7570b3'}; 

[~, datapath, ~, ppi] = get_paths();
figure('Position',[50,50,4.567*ppi,2.4*ppi]);
legsize = 6;
titlesize = 8;
fontsize = 7;
linewidth = 1;
labelsize = 13;

h = gobjects(2,2);
index = [2, 5; 3, 6];
for iexample = 1:2

	filename = ['Model_Component_Example_' num2str(iexample) '.mat'];
	load(fullfile(datapath, filename), 'param', 'IC')

	[~, model_MTF(:,1), model_std(:,1), ~, smoothed(:,1)] = plotMTF(...
		param, IC{iexample}.avBS_lo, 0);
	[~, model_MTF(:,2), model_std(:,2), ~, smoothed(:,2)] = plotMTF(...
		param, IC{iexample}.avBS_on, 0);
	[~, model_MTF(:,3), model_std(:,3), ~, smoothed(:,3)] = plotMTF(...
		param, IC{iexample}.avBS_hi, 0);
	[~, model_MTF(:,4), model_std(:,4), ~, smoothed(:,4)] = plotMTF(...
		param, IC{iexample}.avIC, 0);
	CS_params = IC{iexample}.CS_params;

	% Plot MTF
	names = {'Example BE', 'Example BS'};
	h(iexample, 1) = subplot(2, 3, index(iexample, 1));
	hold on
	yline(model_MTF(1,2),'Color',colors{2}, 'linewidth', linewidth);
	param.all_fms(1) = 1.2;
	plot(param.all_fms, model_MTF(:,2), 'linewidth', linewidth, 'Color',colors{2})
	hold off
	xtick = [1 2 5 20 50 100 200 500];
	xlim([1 530])
	xticklabels([])
	ylabel('Avg. Rate (sp/s)')
	set(gca,'XTick',xtick,'XScale', 'log')
	
	grid on
	axis_set = axis;
	if iexample == 1
		axis_set(3) = 0;
		axis_set(4) = 30;
	else
		axis_set(3) = 0;
		axis_set(4) = 50;
	end

	modifiers = CS_params(1:2);
	hold on
	yline(model_MTF(1,1).*modifiers(1),'Color',colors{1}, 'linewidth', linewidth);
	plot(param.all_fms, model_MTF(:,1).*modifiers(1), ...
		'linewidth', linewidth, 'Color',colors{1})
	yline(model_MTF(1,3).*modifiers(2),'Color',colors{3}, 'linewidth', linewidth);
	plot(param.all_fms, model_MTF(:,3)*modifiers(2), ...
		'linewidth', linewidth, 'Color',colors{3})

	axis(axis_set);
	set(gca, 'FontSize', fontsize)
	title(names{iexample}, 'FontSize',titlesize)
	if iexample == 2
		hleg = legend({'', 'On-CF', '','Low-CF * S_l_o_w','', 'High-CF * S_h_i_g_h'},...
			'Location','best', 'FontSize',legsize, 'box', 'off');
		hleg.ItemTokenSize = [8,8];
	end

	h(iexample, 2) = subplot(2, 3, index(iexample, 2));
	hold on
	yline(model_MTF(1,4),'Color',[0.4 0.4 0.4], 'linewidth', linewidth);
	errorbar(param.all_fms, model_MTF(:,4), model_std(:,4),'k', 'linewidth', linewidth)
	hold off
	xtick = [1 2 5 20 50 100 200 500];
	xlim([1 530])
	xlabel('Mod. Freq. (Hz)')
	ylabel('Avg. Rate (sp/s)')
	set(gca,'XTick',xtick,'XScale', 'log')
	
	grid on
	set(gca, 'FontSize', fontsize)
	if iexample == 2
		ylim([0 50])
		hleg = legend('Unmodulated', 'Location','best', 'FontSize',...
			legsize, 'box', 'off');
		hleg.ItemTokenSize = [8,8];
	else
		ylim([0 30])
	end
end

% Load in images
img = imread(fullfile(datapath, 'BroadInhModel_Schematic.png'));

p(1) = subplot(2, 3, 1);
imshow(img(1:1100,2:end,:))

p(2) = subplot(2, 3, 4);
imshow(img(1100:end,2:end,:))

%% Arrange 
left = [0.04 0.42 0.75];
bottom = [0.57 0.145];
width = 0.22;
height = 0.31;

for iexample = 1:2
	for ii = 1:2
		set(h(iexample,ii), 'Position', [left(iexample+1) bottom(ii) width height])
	end
end

set(p(1), 'Position',[left(1) bottom(1)-0.02 0.27 0.35])
set(p(2), 'Position',[left(1) bottom(2)-0.02 0.27 0.35])

annotation('textbox',[left(1)-0.03 0.96 0.0826 0.0385],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(2)-0.085 0.96 0.0826 0.0385],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[left(3)-0.085 0.96 0.0826 0.0385],'String',{'C'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Export

if save_fig == 1
	saveFigure('Fig10_model_intermediate')
end
end