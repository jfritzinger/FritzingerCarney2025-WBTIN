function supp5_comparison_RM(save_fig)
% supp_comparison_RM plots Figure S2.
%	This function loads in the putative neurons spreadsheet RM comparison
%	spreadsheet and plots two example neurons (V and I type, multiple
%	levels) and box plots of variance explained for different RM types. 
%
% Inputs:
%    save_fig - 1=save .png and .eps files, 0=don't save
%
% Requirements:
%	 m-files required: getPathsWBTIN, analyzeRM, analyzeWBTIN
%	 .mat files required: -
%	 spreadsheets required: Data_Table.xlsx, RM_Comparison.xlsx
%
% Author: J. Fritzinger
% Created: 2024-05-29; Last revision: 2025-03-04
%
% -------------------------------------------------------------------------


%% Load in data

[~, datapath, ~, ppi] = get_paths();
sessions = readtable(fullfile(datapath, 'Data_Table.xlsx'), ...
	'PreserveVariableNames',true);
analysisTable = readtable(fullfile(datapath, 'RM_Comparison.xlsx'), ...
	'PreserveVariableNames',true);

%% Figure Parameters

figure('Position',[100,100,5*ppi,6*ppi]);
legsize = 7;
titlesize = 9;
fontsize = 8;
labelsize = 13;
linewidth = 1;
scattersize = 15;

%% Analysis 

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
	RM_type = sessions.RM_Type{indx};
	CF = sessions.CF(indx);

	% Get data for each stimulus
	params_RM = data{2, 2};
	params_WB = data(6:8,2); % Gets binaural WB-TIN stimuli

	% General analysis
	data_RM = analyzeRM(params_RM);
	datas_WB = analyzeWBTIN(params_WB, []);

	% Find peak rates
	max_RM = max(max(data_RM.rates));
	all_WB = cellfun(@(x) max(max(x.rate)), datas_WB, 'UniformOutput', false);
	max_WB = max([all_WB{:}]);
	max_rate = max([max_WB max_RM]);

	RM_colors = {'k', 'k', '#A49BD0', '#5E50A9', '#20116B'};
	WB_colors = {'#82BB95', '#3F985C', '#034E1C'};

	for iNo = 1:3

		param_WB = params_WB{iNo};
		data_WB = datas_WB{iNo};

		% Find BF nearest to WBTIN overall level
		lbw = param_WB.fpeak_mid * 2^(-param_WB.bandwidth/2);
		hbw = param_WB.fpeak_mid * 2^(param_WB.bandwidth/2);
		overall_spl = param_WB.No + 10*log10(hbw-lbw);
		[~,iBF] = min(abs(data_RM.SPL-overall_spl));
		BF_spl = data_RM.SPL(iBF);


		freqs_interp = logspace(log10(data_WB.fpeaks(2,1)), ...
			log10(data_WB.fpeaks(end,1)), 50);
		rate_WBTIN = data_WB.rate(2:end,end);
		rate_interp_WBTIN = interp1(data_WB.fpeaks(2:end,1),rate_WBTIN,...
			freqs_interp,'pchip'); 

		rate_RM = data_RM.rates(:,iBF);
		freqs_mid = logspace(log10(data_RM.freqs(1)), ...
			log10(data_RM.freqs(end)), 700);
		rate_mid_RM = interp1(data_RM.freqs,rate_RM,freqs_mid,'pchip'); 
		[~, starting] = min(abs(freqs_mid-freqs_interp(1)));
		[~, ending] = min(abs(freqs_mid-freqs_interp(end)));
		rate_interp_RM = interp1(freqs_mid(starting:ending),...
			rate_mid_RM(starting:ending),freqs_interp,'pchip');

		R_int = corrcoef(rate_interp_RM,rate_interp_WBTIN);
		r2 = R_int(1, 2).^2;

		%% Plot

		% Plot RM/WBTIN
		h(iNo+(ineuron-1)*3) = subplot(3, 3, iNo+(ineuron-1)*3);
		hold on
		line([data_WB.fpeaks(2) data_WB.fpeaks(end)], ...
			[1 1]*data_WB.rate(1,end),'Color','#034E1C', 'linewidth', linewidth);
		line([data_RM.freqs(1) data_RM.freqs(end)], ...
			[1 1]*data_RM.spont, 'Color', '#20116B', 'linewidth', linewidth);
		area(data_RM.freqs, data_RM.rates(:,iBF), 'EdgeColor',RM_colors{iBF}, ...
			'FaceColor',RM_colors{iBF}, 'FaceAlpha',0.5);
		errorbar(data_WB.fpeaks(:,1),data_WB.rate(:,end),...
			data_WB.rate_std(:,end)/(sqrt(param_WB.nrep)), ...
			'Color',WB_colors{iNo}, 'linewidth', linewidth);
		xline(CF, '--', 'Color', [0.3 0.3 0.3], 'LineWidth',linewidth)

		set(gca,'FontSize',fontsize)
		grid on
		set(gca, 'XScale', 'log');
		xlim([data_WB.fpeaks(2) data_WB.fpeaks(end)])
		ylim([0 max_rate+10])

		xL=xlim;
		yL=ylim;
		message = ['R^2 = ' num2str(round(r2, 2))];
		text(1.05*xL(1),0.95*yL(2),message,'HorizontalAlignment','left',...
			'VerticalAlignment','top', 'FontSize',legsize)

		if iNo == 3
			title(sprintf('RM Correlation, %s-type', RM_type), ...
				'FontSize',titlesize)
			xticklabels([])
			if ineuron == 2
				hLegend = legend('', '', ['RM, ' num2str(BF_spl) ' dB SPL'], ...
					['WB-TIN, ' num2str(round(param_WB.No_overall)) ' dB SPL'],...
					'', 'fontsize', legsize);
			else
				hLegend = legend('', '', [num2str(BF_spl) ' dB SPL'], ...
					[num2str(round(param_WB.No_overall)) ' dB SPL'], ...
					'fontsize', legsize);
			end
			hLegend.ItemTokenSize = [6,6];
		elseif iNo == 2
			ylabel('Spike rate (sp/s)')
			xticklabels([])
			hLegend = legend('', '', [num2str(BF_spl) ' dB SPL'], ...
				[num2str(round(param_WB.No_overall)) ' dB SPL'],...
				'fontsize', legsize);
			hLegend.ItemTokenSize = [6,6];
		else
			xlabel('Tone Frequency (kHz)')
			xticklabels(xticks/1000)
			hLegend = legend('', '', [num2str(BF_spl) ' dB SPL'], ...
				[num2str(round(param_WB.No_overall)) ' dB SPL'],...
				'fontsize', legsize);
			hLegend.ItemTokenSize = [6,6];
		end
	end
end

%% Scatter plots 

bin = analysisTable.binmode==2;
%contra = analysisTable.binmode==1;
SNR_40 = analysisTable.SNR==40;
No(:,1) = analysisTable.WB_No==3;
No(:,2) = analysisTable.WB_No==23;
No(:,3) = analysisTable.WB_No==43;

RM_types = unique(analysisTable.RM);
RM_type = RM_types([1, 5,8]);
num_type = length(RM_type);

Nos = {'3', '23', '43'};
types = {'I-type', 'Onset','V-type'};
for iNo = 1:3

	type_ind = strcmp(analysisTable.RM, 'V');
	ind = bin & SNR_40 & No(:,iNo) & type_ind;
	num = sum(ind);

	% Seperate into types
	R2 = NaN(num,3);
	for itype = 1:num_type

		type_ind = strcmp(analysisTable.RM, RM_type{itype});
		ind = bin & SNR_40 & No(:,iNo) & type_ind;
		R2(1:sum(ind),itype) = analysisTable.R2(ind);

		meanR2 = round(mean(R2(:,itype), 'omitnan'), 2);
		medianR2 = round(median(R2(:,itype), 'omitnan'), 2);
		num_units =  length(R2(~isnan(R2(:,itype)), :));
		fprintf('No = %s: %s, mean=%g, median=%g, n=%d \n', Nos{iNo},...
			types{itype}, meanR2, medianR2, num_units);

	end
	groups = [ones(num, 1), 2*ones(num, 1), 3*ones(num, 1)];

	h(iNo+6) = subplot(3, 3, iNo+6);
	jitter = 0.5; % Adjust the jitter amount as desired
	hold on
	boxplot(R2,'Colors','k', 'Widths',0.5);
	outliers = isoutlier(R2, 'quartiles');
	R2(outliers==1) = NaN;
	scatter(groups + (rand(size(R2)) - 0.5) * jitter, R2, scattersize, ...
		'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', '#0072BD');
	ylim([0 1])
	xlim([0.5 3.5])
	xticks([1 2 3])
	xticklabels(RM_type)
	if iNo == 1
		ylabel('Variance Explained (R^2)')
	elseif iNo == 2
		xlabel('Response Map Type')
		yticklabels([])
	else
		yticklabels([])
	end
	title([Nos{iNo} 'dB SPL'], 'fontsize', titlesize)
	set(gca, 'FontSize', fontsize)
	box off
end


%% Annotations 

col1 = 0.1;
col2 = 0.56;
width = 0.37;
height = 0.165;

% Example 1
set(h(1), 'Position', [col1 0.45 width height]) 
set(h(2), 'Position', [col1 0.615 width height])
set(h(3), 'Position', [col1 0.78 width height])

% Example 2
set(h(4), 'Position', [col2 0.45 width height])
set(h(5), 'Position', [col2 0.615 width height])
set(h(6), 'Position', [col2 0.78 width height])

% Scatter Plots
set(h(7), 'Position', [0.1 0.07 0.28 0.27])
set(h(8), 'Position', [0.38 0.07 0.28 0.27])
set(h(9), 'Position', [0.66 0.07 0.28 0.27])

% Create textbox
annotation('textbox',[0.0486 0.333 0.0769 0.0545],'String','C',...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.516 0.934 0.0769 0.0545],'String',{'B'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');
annotation('textbox',[0.0532 0.934 0.0769 0.0545],'String',{'A'},...
	'FontWeight','bold','FontSize',labelsize,'EdgeColor','none');

%% Export Figure

if save_fig == 1
	saveFigure('Supp5_comparison_RM')
end
end


